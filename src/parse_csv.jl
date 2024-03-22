function parse_pipes(pipes_path::String,
                     pipes_ic_path::String,
                     dx,
                     density_from_pressure;
                     T=Float32)::Array{Edge,1}
    pipes_parse_dict = Dict("length"=>T,
                            "diameter"=>T,
                            "friction_factor"=>T);
    ic_parse_dict = Dict("initial_pipe_flow"=>T,
                         "initial_pipe_pressure_in"=>T,
                         "initial_pipe_pressure_out"=>T);
    
    pipes_file = CSV.File(pipes_path, comment="#", types=pipes_parse_dict);
    ic_file = CSV.File(pipes_ic_path, comment="#", types=ic_parse_dict);

    num_pipes = length(pipes_file);
    output = [];
    for i in 1:num_pipes
        L = pipes_file[i].length;
        rho_out = density_from_pressure.(ic_file[i].initial_pipe_pressure_out);
        rho_in = density_from_pressure.(ic_file[i].initial_pipe_pressure_in);
        A = π*(pipes_file[i].diameter^2)/4;
        initial_density(x) = (1-x/L)*rho_in + (x/L)*rho_out;
        initial_mass_flux(x) = ic_file[i].initial_pipe_flow/A
        param_dict = Dict(:L => L,
                          :from_node => pipes_file[i].start_node,
                          :to_node => pipes_file[i].end_node,
                          :D => pipes_file[i].diameter,
                          :λ => pipes_file[i].friction_factor,
                          :initial_density => initial_density,
                          :initial_mass_flux => initial_mass_flux,
                          :dx => dx,
                          :type => Pipe)
        push!(output, Edge(i, param_dict));
    end
    return convert(Array{Edge, 1}, output);
end

function parse_nodes(node_filepath::String,
                     node_ic_filepath::String,
                     pipes::Array{Edge, 1},
                     density_from_pressure;
                     T=Float32)::Array{Vertex,1}
    nodes_parse_dict = Dict("height"=>T,
                            "min_pressure"=>T,
                            "max_pressure"=>T,
                            "min_flow"=>T,
                            "max_flow"=>T);
    nodes_ic_parse_dict = Dict("initial_nodal_flow"=>T,
                               "initial_nodal_pressure"=>T);
    node_file = CSV.File(node_filepath, comment="#", types=nodes_parse_dict, validate=false);
    node_ic_file = CSV.File(node_ic_filepath, comment="#", types=nodes_ic_parse_dict, validate=false);

    output = [];
    for i in 1:length(node_file)
        id = node_file[i].number;
        is_slack = (node_file[i].type == 1);

        incoming_edges = [];
        outgoing_edges = [];
        for i in 1:length(pipes)
            if (get_metadata(pipes[i], :to_node) == id)
                push!(incoming_edges, pipes[i].id);
            elseif (get_metadata(pipes[i], :from_node) == id)
                push!(outgoing_edges, pipes[i].id);
            end
        end
        param_dict = Dict(:incident_edges => [incoming_edges; outgoing_edges],
                          :initial_node_density => density_from_pressure(node_ic_file[i].initial_nodal_pressure),
                          :initial_node_mass_flow =>  node_ic_file[i].initial_nodal_flow,
                          :min_pressure => node_file[i].min_pressure,
                          :max_pressure => node_file[i].max_pressure,
                          :min_flow => node_file[i].min_flow,
                          :max_flow => node_file[i].max_flow,
                          :type => FluxNode,
                          :longitude => node_file[i].lon,
                          :latitude => node_file[i].lat);
        push!(output, Vertex(id, param_dict));
    end
    return convert(Array{Vertex, 1}, output);
end

function parse_network_params(params_path::String; T=Float32)
    params_parse_dict = Dict("temperature"=>T,
                             "gas_specific_gravity"=>T,
                             "specific_heat_capacity_ratio"=>T,
                             "initial_time"=>T,
                             "final_time"=>T,
                             "discretization_time_step"=>T,
                             "courant_number"=>T,
                             "output_dt"=>T,
                             "output_dx"=>T,
                             "save_final_state"=>Bool);
    file = CSV.File(params_path, comment="#", types=params_parse_dict)[1];
    return Dict{Symbol, Any}(:temperature=>file.temperature,
                :gas_specific_gravity=>file.gas_specific_gravity,
                :specific_heat_capacity_ratio=>file.specific_heat_capacity_ratio,
                :initial_time=>file.initial_time,
                :final_time=>file.final_time,
                :discretization_time_step=>file.discretization_time_step,
                :courant_number=>file.courant_number,
                :output_dt=>file.output_dt,
                :output_dx=>file.output_dx,
                :save_final_state=>file.save_final_state);
end

function parse_network_from_files(dx,
                                  pipe_path,
                                  node_path,
                                  pipe_ic_path,
                                  node_ic_path,
                                  param_path,
                                  density_from_pressure,
                                  pressure_from_density;
                                  T=Float32)
    network_params = parse_network_params(param_path, T=T);
    push!(network_params, :density_from_pressure => density_from_pressure);
    push!(network_params, :pressure_from_density => pressure_from_density);
    pipes = parse_pipes(pipe_path,
                        pipe_ic_path,
                        dx,
                        density_from_pressure, T=T);
    nodes = parse_nodes(node_path,
                        node_ic_path,
                        pipes,
                        density_from_pressure,
                        T=T);
#    make_pipe_ic_consistent!(pipes, nodes)
    return Network(pipes, nodes, network_params);
end

# function make_pipe_ic_consistent!(pipes, nodes)
#     for i in 1:length(pipes)
#         p = pipes[i];
#         start_node = nodes[findfirst(n->n.id==p.from_node, nodes)]
#         end_node = nodes[findfirst(n->n.id==p.to_node, nodes)]

#         # fix density
#         in_density = start_node.init_density;
#         out_density = end_node.init_density;
#         L = length(p.init_density)-1
#         p.init_density .= [(1-j/L)*in_density + (j/L)*out_density for j in 0:L];

#         # fix mass flux
#         cross_section = (π*(p.diameter/2)^2);
#         in_flux = -start_node.init_mass_flow/cross_section;
#         out_flux = end_node.init_mass_flow/cross_section;
#         L = length(p.init_mass_flux)-1
#         p.init_mass_flux .= [(1-j/L)*in_flux + (j/L)*out_flux for j in 0:L];
#     end
# end


function get_params(nw, flux_node_bc_path, tspan::Tuple{T, T}) where {T}

    file = CSV.File(flux_node_bc_path, comment="#", types=T);

    (length(file[1]) != length(nw.nodes)+1) &&
        @error "CSV boundary value file has $(length(file[1])-1) entries, expecting num_nodes = $(length(nw.nodes))"

    file_dt = file.time[2] - file.time[1];
    num_tsteps = ceil(Int, tspan[end]/file_dt);
    
    ps = zeros(T, length(nw.nodes)*num_tsteps);
    for i in 1:length(nw.nodes)
        for j in 1:num_tsteps
            ps[(i-1)*num_tsteps + j] = file[j][i+1];
        end
    end

    return ps;
end
