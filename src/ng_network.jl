using NetworkPDE, ModelingToolkit, DomainSets, Zygote, CSV, OrdinaryDiffEq,
    SymbolicIndexingInterface

include("./parse_csv.jl")
#@register_symbolic my_floor(T::Type, x)::UInt64
sym_int_floor(x) = floor(Int, x);
@register_symbolic sym_int_floor(x)::Int
global N = 1;
node_tstep(t) = sym_int_floor(t/bc_dt)+1
@register_symbolic node_tstep(t)::Int
abstract type Pipe <: AbstractNetworkComponent end
abstract type FluxNode <: AbstractNetworkComponent end

const Wave_Speed = 454.55;
pressure_from_density(density) = Wave_Speed^2 * density
density_from_pressure(pressure) = pressure / Wave_Speed^2


function create_ode_problem!(net, dx, tspan; bc_dt=3600)
    @parameters t
    Dt = Differential(t)

    #CFL ~ min(dx)/dt < sqrt(P'(ρ))
    # dt = 0.9*minimum([get_metadata(net.edges[i], :dx) for i in 1:length(net.edges)]) /
    #     sqrt(gradient(pressure_from_density, 1.0)[1])
    dt = dx/gradient(pressure_from_density, 1.0)[1];
    add_metadata(net, :dt, dt);

    function dρ(ρ, ϕ, dx, dt)
        return @. -(1 / dx) * (ϕ[2:end] - ϕ[1:(end - 1)])
    end

    function dϕ(ρ, ϕ, dx, dt, β, pressure_from_density)
        return -(1/dx).*(pressure_from_density(ρ[2:end]) .-
            pressure_from_density(ρ[1:(end - 1)])) .-
            β.*(ϕ[2:(end - 1)] .* abs.(ϕ[2:(end - 1)]))./ ((ρ[2:end] .+ ρ[1:end-1])/2)
    end
    function dϕ_bdy(ρ, ϕ, dx, dt, β, pressure_from_density)
        return -(1/dx)*(pressure_from_density(ρ[2]) -
            pressure_from_density(ρ[1])) -
            (β*(ϕ * abs.(ϕ))/((ρ[1] + ρ[2])/2))
    end

    function instantiate(::Type{Pipe}, comp::AbstractNetworkComponent; name)
        L = get_metadata(comp, :L)
        dx = get_metadata(comp, :dx)
        β = get_metadata(comp, :λ) / (2 * get_metadata(comp, :D))
        nx = floor(Int, L / dx)
        add_metadata(comp, :ρ_xs, dx/2:dx:L-dx/2)
        add_metadata(comp, :ϕ_xs, 0:dx:L)
        @variables ρ(t)[1:nx] ϕ(t)[1:(nx + 1)] #staggered grid
        add_metadata(comp, :ρ_xs, dx/2:dx:L-dx/2)
        add_metadata(comp, :ϕ_xs, 0:dx:L)

        drho = dρ(ρ, ϕ, dx, dt)
        internal_eqs = vcat(Symbolics.scalarize.([Dt(ρ) ~ drho;
                                                  Dt(ϕ[2:(end - 1)]) ~ dϕ(ρ, ϕ, dx, dt, β,
                                                                          pressure_from_density)])...)
        ics = [ρ .=> get_metadata(comp, :initial_density).(get_metadata(comp, :ρ_xs)),
               ϕ .=> get_metadata(comp, :initial_mass_flux).(get_metadata(comp, :ϕ_xs))]
        bcs = [] #fill in later
        return ODESystem(internal_eqs, t, name = name, defaults = Dict(vcat(Symbolics.scalarize.(ics)...)))
    end

    function instantiate(::Type{FluxNode}, comp::AbstractNetworkComponent; name)
#        @parameters q[1:get_metadata(comp, :num_bc_tsteps)]
        @variables ρ(t)
        ics = [ρ => get_metadata(comp, :initial_node_density)]#,
#               q .=> get_metadata(comp, :initial_node_mass_flow)]
        return ODESystem(Equation[], t, Num[ρ], Num[], name = name, defaults=Dict(ics))
    end

    function connect(net::Network; name)
        coupling_eqs = Equation[]
        @parameters q[1:get_metadata(net, :num_bc_tsteps)*num_vertices(net)]
        for v in net.vertices
            inc_edges = [net.edges[i] for i in get_metadata(v, :incident_edges)]
            num_inc_edges = length(inc_edges)
            sgns = [get_metadata(e, :to_node) == v.id ? 1 : -1 for e in inc_edges]
            cross_sections = [π * (get_metadata(e, :D) / 2)^2 for e in inc_edges]
            edge_dxs = [get_metadata(e, :dx) for e in inc_edges]
            edge_ρs = []
            edge_ϕs = []
            node_ρ = unknowns(v.metadata[:model], v.metadata[:model].var_to_name[:ρ])
            #node_q = unknowns(v.metadata[:model], v.metadata[:model].var_to_name[:q])
            for i in 1:num_inc_edges
                model = inc_edges[i].metadata[:model]
                if (sgns[i] == 1)
                    push!(edge_ρs,
                          unknowns(
                              model, model.var_to_name[:ρ][length(model.var_to_name[:ρ])]))
                    tmp = [unknowns(model, model.var_to_name[:ϕ][i]) for i in length(model.var_to_name[:ϕ]):-1:length(model.var_to_name[:ϕ])-2]
                    push!(edge_ϕs, tmp)
                else
                    push!(edge_ρs, unknowns(model, model.var_to_name[:ρ][1]))
                    tmp = [unknowns(model, model.var_to_name[:ϕ][i]) for i in 1:3]
                    push!(edge_ϕs, tmp)
                end
            end
            # define equation terms for readability
            t_node_flow = q[get_metadata(net, :get_node_bc)(v.id, t)]
            t_pipe_flux_approx = [edge_ϕs[i][1] + sgns[i]*(-3*edge_ϕs[i][1]+4*edge_ϕs[i][2]-edge_ϕs[i][3])/(dx/2) for i in 1:num_inc_edges]
            t_pipe_flow = sum([sgns[i]*cross_sections[i]*t_pipe_flux_approx[i] for i in 1:num_inc_edges])
            t_node_vol = sum([(edge_dxs[i]/2)*cross_sections[i] for i in 1:num_inc_edges])
            dρ = (t_node_flow + t_pipe_flow)/t_node_vol;
            push!(coupling_eqs,
                  Dt(node_ρ) ~ dρ)
            for i in 1:num_inc_edges
                tmp_ρ = sgns[i] == 1 ? [edge_ρs[i]; node_ρ] : [node_ρ; edge_ρs[i]];
                push!(coupling_eqs,
                    Dt(edge_ϕs[i][1]) ~ dϕ_bdy(tmp_ρ, edge_ϕs[i][1], dx, dt,
                        get_metadata(net.edges[i], :λ) / (2 * get_metadata(net.edges[i], :D)),
                        pressure_from_density))
            end
        end
        vars = Symbolics.wrap.([unknowns(v.metadata[:model], v.metadata[:model].var_to_name[:ρ]) for v in net.vertices]);
        parent_system = ODESystem(coupling_eqs, t, vars, [q...], name = name); 
        subsystems = vcat([get_metadata(e, :model) for e in net.edges]..., [get_metadata(v, :model) for v in net.vertices]...)
        return compose(parent_system, subsystems)
    end

    function instantiate(network::Network; name)
        connection_eqs = []
        for (i, e) in enumerate(network.edges)
            add_metadata(e, :model,
                         instantiate(get_metadata(e, :type),
                                     e, name = Symbol("edge_$(i)")))
        end
        for (i, v) in enumerate(network.vertices)
#            add_metadata(v, :num_bc_tsteps, get_metadata(network, :num_bc_tsteps))
            add_metadata(v, :model,
                         instantiate(get_metadata(v, :type),
                                     v, name = Symbol("node_$(i)")))
        end
        return connect(network, name = name)
    end
    num_bc_tsteps = Base.floor(Int, (tspan[end]-tspan[1] - 0.1)/bc_dt) + 1
    get_node_bc(id, t) = num_bc_tsteps*(id-1) + sym_int_floor(t/bc_dt)+1
    add_metadata(net, :get_node_bc, get_node_bc)
    add_metadata(net, :num_bc_tsteps, num_bc_tsteps)
    sys = instantiate(net, name=:sys)
    sys = structural_simplify(sys)
#    new_sys = ODESystem(equations(sys), independent_variable(sys), unknowns(sys), parameters(sys)[good_param_inds], defaults=default_values(sys), name=:new_sys)
#    new_sys = structural_simplify(new_sys)

    prob = ODEProblem(sys, nothing, tspan, strict_mode=false)
    return sys, prob
end

