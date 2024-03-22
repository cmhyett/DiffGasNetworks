using NetworkPDE, ModelingToolkit, DomainSets, Zygote, CSV, OrdinaryDiffEq, SymbolicIndexingInterface
using Test, SafeTestsets, Plots

include("../src/parse_csv.jl")
include("../src/ng_network.jl")

@testset "no dynamics" begin
    pressure_from_density(density) = Wave_Speed^2 * density;
    density_from_pressure(pressure) = pressure / Wave_Speed^2;
    dx = 1_000; #m
    basepath = "./network_specs/single_pipe_test/";
    net = parse_network_from_files(dx,
                                   basepath * "network_data/pipes.csv",
                                   basepath * "network_data/nodes.csv",
                                   basepath * "initial_conditions/pipe_ic.csv",
                                   basepath * "initial_conditions/node_ic.csv",
                                   basepath * "network_data/params.csv",
                                   density_from_pressure,
                                   pressure_from_density);
    T = 1*3600; #s
    tspan = (0.0, T); #s
    sys, prob = create_ode_problem!(net, dx, tspan);
    dt = get_metadata(net, :dt) #appropriate dt is set when parsing
    params = [sys.node_1.q => 0.0,
              sys.node_2.q => 0.0]
    prob = remake(prob, p=params);
    sol = solve(prob, Euler(), dt=dt, saveat=prob.tspan[1]:dt:prob.tspan[2]); #use Euler() for now, generalize later
    node1_ρ_inds = variable_index(sys, sys.node_1.ρ);
    node2_ρ_inds = variable_index(sys, sys.node_2.ρ);
    edge1_ρ_inds = variable_index(sys, sys.edge_1.ρ);
    @test sol[1][node1_ρ_inds] == sol[end][node1_ρ_inds]
    @test sol[1][node2_ρ_inds] == sol[end][node2_ρ_inds]
    @test sol[1][edge1_ρ_inds] == sol[end][edge1_ρ_inds]
    mass_t = [sum(sol[i][edge1_ρ_inds]) for i in 1:length(sol)]
    @test all(x->x≈mass_t[end], mass_t)

    sol = solve(prob, RK4(), dt=dt, saveat=prob.tspan[1]:dt:prob.tspan[2]); 
    @test sol[1][node1_ρ_inds] == sol[end][node1_ρ_inds]
    @test sol[1][node2_ρ_inds] == sol[end][node2_ρ_inds]
    @test sol[1][edge1_ρ_inds] == sol[end][edge1_ρ_inds]
    mass_t = [sum(sol[i][edge1_ρ_inds]) for i in 1:length(sol)]
    @test all(x->x≈mass_t[end], mass_t)

end

@testset "settling initial condition" begin
    pressure_from_density(density) = Wave_Speed^2 * density;
    density_from_pressure(pressure) = pressure / Wave_Speed^2;
    dx = 1_000; #m
    basepath = "./network_specs/single_pipe_test/";
    net = parse_network_from_files(dx,
                                   basepath * "network_data/pipes.csv",
                                   basepath * "network_data/nodes.csv",
                                   basepath * "initial_conditions/pipe_ic.csv",
                                   basepath * "initial_conditions/node_ic.csv",
                                   basepath * "network_data/params.csv",
                                   density_from_pressure,
                                   pressure_from_density);
    T = 3*3600; #s
    tspan = (0.0, T); #s
    sys, prob = create_ode_problem!(net, dx, tspan);
    dt = get_metadata(net, :dt) #appropriate dt is set when parsing
    params = [sys.node_1.q => 0.0,
              sys.node_2.q => 0.0]

    #set initial pipe density
    node1_ρ_inds = variable_index(sys, sys.node_1.ρ);
    node2_ρ_inds = variable_index(sys, sys.node_2.ρ);
    edge1_ρ_inds = [variable_index(sys, sys.edge_1.ρ[i]) for i in 1:length(sys.edge_1.ρ)];
    L = get_metadata(net.edges[1], :L);
    initial_density = density_from_pressure(7000000) .+ density_from_pressure(500000)*[-tanh(x/L-1/2) for x in get_metadata(net.edges[1], :ρ_xs)]

    u0 = [(i in edge1_ρ_inds) ? initial_density[findfirst(x->x==i, edge1_ρ_inds)] : prob.u0[i] for i in 1:length(prob.u0)]
    u0[node1_ρ_inds] = initial_density[1]
    u0[node2_ρ_inds] = initial_density[end]
    prob = remake(prob, p=params, u0 = u0);
    sol = solve(prob, Euler(), dt=dt, saveat=prob.tspan[1]:dt:prob.tspan[2]); #use Euler() for now, generalize later
    mass_t = [sum(sol[i][edge1_ρ_inds])+sol[i][node1_ρ_inds]+sol[i][node2_ρ_inds] for i in 1:length(sol)]
    @test all(x->isapprox(x, mass_t[end], rtol=1e-4), mass_t)
    @test all(isapprox.(sol[end][edge1_ρ_inds], sum(sol[1][edge1_ρ_inds])/length(edge1_ρ_inds), atol=dx/L))

    sol = solve(prob, RK4(), dt=dt, saveat=prob.tspan[1]:dt:prob.tspan[2]);
    mass_t = [sum(sol[i][edge1_ρ_inds])+sol[i][node1_ρ_inds]+sol[i][node2_ρ_inds] for i in 1:length(sol)]
    @test all(x->isapprox(x, mass_t[end], rtol=1e-4), mass_t)
    @test all(isapprox.(sol[end][edge1_ρ_inds], sum(sol[1][edge1_ρ_inds])/length(edge1_ρ_inds), atol=dx/L))

end

@safetestset "3 pipe" begin
    using NetworkPDE, ModelingToolkit, DomainSets, Zygote, CSV, OrdinaryDiffEq, SymbolicIndexingInterface

    include("../src/parse_csv.jl")
    include("../src/ng_network.jl")

    pressure_from_density(density) = Wave_Speed^2 * density;
    density_from_pressure(pressure) = pressure / Wave_Speed^2;
    dx = 1_000; #m
    basepath = "./network_specs/three_pipe_test/";
    net = parse_network_from_files(dx,
                                   basepath * "network_data/pipes.csv",
                                   basepath * "network_data/nodes.csv",
                                   basepath * "initial_conditions/pipe_ic.csv",
                                   basepath * "initial_conditions/node_ic.csv",
                                   basepath * "network_data/params.csv",
                                   density_from_pressure,
                                   pressure_from_density);

    T = 3*3600; #s
    tspan = (0.0, T); #s
    global sys, prob = create_ode_problem!(net, dx, tspan); #'global' gets around main scope of 'eval' later
    dt = get_metadata(net, :dt) #appropriate dt is set when parsing
    params = [sys.node_1.q => 0.0,
              sys.node_2.q => 0.0,
              sys.node_3.q => 0.0]
    prob = remake(prob, p=params);
    edge_ρ_inds = vcat([[variable_index(sys, eval(Meta.parse("sys.edge_$(e_num).ρ[$(i)]"))) for i in 1:length(sys.edge_1.ρ)] for e_num in 1:length(net.edges)]...);
    node_ρ_inds = vcat([variable_index(sys, eval(Meta.parse("sys.node_$(n_num).ρ"))) for n_num in 1:length(net.vertices)]...);
    # no dynamics
    sol = solve(prob, RK4(), dt=dt);
    mass_t = [sum(sol[i][edge_ρ_inds])+sum(sol[i][node_ρ_inds]) for i in 1:length(sol)]
    @test all(x->isapprox(x, mass_t[end], rtol=1e-4), mass_t)

    # settling dynamics
    string_to_var(string) = eval(Meta.parse(string))
    function set_initial_condition(_u0, sys, var, vals)
        var_len = length(var)
        @assert var_len == length(vals)
        inds = [variable_index(sys, var[i]) for i in 1:var_len];
        new_u0 = [(i in inds) ? vals[findfirst(x->x==i, inds)] : _u0[i] for i in 1:length(_u0)]
        return new_u0;
    end

    global p1_ic = [density_from_pressure(7e6)*(1-x/get_metadata(net.edges[1], :L)) + density_from_pressure(7.5e6)*(x/get_metadata(net.edges[1], :L)) for x in get_metadata(net.edges[1], :ρ_xs)]
    global p2_ic = [density_from_pressure(7.5e6)*(1-x/get_metadata(net.edges[2], :L)) + density_from_pressure(8.0e6)*(x/get_metadata(net.edges[2], :L)) for x in get_metadata(net.edges[2], :ρ_xs)]
    global p3_ic = [density_from_pressure(8e6)*(1-x/get_metadata(net.edges[3], :L)) + density_from_pressure(7e6)*(x/get_metadata(net.edges[3], :L)) for x in get_metadata(net.edges[3], :ρ_xs)]
    global n1_ic = p1_ic[1]
    global n2_ic = p2_ic[1]
    global n3_ic = p3_ic[1]
    u0 = prob.u0;
    for i in 1:3
        u0 = set_initial_condition(u0, sys, string_to_var("sys.edge_$(i).ρ"), string_to_var("p$(i)_ic"))
    end
    for i in 1:3
        u0 = set_initial_condition(u0, sys, string_to_var("sys.node_$(i).ρ"), string_to_var("n$(i)_ic"))
    end
    prob = remake(prob, u0=u0);
    sol = solve(prob, RK4(), dt=dt)
    mass_t = [sum(sol[i][edge_ρ_inds])+sum(sol[i][node_ρ_inds]) for i in 1:length(sol)]
    @test all(x->isapprox(x, mass_t[end], rtol=1e-4), mass_t)

    #settling dynamics + matched withdrawals
    params = [sys.node_1.q => -50.0,
              sys.node_2.q => 25.0,
              sys.node_3.q => 25.0]
    prob = remake(prob, u0=u0, p=params);
    sol = solve(prob, RK4(), dt=dt)
    mass_t = [sum(sol[i][edge_ρ_inds])+sum(sol[i][node_ρ_inds]) for i in 1:length(sol)]
    @test all(x->isapprox(x, mass_t[end], rtol=1e-4), mass_t)
end

@safetestset "numerical accuracy" begin
    using NetworkPDE, ModelingToolkit, DomainSets, Zygote, CSV, OrdinaryDiffEq, SymbolicIndexingInterface, LinearAlgebra
    using DataInterpolations, Optim
    
    include("../src/parse_csv.jl")
    include("../src/ng_network.jl")

    pressure_from_density(density) = Wave_Speed^2 * density;
    density_from_pressure(pressure) = pressure / Wave_Speed^2;
    results = []
    dxs = reverse([4_000*2.0^(i) for i in 0:-1:-7])
    small_dt = 0.0;
    for dx in dxs
        println("dx = $(dx)")
        basepath = "./network_specs/three_pipe_test/";
        net = parse_network_from_files(dx,
                                       basepath * "network_data/pipes.csv",
                                       basepath * "network_data/nodes.csv",
                                       basepath * "initial_conditions/pipe_ic.csv",
                                       basepath * "initial_conditions/node_ic.csv",
                                       basepath * "network_data/params.csv",
                                       density_from_pressure,
                                       pressure_from_density);

        T = 3*3600; #s
        tspan = (0.0, T); #s
        global sys, prob = create_ode_problem!(net, dx, tspan); #'global' gets around main scope of 'eval' later
        dt = get_metadata(net, :dt) #appropriate dt is set when parsing
        small_dt = 0.0
        if (dx == dxs[1])
            small_dt = dt/10;
        end
        params = [sys.node_1.q => 0.0,
                  sys.node_2.q => 0.0,
                  sys.node_3.q => 0.0]
        prob = remake(prob, p=params);
        edge_ρ_inds = [[variable_index(sys, eval(Meta.parse("sys.edge_$(e_num).ρ[$(i)]"))) for i in 1:length(sys.edge_1.ρ)] for e_num in 1:length(net.edges)]
        edge_ϕ_inds = [[variable_index(sys, eval(Meta.parse("sys.edge_$(e_num).ϕ[$(i)]"))) for i in 1:length(sys.edge_1.ϕ)] for e_num in 1:length(net.edges)]
        node_ρ_inds = [variable_index(sys, eval(Meta.parse("sys.node_$(n_num).ρ"))) for n_num in 1:length(net.vertices)];
        
        #settling dynamics + matched withdrawals
        string_to_var(string) = eval(Meta.parse(string))
        function set_initial_condition(_u0, sys, var, vals)
            var_len = length(var)
            @assert var_len == length(vals)
            inds = [variable_index(sys, var[i]) for i in 1:var_len];
            new_u0 = [(i in inds) ? vals[findfirst(x->x==i, inds)] : _u0[i] for i in 1:length(_u0)]
            return new_u0;
        end

        global p1_ic = [density_from_pressure(7e6)*(1-x/get_metadata(net.edges[1], :L)) + density_from_pressure(7.5e6)*(x/get_metadata(net.edges[1], :L)) for x in get_metadata(net.edges[1], :ρ_xs)]
        global p2_ic = [density_from_pressure(7.5e6)*(1-x/get_metadata(net.edges[2], :L)) + density_from_pressure(8.0e6)*(x/get_metadata(net.edges[2], :L)) for x in get_metadata(net.edges[2], :ρ_xs)]
        global p3_ic = [density_from_pressure(8e6)*(1-x/get_metadata(net.edges[3], :L)) + density_from_pressure(7e6)*(x/get_metadata(net.edges[3], :L)) for x in get_metadata(net.edges[3], :ρ_xs)]
        global n1_ic = p1_ic[1]
        global n2_ic = p2_ic[1]
        global n3_ic = p3_ic[1]
        u0 = prob.u0;
        for i in 1:3
            u0 = set_initial_condition(u0, sys, string_to_var("sys.edge_$(i).ρ"), string_to_var("p$(i)_ic"))
        end
        for i in 1:3
            u0 = set_initial_condition(u0, sys, string_to_var("sys.node_$(i).ρ"), string_to_var("n$(i)_ic"))
        end

        params = [sys.node_1.q => -500.0,
                  sys.node_2.q => 250.0,
                  sys.node_3.q => 250.0]
        prob = remake(prob, u0=u0, p=params, tspan=(0.0, 100.0));
        @time push!(results, Dict(:sol => solve(prob, SSPRK432(), adaptive=true, dt=dt, save_end=true),
                                  :sys => sys,
                                  :net => net,
                                  :edge_ρ_inds => edge_ρ_inds,
                                  :edge_ϕ_inds => edge_ϕ_inds,
                                  :edge_ϕ_xs => [get_metadata(net.edges[i], :ϕ_xs) for i in 1:length(net.edges)],
                                  :edge_ρ_xs => [get_metadata(net.edges[i], :ρ_xs) for i in 1:length(net.edges)],
                                  :node_ρ_inds => node_ρ_inds,
                                  :prob => prob));
    end

    function create_interpolations(results; t_step=small_dt, edge_num=1)
        interpolations = []
        for i in 1:length(results)
            #create linear interpolations
            xs = vcat([0.0], results[i][:edge_ρ_xs][1], [results[i][:net].edges[1].metadata[:L]])
            pipe_ρ_inds = results[i][:edge_ρ_inds][1]
            node_ρ_inds = results[i][:node_ρ_inds];
            ys = vcat([results[i][:sol](t_step)[node_ρ_inds[1]]], results[i][:sol](t_step)[pipe_ρ_inds], [results[i][:sol](t_step)[node_ρ_inds[2]]])
            push!(interpolations, LinearInterpolation(ys, xs))
        end
        return interpolations
    end

    function calc_errors(results; t_step=small_dt, p=2)
        interpolations = create_interpolations(results; t_step=small_dt, edge_num=1)
        highest_refinement_xs = results[1][:edge_ρ_xs][1];
        highest_refinement_int = interpolations[1]
        errors = [norm(interpolations[i].(highest_refinement_xs) .- highest_refinement_int.(highest_refinement_xs), p)
                  for i in 1:length(results)];
        return errors
    end

    l1_errors = calc_errors(results, p=1);
    l2_errors = calc_errors(results, p=2);
    model(x, ps) = ps[1].*x .+ ps[2]
    guess = [2.0, 0.0]
    A = Curvefit(log.(l1_errors[2:end]), log.(dxs[2:end]), model, guess, LBFGS())
    @test isapprox(A.pmin[1], 2.0, atol=0.1)
    B = Curvefit(log.(l2_errors[2:end]), log.(dxs[2:end]), model, guess, LBFGS())
    @test isapprox(B.pmin[1], 1.5, atol=0.1)
    # plot(dxs[2:end], errors[2:end]/errors[end], xaxis=:log, yaxis=:log)
    # plot!(dxs[2:end], [x^2 for x in dxs[2:end]]/dxs[end]^2)
    # plot!(dxs[2:end], [x^2 for x in dxs[2:end]]/dxs[end]^2)
    # scatter!(dxs[2:end], errors[2:end]/errors[end], xaxis=:log, yaxis=:log)
    #plot(dxs[2:end], errors, xaxis=:log, yaxis=:log)
end

@safetestset "computational complexity" begin
    using NetworkPDE, ModelingToolkit, DomainSets, Zygote, CSV, OrdinaryDiffEq, SymbolicIndexingInterface, LinearAlgebra
    using DataInterpolations, Optim
    
    include("../src/parse_csv.jl")
    include("../src/ng_network.jl")

    pressure_from_density(density) = Wave_Speed^2 * density;
    density_from_pressure(pressure) = pressure / Wave_Speed^2;
    results = []
    dxs = reverse([4_000*2.0^(i) for i in 0:-1:-7])
    small_dt = 0.0;
    for dx in dxs
        println("dx = $(dx)")
        basepath = "./network_specs/three_pipe_test/";
        net = parse_network_from_files(dx,
                                       basepath * "network_data/pipes.csv",
                                       basepath * "network_data/nodes.csv",
                                       basepath * "initial_conditions/pipe_ic.csv",
                                       basepath * "initial_conditions/node_ic.csv",
                                       basepath * "network_data/params.csv",
                                       density_from_pressure,
                                       pressure_from_density);

        T = 3*3600; #s
        tspan = (0.0, T); #s
        global sys, prob = create_ode_problem!(net, dx, tspan); #'global' gets around main scope of 'eval' later
        dt = get_metadata(net, :dt) #appropriate dt is set when parsing
        small_dt = 0.0
        if (dx == dxs[1])
            small_dt = dt/10;
        end
        params = [sys.node_1.q => 0.0,
                  sys.node_2.q => 0.0,
                  sys.node_3.q => 0.0]
        prob = remake(prob, p=params);
        edge_ρ_inds = [[variable_index(sys, eval(Meta.parse("sys.edge_$(e_num).ρ[$(i)]"))) for i in 1:length(sys.edge_1.ρ)] for e_num in 1:length(net.edges)]
        edge_ϕ_inds = [[variable_index(sys, eval(Meta.parse("sys.edge_$(e_num).ϕ[$(i)]"))) for i in 1:length(sys.edge_1.ϕ)] for e_num in 1:length(net.edges)]
        node_ρ_inds = [variable_index(sys, eval(Meta.parse("sys.node_$(n_num).ρ"))) for n_num in 1:length(net.vertices)];
        
        #settling dynamics + matched withdrawals
        string_to_var(string) = eval(Meta.parse(string))
        function set_initial_condition(_u0, sys, var, vals)
            var_len = length(var)
            @assert var_len == length(vals)
            inds = [variable_index(sys, var[i]) for i in 1:var_len];
            new_u0 = [(i in inds) ? vals[findfirst(x->x==i, inds)] : _u0[i] for i in 1:length(_u0)]
            return new_u0;
        end

        global p1_ic = [density_from_pressure(7e6)*(1-x/get_metadata(net.edges[1], :L)) + density_from_pressure(7.5e6)*(x/get_metadata(net.edges[1], :L)) for x in get_metadata(net.edges[1], :ρ_xs)]
        global p2_ic = [density_from_pressure(7.5e6)*(1-x/get_metadata(net.edges[2], :L)) + density_from_pressure(8.0e6)*(x/get_metadata(net.edges[2], :L)) for x in get_metadata(net.edges[2], :ρ_xs)]
        global p3_ic = [density_from_pressure(8e6)*(1-x/get_metadata(net.edges[3], :L)) + density_from_pressure(7e6)*(x/get_metadata(net.edges[3], :L)) for x in get_metadata(net.edges[3], :ρ_xs)]
        global n1_ic = p1_ic[1]
        global n2_ic = p2_ic[1]
        global n3_ic = p3_ic[1]
        u0 = prob.u0;
        for i in 1:3
            u0 = set_initial_condition(u0, sys, string_to_var("sys.edge_$(i).ρ"), string_to_var("p$(i)_ic"))
        end
        for i in 1:3
            u0 = set_initial_condition(u0, sys, string_to_var("sys.node_$(i).ρ"), string_to_var("n$(i)_ic"))
        end

        params = [sys.node_1.q => -500.0,
                  sys.node_2.q => 250.0,
                  sys.node_3.q => 250.0]
        prob = remake(prob, u0=u0, p=params, tspan=(0.0, 100.0));
        @time push!(results, Dict(:sol => solve(prob, SSPRK432(), adaptive=true, dt=dt, save_end=true),
                                  :sys => sys,
                                  :net => net,
                                  :edge_ρ_inds => edge_ρ_inds,
                                  :edge_ϕ_inds => edge_ϕ_inds,
                                  :edge_ϕ_xs => [get_metadata(net.edges[i], :ϕ_xs) for i in 1:length(net.edges)],
                                  :edge_ρ_xs => [get_metadata(net.edges[i], :ρ_xs) for i in 1:length(net.edges)],
                                  :node_ρ_inds => node_ρ_inds,
                                  :prob => prob));
    end

end

