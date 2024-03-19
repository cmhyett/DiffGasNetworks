include("ng_network.jl")

dx = 1_000; #m
T = 15 * 60; # time [s]

basepath = "./network_specs/single_pipe_test/";

net = parse_network_from_files(dx,
    basepath * "network_data/pipes.csv",
    basepath * "network_data/nodes.csv",
    basepath * "initial_conditions/pipe_ic.csv",
    basepath * "initial_conditions/node_ic.csv",
    basepath * "network_data/params.csv",
    density_from_pressure,
    pressure_from_density);

tspan = (0.0, T);
bc_dt=100
sys, prob = create_ode_problem!(net, dx, tspan, bc_dt=bc_dt);

println("finished")
