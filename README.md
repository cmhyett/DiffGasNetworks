# NG Modeling setup

Clone this repository, then run
```
git submodule init
git submodule update
```
to clone the submodule `NetworkPDE.jl` dependency. Then startup julia in the `SecondaryFuel` directory and run
```
]activate .; instantiate
]dev ./NetworkPDE.jl
```
to point your julia install to the local `NetworkPDE.jl`

## Creating Problem
```
dx = 1_000; #m
tspan = (0.0, 15*60); #s
basepath = "./network_specs/israel_reduced/";
net = parse_network_from_files(dx,
                             basepath * "network_data/pipes.csv",
                             basepath * "network_data/nodes.csv",
                             basepath * "initial_conditions/pipe_ic.csv",
                             basepath * "initial_conditions/node_ic.csv",
                             basepath * "network_data/params.csv",
                             density_from_pressure,
                             pressure_from_density);
sys, prob = create_ode_problem(net, dx, tspan);

#sys isa ODESystem - needed for namespacing
#net isa Network - needed for metadata/structure
#prob isa ODEProblem - needed to solve/remake problem
```

## Solving the problem
```
dt = get_metadata(net, :dt) #appropriate dt is set when parsing
params = [sys.node_1.q => 0.0,
      sys.node_2.q => 0.0,
      sys.node_3.q => 0.0,
      sys.node_4.q => 0.0,
      sys.node_5.q => 0.0,
      sys.node_6.q => 0.0,
      sys.node_7.q => 0.0,
      sys.node_8.q => 0.0,
      sys.node_9.q => 0.0,
      sys.node_10.q => 0.0,
      sys.node_11.q => 0.0,
      sys.node_12.q => 0.0,
      sys.node_13.q => 0.0]
prob = remake(prob, p=params);
sol = solve(prob, RK4(), dt=dt, saveat=tstops)
```

## Using the solution
You can query the solution using symbolic indexing, e.g., to check if a pressure constraint is crossed
```
pressure_crossed =
    ( get_metadata(net.vertices[1], :min_pressure) >
      get_metadata(net, :pressure_from_density)(sol[sys.node_1.ρ][end]) )
if (pressure_crossed)
    # something bad
end
```
If `sol` is long, e.g. many timepoints, this symbolic indexing `sol[sys.node_1.ρ]` is very slow. Instead, use utilities from the `SymbolicIndexingInterface` package:
```
node1_ρ_inds = variable_index(sys, sys.node_1.ρ)
sol[node1_ρ_inds]
edge1_ρ_inds = [variable_index(sys, sys.edge_1.ρ[i]) for i in 1:length(sys.edge_1.ρ)]
```

## Re-running the simulation
If you want to iteratively run the simulation, just remake the problem
```
params = [sys.node_1.q => 0.0,
      sys.node_2.q => 0.0,
      sys.node_3.q => 0.0,
      sys.node_4.q => 0.0,
      sys.node_5.q => 0.0,
      sys.node_6.q => 0.0,
      sys.node_7.q => 0.0,
      sys.node_8.q => 0.0,
      sys.node_9.q => 0.0,
      sys.node_10.q => 0.0,
      sys.node_11.q => 0.0,
      sys.node_12.q => 0.0,
      sys.node_13.q => 0.0]

new_prob = remake(prob, u0 = sol[end], p=params)
```
