using OrdinaryDiffEq
using Trixi

equations = PerturbationMomentSystem2D(2.0, 2.0, 2.0)


initial_condition = initial_condition_constant 

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (0, 0)
coordinates_max = (1, 1)
mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=4, n_cells_max=10_000)
#semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms = source_terms_convergence_test)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

tspan = (0.0, 400.0)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval=100)

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval=100,solution_variables=cons2prim)

# The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl=0.8)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, save_solution, analysis_callback, stepsize_callback)



###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=1.0, save_everystep=false, callback=callbacks);
#sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=1.0, save_everystep=false);
# Print the timer summary
summary_callback()