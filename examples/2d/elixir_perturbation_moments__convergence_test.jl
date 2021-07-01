using OrdinaryDiffEq
using Trixi

equations = PerturbationMomentSystem2D(2.0, 2.0, 4.0)

initial_condition = initial_condition_convergence_test

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (-2, -2)
coordinates_max = ( 2,  2)
mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=3, n_cells_max=10_000)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms = source_terms_convergence_test)

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)


# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

analysis_interval = 100

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval=100)

# das hier sorgt für die Ausgabe!
alive_callback = AliveCallback(analysis_interval=analysis_interval)


# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval=100,solution_variables=cons2prim)

# The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl=0.1)

save_restart = SaveRestartCallback(interval=100,
                                   save_final_restart=true)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution, stepsize_callback, save_restart)



###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=1.0, save_everystep=false, callback=callbacks);

# Print the timer summary
summary_callback()