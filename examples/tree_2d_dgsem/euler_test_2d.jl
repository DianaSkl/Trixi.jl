using OrdinaryDiffEq
using Trixi
using Plots

###############################################################################
# semidiscretization of the compressible Euler equations

equations = EulerEquations2D(5/3)

initial_condition = initial_condition_convergence_test



solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)


coordinates_min = (0.0, 0.0)
coordinates_max = (6.0, 6.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=4,
                n_cells_max=10_000)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms=source_terms_convergence_test)


###############################################################################
# ODE solvers, callbacks etc.
t = 0.1
tspan = (0.0, t)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_restart = SaveRestartCallback(interval=100,
                                   save_final_restart=true)

save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=0.1)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_restart, save_solution,
                        stepsize_callback)
###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary

pd1 = PlotData1D(sol; solution_variables=cons2prim)
plot!(pd1.x, pd1.data[:,4], xlims = (0.0,6.0), label = "2D Euler t = " *string(t), title ="rho, t = 0.1")