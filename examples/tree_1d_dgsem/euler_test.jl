using OrdinaryDiffEq
using Trixi
using Plots

###############################################################################
# semidiscretization of the compressible Euler equations

# tau = 0.01
# equations = PerturbationMomentSystem1D(0.1, 2/3, tau)
# initial_condition = initial_condition_convergence_test
# surface_flux = flux_lax_friedrichs

# boundary_condition = BoundaryConditionDirichlet(initial_condition)
# boundary_conditions = (x_neg=boundary_condition, x_pos=boundary_condition)
# volume_flux  = flux_lax_friedrichs
# basis = LobattoLegendreBasis(3)

# volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
# solver = DGSEM(basis, surface_flux, volume_integral)
# coordinates_min = 0
# coordinates_max = 2

# mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=6, n_cells_max=10_000)
# semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, 
#                                     boundary_conditions=boundary_conditions, source_terms=source_terms_convergence_test)


# ###############################################################################
# # ODE solvers, callbacks etc.

# t = 0.2
# tspan = (0.0, t)
# ode = semidiscretize(semi, tspan)
# # At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# # and resets the timers

# summary_callback = SummaryCallback()

# # The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
# analysis_interval = 100
# analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

# # The AliveCallback prints short status information in regular intervals
# alive_callback = AliveCallback(analysis_interval=analysis_interval)

# # The SaveRestartCallback allows to save a file from which a Trixi simulation can be restarted
# save_restart = SaveRestartCallback(interval=100,
#                                    save_final_restart=true)

# # The SaveSolutionCallback allows to save the solution to a file in regular intervals
# save_solution = SaveSolutionCallback(interval=100, solution_variables=cons2prim)

# # The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
# stepsize_callback = StepsizeCallback(cfl=0.4)

# # Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
# callbacks = CallbackSet(summary_callback,
#                         analysis_callback, alive_callback,
#                         save_restart, save_solution, stepsize_callback )


# ###############################################################################
# # run the simulation

# # OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
# sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
#             dt=0.001, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep=false, callback=callbacks);

# # Print the timer summary
# summary_callback()

# pd = PlotData1D(sol; solution_variables=cons2prim)
# plot(pd.x, pd.data[:,1], xlims = (coordinates_min, coordinates_max), label = "Momentensystem", title ="rho, t = "*string(t)*", tau = "*string(tau), seriestype = :scatter, markersize=2)


#######################################################################################################################################################
####################################### Euler #########################################################################################################


# equations = EulerEquations1D(5/3)

# initial_condition = initial_condition_convergence_test
# surface_flux = flux_lax_friedrichs

# boundary_condition = BoundaryConditionDirichlet(initial_condition)
# boundary_conditions = (x_neg=boundary_condition, x_pos=boundary_condition)
# #volume_flux  = flux_lax_friedrichs
# #basis = LobattoLegendreBasis(3)

# #volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
# #solver = DGSEM(basis, surface_flux, volume_integral)

# solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)
# coordinates_min = 0
# coordinates_max = 2

# mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=5, n_cells_max=10_000, periodicity=true)
# semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms=source_terms_convergence_test, boundary_conditions=boundary_conditions)

# ###############################################################################
# # ODE solvers, callbacks etc.

# t = 0.5
# tspan = (0.0, t)
# ode = semidiscretize(semi, tspan)
# # At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# # and resets the timers

# summary_callback = SummaryCallback()

# # The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
# analysis_interval = 100
# analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

# # The AliveCallback prints short status information in regular intervals
# alive_callback = AliveCallback(analysis_interval=analysis_interval)

# # The SaveRestartCallback allows to save a file from which a Trixi simulation can be restarted
# save_restart = SaveRestartCallback(interval=100,
#                                    save_final_restart=true)

# # The SaveSolutionCallback allows to save the solution to a file in regular intervals
# save_solution = SaveSolutionCallback(interval=100, solution_variables=cons2prim)

# # The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
# stepsize_callback = StepsizeCallback(cfl=0.01)

# # Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
# callbacks = CallbackSet(summary_callback,
#                         analysis_callback, alive_callback,
#                         save_restart, save_solution, stepsize_callback )


# ###############################################################################
# # run the simulation

# # OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
# sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
#             dt=0.001, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep=false, callback=callbacks);

# # Print the timer summary
# summary_callback()


# pd = PlotData1D(sol; solution_variables=cons2prim)


# #plot!(pd.x, pd.data[:,2], xlims = (coordinates_min, coordinates_max), label = "Euler t= "*string(t), markersize=3)
# plot(pd)
using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations1D(1.4)

initial_condition = initial_condition_convergence_test

# Note that the expected EOC of 5 is not reached with this flux.
# Using flux_hll instead yields the expected EOC.
solver = DGSEM(polydeg=3, surface_flux=flux_hll)

coordinates_min = 0.0
coordinates_max = 2.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=1,
                n_cells_max=10_000)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms=source_terms_convergence_test)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.5)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
                                     extra_analysis_errors=(:l2_error_primitive,
                                                            :linf_error_primitive))

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=0.8)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)


###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary

pd = PlotData1D(sol; solution_variables=cons2prim)
plot(pd)