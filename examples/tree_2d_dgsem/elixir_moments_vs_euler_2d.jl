using OrdinaryDiffEq
using Trixi
using Plots

t = 0.4
min = -1.0
max =  1.0
coordinates_min = (min, min)
coordinates_max = (max, max)
###############################################################################
# semidiscretization of the compressible Euler equations

equations = EulerEquations2D(5/3)

initial_condition = initial_condition_constant

surface_flux = flux_lax_friedrichs

boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = (x_neg=boundary_condition,
                       x_pos=boundary_condition,
                       y_neg=boundary_condition,
                       y_pos=boundary_condition,)

volume_flux  = flux_ranocha
basis = LobattoLegendreBasis(3)                                  
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)


mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=4, n_cells_max=10_000, periodicity=false)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions, source_terms=source_terms_convergence_test)



tspan = (0.0, t)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

analysis_interval = 100

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval=100)
alive_callback = AliveCallback(analysis_interval=analysis_interval)


# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval=100,solution_variables=cons2prim)

# The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl=0.3)

save_restart = SaveRestartCallback(interval=100,save_final_restart=true)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution, stepsize_callback, save_restart)



# ###############################################################################
# # run the simulation
sol = solve(ode, SSPRK43(), save_everystep=false, callback=callbacks);


summary_callback() # print the timer summary


pde = PlotData1D(sol; solution_variables=cons2prim)


########################### Perturbation System ###################################


# vxr = 0.0
# vyr = 0.0
# theta_r = 1.0
# rho_r = 0.1 
# tau = 0.001
# equations = MomentSystem2D(vxr, vyr,  theta_r, rho_r, tau)

# initial_condition = initial_condition_constant


# boundary_condition = BoundaryConditionDirichlet(initial_condition)
# boundary_conditions = (x_neg=boundary_condition,
#                        x_pos=boundary_condition,
#                        y_neg=boundary_condition,
#                        y_pos=boundary_condition,)

# surface_flux = flux_lax_friedrichs
# volume_flux  = flux_ranocha
# basis = LobattoLegendreBasis(3)


# volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

# solver = DGSEM(basis, surface_flux, volume_integral)




# mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=5, n_cells_max=10_000, periodicity=false)
# semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions, source_terms=source_terms_convergence_test)


# tspan = (0.0, t)
# ode = semidiscretize(semi, tspan)

# # At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# # and resets the timers
# summary_callback = SummaryCallback()

# analysis_interval = 100

# # The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
# analysis_callback = AnalysisCallback(semi, interval=100)
# alive_callback = AliveCallback(analysis_interval=analysis_interval)


# # The SaveSolutionCallback allows to save the solution to a file in regular intervals
# save_solution = SaveSolutionCallback(interval=100,solution_variables=cons2prim)

# # The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
# #stepsize_callback = StepsizeCallback(cfl=0.3)

# save_restart = SaveRestartCallback(interval=100,save_final_restart=true)

# # Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
# callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution, save_restart)


# #plot(pd.x, pd.data[:,1], xlims = (-1.0, 1.0), label = "Euler", title ="rho, t = "*string(t)*", tau = "*string(tau), seriestype = :scatter, markersize=2)
# ###############################################################################
# # run the simulation
# sol = solve(ode, SSPRK43(), save_everystep=false, callback=callbacks);



# # # Print the timer summary
# # summary_callback()

# pd2 = PlotData1D(sol; solution_variables=cons2prim)

# plot(pde.x, pde.data[:,1], xlims = (-1, 1), label = "Euler", seriestype = :scatter, markersize=3)
# # #plot(pd.x, pd.data[:,1], label = "Momentensystem", title ="ρ, t = "*string(t)*", τ = "*string(tau), linewidth = 3, guidefont=font(28))
# # #plot(pd.x, pd.data[:,1], label = "Momentensystem", title ="ρ, t = "*string(t)*", τ \u2192 ∞", linewidth = 3)

plot!(pde)
#plot!(pd2.x, pd2.data[:,1], label = "Momentensystem", title ="ρ, t = "*string(t)*", τ = "*string(tau), linewidth = 3, guidefont=font(28))

