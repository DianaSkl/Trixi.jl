using OrdinaryDiffEq
using Trixi
using Plots

tau = 0.001
equations = PerturbationMomentSystem1D(0.0, 1.0, tau)

initial_condition = initial_condition_constant
surface_flux = flux_lax_friedrichs

boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = (x_neg=boundary_condition, x_pos=boundary_condition)
volume_flux  = flux_lax_friedrichs
basis = LobattoLegendreBasis(3)
shock_indicator_variable = density_pressure
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.5,
                                         alpha_min=0.001,
                                         alpha_smooth=true,
                                         variable=shock_indicator_variable)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)
volume_integral = VolumeIntegralFluxDifferencing(flux_central)
solver = DGSEM(basis, surface_flux, volume_integral)
#solver = DGSEM(basis, surface_flux)
coordinates_min = (-1.0,)
coordinates_max = ( 1.0,)

mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=6, n_cells_max=10_000, periodicity=false)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions, source_terms=source_terms_convergence_test)

t = 0.1
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
#stepsize_callback = StepsizeCallback(cfl=0.9)

save_restart = SaveRestartCallback(interval=100,save_final_restart=true)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution, save_restart)



###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(), save_everystep=false, callback=callbacks);



# Print the timer summary

summary_callback()

pd = PlotData1D(sol; solution_variables=cons2prim)
pd2 = PlotData1D(sol; solution_variables=cons2cons)
plot(pd2)
#savefig("schwacheform.png")
#plot(pd.x, pd.data[:,1], xlims = (-1.0, 1.0), label = "DG schwach", title ="rho, t = "*string(t)*", tau = "*string(tau), seriestype = :scatter, markersize=3)


#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)
#solver = DGSEM(basis, surface_flux)
coordinates_min = (-1.0,)
coordinates_max = ( 1.0,)

mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=6, n_cells_max=10_000, periodicity=false)
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
#stepsize_callback = StepsizeCallback(cfl=0.9)

save_restart = SaveRestartCallback(interval=100,save_final_restart=true)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution, save_restart)



###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(), save_everystep=false, callback=callbacks);



# Print the timer summary

summary_callback()

plot(pd2)

#savefig("schwacheform.png")

pd = PlotData1D(sol; solution_variables=cons2prim)
pd2 = PlotData1D(sol; solution_variables=cons2cons)
plot!(pd.x, pd.data[:,1], xlims = (-1.0, 1.0), label = "DG stark", title ="rho, t = "*string(t)*", tau = "*string(tau), seriestype = :scatter, markersize=3)
savefig("stark-schwach.png")

plot(pd2)
savefig("starkeform.png")