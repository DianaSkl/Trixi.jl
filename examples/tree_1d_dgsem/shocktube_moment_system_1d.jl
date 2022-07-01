using OrdinaryDiffEq
using Trixi
using Plots

t = 0.6
coordinates_min = (-1.0,)
coordinates_max = ( 1.0,)


##############################################################################
# semidiscretization of the Moment System


vxr = 0.0
theta_r = 1.0
rho_r = 1.0 
tau = 0.001
equations = MomentSystem1D(vxr, theta_r, rho_r, tau)

initial_condition = initial_condition_constant


boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = (x_neg=boundary_condition, x_pos=boundary_condition)

surface_flux = flux_lax_friedrichs
volume_flux  = flux_ds

basis = LobattoLegendreBasis(3)

shock_indicator_variable = density_pressure

indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.5,
                                         alpha_min=0.0001,
                                         alpha_smooth=true,
                                         variable=shock_indicator_variable)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc; volume_flux_dg=volume_flux, volume_flux_fv=surface_flux)


solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=6, n_cells_max=10_000, periodicity=false)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions,source_terms=source_terms_convergence_test)



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


###############################################################################
# run the simulation
#sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=1.0, save_everystep=false, callback=callbacks);
sol = solve(ode, SSPRK43(), save_everystep=false, callback=callbacks)

# Print the timer summary
summary_callback()

###############################################################################
# plot the simulation

pd = PlotData1D(sol; solution_variables=cons2prim)
#  f1 = plot(pd.x, pd.data[:,1], title = "ρ")
# f2 = plot(pd.x, pd.data[:,2], title = "v\u2093")
#  f3 = plot(pd.x, pd.data[:,3], title = "p")
#plot!(f1, f2, f3, legend = false, layout=(1,3),titlefontsize = 21, tickfontsize=12, linewidth = 2, size=(900,500))
#plot!(pd,  size=(1000,1000),  linewidth = 3)

plot!(pd.x, pd.data[:,4] , title= L"\sigma_{xx}", label="t = 0.6", legend=:topleft, guidefontsize = 20, legendfontsize= 20, titlefontsize = 30, tickfontsize=18, linewidth = 3, size=(900,500))
#plot!(pd.x, pd.data[:,5] , title= L"q_x", label="t = 0.6", legend=:bottomleft, guidefontsize = 20, legendfontsize= 20, titlefontsize = 30, tickfontsize=18, linewidth = 3, size=(900,500))
#plot(pd.x, pd.data[:,1] , title= "ρ", legendfontsize= 18, titlefontsize = 25, tickfontsize=12, linewidth = 3, size=(900,500))
