using OrdinaryDiffEq
using Trixi
using Plots



t = 0.4
coordinates_min = (-1.0,)
coordinates_max = ( 1.0,)

##############################################################################
# semidiscretization of the compressible Euler equations

equations = EulerEquations1D(5/3)

initial_condition = initial_condition_constant

boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = (x_neg=boundary_condition, x_pos=boundary_condition)


surface_flux = flux_lax_friedrichs
volume_flux  = flux_kennedy_gruber

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

solver = DGSEM(basis, surface_flux, volume_integral)



mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=6, n_cells_max=10_000, periodicity=false)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions)

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
sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=1.0, save_everystep=false, callback=callbacks);

summary_callback() # print the timer summary

pds = PlotData1D(sol; solution_variables=cons2prim)

plot(pds, linecolor ="red", guidefontsize = 10, legendfontsize= 10, titlefontsize = 25, tickfontsize=11, linewidth = 3, guidefont=font(19),size=(700,500))