using OrdinaryDiffEq
using Trixi
using Plots

tau = 200.0
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
solver = DGSEM(basis, surface_flux, volume_integral)
coordinates_min = (-1.0,)
coordinates_max = ( 1.0,)

mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=8, n_cells_max=10_000, periodicity=false)
#semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions=boundary_conditions, source_terms=source_terms_convergence_test)

t = 0.4
tspan = (0.0, t)
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

analysis_interval = 100

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval=100)

# das hier sorgt für die Ausgabe!
alive_callback = AliveCallback(analysis_interval=analysis_interval)

# amr
amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable=first),
                                      base_level=4,
                                      med_level=5, med_threshold=0.1,
                                      max_level=6, max_threshold=0.3)
amr_callback = AMRCallback(semi, amr_controller,
                           interval=3,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)


# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval=100,solution_variables=cons2prim)

# The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl=0.9)

save_restart = SaveRestartCallback(interval=100,save_final_restart=true)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution, stepsize_callback,amr_callback, save_restart)



###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=1.0, save_everystep=false, callback=callbacks);

# Print the timer summary
summary_callback()

pd = PlotData1D(sol; solution_variables=cons2prim)

plot(pd.x, pd.data[:,1], xlims = (-1.8,), ylims = (1.0, 5.1), label = "DGSEM", title ="rho, t = "*string(t)*", tau = "*string(tau)*", B ="*string(1/tau), seriestype = :scatter, markersize=1)
# bildx =[-1.70073327222731,	-1.63473877176902,	-1.55499541704858,	-1.480751604033,	-1.40650779101742,	-1.33776351970669,	-1.26076993583868,	-1.17277726856095,	-1.07103574702108,	-0.972043996333639,	-0.867552703941338,	-0.774060494958754,	-0.67781851512374,	-0.578826764436297,	-0.507332722273144,	-0.441338221814849,	-0.378093492208983,	-0.325847846012832,	-0.270852428964253,	-0.210357470210816,	-0.152612282309808,	-0.0976168652612281,	-0.037121906507791,	0.00412465627864345,	0.0371219065077908,	0.10724106324473,	0.185609532538955,	0.262603116406966,	0.336846929422548,	0.400091659028414,	0.457836846929423,	0.488084326306141,	0.504582951420715,	0.559578368469294,	0.644821264894593,	0.768560953253895,	0.897800183318057,	1.04353803849679,	1.06278643446379,	1.07378551787351,	1.0820348304308,	1.16452795600367,	1.24977085242896,	1.34601283226398,	1.41200733272227,	1.50274977085243,	1.57149404216315,	1.64573785517874,	1.74472960586618,	1.76397800183318]
# bildy =[5.00679873273824,	5.01254049650017,	4.99420318577056,	4.98189390123396,	4.94556059267333,	4.8912312862817,	4.84287495561556,	4.7824625725634,	4.69196510901369,	4.60748466156349,	4.52898819993229,	4.44452977266909,	4.32402429927636,	4.22152583380814,	4.13715548729298,	4.07082517898283,	3.96246383872415,	3.7760684607614,	3.61368609672918,	3.42725768848592,	3.21081026030614,	3.00037984822586,	2.78992741595858,	2.60357607836985,	2.48332383712769,	2.42298301968329,	2.34759415694338,	2.28422281126222,	2.22686848168058,	2.19057921349397,	2.16932698050938,	2.03707373734871,	1.83280347258349,	1.82057125870141,	1.78419390976678,	1.71763238949307,	1.65705485503836,	1.57839324200461,	1.38612397915789,	1.19989375259769,	1.1157766382331,	1.10343432341599,	1.10009000751439,	1.09670165123877,	1.09043140298869,	1.09006806990308,	1.0867898145625,	1.08348953903491,	1.07708716966279,	1.02896205096022]
# plot!(bildx, bildy)