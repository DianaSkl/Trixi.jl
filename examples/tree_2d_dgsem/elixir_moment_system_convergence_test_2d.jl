# This elixir was one of the setups used in the following master thesis:
# - Diana Sklema (2022)
#   "Untersuchung eines gestoerten Momentensystems in der kompressiblen Stroeomungsmechanik"
#   University of Cologne, advisors: Gregor Gassner, Michael Schlottke-Lakemper
#  Run "convergence_test(joinpath(examples_dir(), "tree_2d_dgsem", "elixir_moment_system_convergence_test_2d.jl"), 4)"
#  for the 2D convergence test starting with refinement 2

using OrdinaryDiffEq
using Trixi
using Plots

###############################################################################
# semidiscretization of the perturbated moment system

vxr = 0.5
vyr = 0.5
theta_r = 1/3
rho_r = 2.0
tau = 0.1
equations = MomentSystem2D(vxr, vyr, theta_r, rho_r, tau)

initial_condition = initial_condition_convergence_test


solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

min = 0
max = 2
coordinates_min = (min, min)
coordinates_max = (max, max)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=2,
                n_cells_max=100_000)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms=source_terms_convergence_test)


###############################################################################
# ODE solvers, callbacks etc.
t = 0.5
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

stepsize_callback = StepsizeCallback(cfl=0.8)

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


