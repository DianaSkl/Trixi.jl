using OrdinaryDiffEq
using Trixi

equation = Trixi.IsothermalEulerEquations1D(1.4)
initial_condition = initial_condition_constant
 surface_flux = flux_lax_friedrichs

solver = DGSEM(3, surface_flux)
coordinates_min = (-1,)
coordinates_max = ( 1,)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=4,
                n_cells_max=30_000)
semi = SemidiscretizationHyperbolic(mesh, equation, initial_condition, solver)
tspan = (0.0, 0.4)
ode = semidiscretize(semi, tspan)

# # ###############################################################################
# # # run the simulation

# sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
#             dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep=false, maxiters=1e5);
