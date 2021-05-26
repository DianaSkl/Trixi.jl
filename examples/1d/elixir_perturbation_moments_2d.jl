using OrdinaryDiffEq
using Trixi

equations = PerturbationMomentSystem2D(2.0, 2.0, 2.0)

source_terms = source_terms_convergence_test

initial_condition = initial_condition_constant 

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (-5, -5)
coordinates_max = ( 5,  5)
mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=4, n_cells_max=30_000, periodicity=false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

tspan = (0.0, 5.0)
ode = semidiscretize(semi, tspan)