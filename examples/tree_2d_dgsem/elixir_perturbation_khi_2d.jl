using OrdinaryDiffEq
using Trixi
using Plots 
using TickTock

###############################################################################
# semidiscretization of the compressible Euler equations
vxr = 0.2
vyr = 0.1
theta_r = 1.18
rho_r = 1.24
tau = 0.001
equations = MomentSystem2D(vxr, vyr, theta_r, rho_r, tau)


tick()

"""
    initial_condition_kelvin_helmholtz_instability(x, t, equations::CompressibleEulerEquations2D)
A version of the classical Kelvin-Helmholtz instability based on
- Andrés M. Rueda-Ramírez, Gregor J. Gassner (2021)
  A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations
  of the Euler Equations
  [arXiv: 2102.06017](https://arxiv.org/abs/2102.06017)
"""
function initial_condition_kelvin_helmholtz_instability(x, t, equations::MomentSystem2D)
    # change discontinuity to tanh
    # typical resolution 128^2, 256^2
    # domain size is [-1,+1]^2
    slope = 15
    amplitude = 0.02
    B = tanh(slope * x[2] + 7.5) - tanh(slope * x[2] - 7.5)
    rho = 0.5 + 0.75 * B
    vx = 0.5 * (B - 1)
    vy = 0.1 * sin(2 * pi * x[1])
    p = 1.0
    return prim2cons(SVector(rho, vx, vy, p), equations)
end
  

initial_condition = initial_condition_kelvin_helmholtz_instability

surface_flux = flux_lax_friedrichs
volume_flux  = flux_kennedy_gruber
polydeg = 3
basis = LobattoLegendreBasis(polydeg)


indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.002,
                                         alpha_min=0.0001,
                                         alpha_smooth=true,
                                         variable=density_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)

#volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (-1.0, -1.0)
coordinates_max = ( 1.0,  1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=6,
                n_cells_max=400_000)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms=source_terms_convergence_test)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.5)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=20,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)                                   

stepsize_callback = StepsizeCallback(cfl=0.001)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution, stepsize_callback)


###############################################################################
# run the simulation

# sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
#             dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep=false, callback=callbacks);
sol = solve(ode, SSPRK43(), save_everystep=false, callback=callbacks, maxiters = 1e25);


summary_callback() # print the timer summary

# pdt = PlotData1D(sol; solution_variables=cons2prim)
# plot(pdt)


tock()


# pdt1 = PlotData1D(sol; solution_variables=cons2prim)
 pdt2 = PlotData2D(sol; solution_variables=cons2prim)
 plot(pdt2, size = (1900,1200),  titlefontsize = 21, tickfontsize=12)

# #plot(pdt2["ρ"], title = "ρm", size = (1000,800))
# plot(pdt2["v1"], title = L"v_x")
# plot(pdt2["v2"], title = L"v_y")
# plot(pdt2["p"], title = L"p")
