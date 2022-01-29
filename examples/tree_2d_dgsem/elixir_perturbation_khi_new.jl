using OrdinaryDiffEq
using Trixi
using Plots 


###############################################################################
# semidiscretization of the compressible Euler equations
vxr = -0.002
vyr = 0.0
theta_r = 1.18
rho_r = 1.24
tau = 0.0001


equations = PerturbationMomentSystem2D(vxr, vyr, theta_r, rho_r, tau)

"""
    initial_condition_kelvin_helmholtz_instability(x, t, equations::CompressibleEulerEquations2D)
A version of the classical Kelvin-Helmholtz instability based on
- Andrés M. Rueda-Ramírez, Gregor J. Gassner (2021)
  A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations
  of the Euler Equations
  [arXiv: 2102.06017](https://arxiv.org/abs/2102.06017)
"""
function initial_condition_kelvin_helmholtz_instability(x, t, equations::PerturbationMomentSystem2D)
    @unpack vxr, vyr, theta_r, rho_r, tau = equations

    slope = 15
    amplitude = 0.02
    B = tanh(slope * x[2] + 7.5) - tanh(slope * x[2] - 7.5)
    rho = 0.5 + 0.75 * B
    vx = 0.5 * (B - 1)
    vy = 0.1 * sin(2 * pi * x[1])
    p = 1.0
    theta = p/rho



    drho = rho - rho_r
    dv_x = vx - vxr 
    dv_y = vy - vyr 
    dtheta = theta - theta_r
    
    diff_theta_y = -(0.75*(15*sech(15*x[2]+7.5)^2-15*sech(15*x[2]-7.5)^2))/(0.75*(tanh(15*x[2]+7.5)-tanh(15*x[2]-7.5))+0.5)^2
    diff_vx_x = 0
    diff_vx_y = 0.5*(15*sech(15*x[2]+7.5)^2-15*sech(15*x[2]-7.5)^2)
    diff_vy_x = 0.2*pi*cos(2*pi*x[1])
    diff_vy_y = 0
  
    # sigma_xx = -2*(5*diff_vx_x/3 - diff_vy_y/3)*tau
    # sigma_yy = -2*(5*diff_vy_y/3 - diff_vx_x/3)*tau
    # sigma_xy = -2*(diff_vy_x + diff_vx_y - (diff_vx_x + diff_vy_y)/3)*tau
    sigma_xx = (5*diff_vx_x/3 - diff_vy_y/3)*tau
    sigma_yy = (5*diff_vy_y/3 - diff_vx_x/3)*tau
    sigma_xy = (diff_vy_x + diff_vx_y - (diff_vx_x + diff_vy_y)/3)*tau

    q_x = 0
    #q_y = tau*(-15*diff_theta_y)/(4)
    q_y = 15*tau*diff_theta_y/4
   

    w0 = 1 + drho / rho_r
    w0x = dv_x / sqrt(theta_r) + (drho * dv_x)/(rho_r * sqrt(theta_r))
    w0y = dv_y / sqrt(theta_r) + (drho * dv_y)/(rho_r * sqrt(theta_r))
    w1 = - (rho * (dv_x *dv_x + dv_y * dv_y) )/(3.0 * (rho_r * theta_r)) - (drho * dtheta)/(rho_r * theta_r) - dtheta / theta_r
    w0xx = 0.5 * sigma_xx/(rho_r * theta_r) + (2.0 * dv_x * dv_x - dv_y * dv_y)/(6.0 * theta_r) + (2.0 * drho * dv_x * dv_x - drho * dv_y * dv_y)/(6.0 * rho_r * theta_r)  
    w0yy = 0.5 * sigma_yy/(rho_r * theta_r) + (2.0 * dv_y * dv_y - dv_x * dv_x)/(6.0 * theta_r) + (2.0 * drho * dv_y * dv_y - drho * dv_x * dv_x)/(6.0 * rho_r * theta_r)
    w0xy = 0.5 * sigma_xy/(rho_r * theta_r) + 0.5 * ((rho_r + drho) * (dv_x * dv_y))/(rho_r * theta_r)
    w1x = - 2.0 * q_x / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xx * dv_x + sigma_xy * dv_y))/(5.0*rho_r* sqrt(theta_r).^3.0) - (dtheta * dv_x * rho)/ (rho_r * sqrt(theta_r).^3.0) - rho * (dv_x * dv_y^2 + dv_x^3)/(5.0 * rho_r * sqrt(theta_r).^3.0)
    w1y = - 2.0 * q_y / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xy * dv_x + sigma_yy * dv_y))/(5.0 * rho_r * sqrt(theta_r).^3.0) - (dtheta * dv_y * rho)/ (rho_r * sqrt(theta_r).^3.0) - rho * (dv_y * dv_x^2 + dv_y^3)/(5.0 * rho_r * sqrt(theta_r).^3.0) 
    
    return SVector(w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y)
end
  

initial_condition = initial_condition_kelvin_helmholtz_instability

surface_flux = flux_lax_friedrichs
volume_flux  = flux_lax_friedrichs
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

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (-1.0, -1.0)
coordinates_max = ( 1.0,  1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=4,
                n_cells_max=100_000)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms=source_terms_convergence_test)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=20,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=0.4)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)


###############################################################################
# run the simulation

# sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
#             dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep=false, callback=callbacks);
sol = solve(ode, SSPRK43(), save_everystep=false, callback=callbacks);


summary_callback() # print the timer summary


pdt1 = PlotData1D(sol; solution_variables=cons2prim)
pdt2 = PlotData2D(sol; solution_variables=cons2prim)
plot(pdt2)
