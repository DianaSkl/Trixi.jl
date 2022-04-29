using OrdinaryDiffEq
using Trixi
using Plots 


###############################################################################
# semidiscretization of the compressible Euler equations
vxr = 1.0
vyr = 0.004
theta_r = 7.4
rho_r = 2.0
tau = 0.02


equations = MomentSystem2D(vxr, vyr, theta_r, rho:r,tau)

function initial_condition_kelvin_helmholtz_instability(x, t, equations::MomentSystem2D)

    @unpack vxr, vyr, theta_r = equations

    a = 0.05
    z1 = 0.5
    z2 = 1.5
    s=0.2
    tmp = tanh((x[2]-z1)/a)- tanh((x[2]-z2)/a)
    rho = 1 + 0.5*tmp
    vx = tmp - 1
    vy = 0.01*sin(2 * pi * x[1])*(exp(-((x[2]-z1)^2)/s^2)+ exp(-((x[2]-z2)^2)/s^2))
    p = 10
    theta = p/rho

    rho_r = rho - 0.5

    drho = rho - rho_r
    dv_x = vx - vxr 
    dv_y = vy - vyr 
    dtheta = theta - theta_r
    
    mu = 2/10^5
    lambda = mu*rho

    dvy = 0.01*sin(2 * pi * x[1])*(exp(-((x[2]-z2)^2)/s^2)*(-2*(x[2]-z2)/s^2) + exp(-((x[2]-z1)^2)/s^2)*(-2*(x[2]-z1)/s^2))
    dvxy = (sech((x[2]-z1)/a))^2 - (sech((x[2]-z2)/a))^2 
    dvyx = 0.02*pi*cos(2*pi*x[1])*(exp(-(x[2]-z1)^2/s^2) + exp(-(x[2]-z2)^2/s^2))
    sigma_xx = mu*(-dvy/3)
    sigma_yy = mu*(dvyx + dvxy - dvy/3)
    sigma_xy = mu*(5*dvy/3)

    q_x = 0
    q_y = -lambda*(5.0*(sech((x[2]-z1)/a)^2/a-sech((x[2]-z2)/a)^2/a))/(0.5*(tanh((x[2]-z1)/a)-tanh((x[2]-z2)/a))+1)^2

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

#volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=5,
                n_cells_max=100_000)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.5)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=20,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=0.1)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)


###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
#sol = solve(ode, SSPRK43(), save_everystep=false, callback=callbacks);


summary_callback() # print the timer summary

# pdt = PlotData1D(sol; solution_variables=cons2prim)
# plot(pdt)

pdt1 = PlotData1D(sol; solution_variables=cons2prim)
pdt2 = PlotData2D(sol; solution_variables=cons2prim)
plot(pdt2)
#plot!(pdt1.x, pdt1.data[:,1], xlims = (-1.0, 1.0),  title ="rho")
#label = "MS rho_r= 2",*string(theta_r)