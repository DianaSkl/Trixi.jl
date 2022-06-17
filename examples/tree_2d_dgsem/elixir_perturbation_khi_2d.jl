using OrdinaryDiffEq
using Trixi
using Plots 
using TickTock

###############################################################################

#tau
tau = 1.0e-2
e = 2
#initial refinement
r = 5
#time
t = 1.5
#Shock Capturing Blending Factor
alpha = 0.008
###############################################################################

vxr = 0.2
vyr = 0.1
theta_r = 1.18
rho_r = 1.24

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
    slope = 15
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
                                         alpha_max=alpha,
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
                initial_refinement_level=r,
                n_cells_max=400_000)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms=source_terms_convergence_test)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, t)
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


#  fig1 = plot(pdt2["ρ"], title = "ρ", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig1, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_rho.png")
#  fig2 = plot(pdt2["vx"], title = L"v_x", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig2, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_vx.png")
#  fig3 = plot(pdt2["vy"], title = L"v_y", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig3, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_vy.png")
#  fig4 = plot(pdt2["p"], title = L"p", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig4, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_p.png")
#  fig5 = plot(pdt2["σxx"], title = L"\sigma_{xx}", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig5, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_sxx.png")
#  fig6 = plot(pdt2["σxy"], title = L"\sigma_{xy}", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig6, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_sxy.png")
#  fig7 = plot(pdt2["σyy"], title = L"\sigma_{yy}", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig7, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_syy.png")
#  fig8 = plot(pdt2["qx"], title = L"q_x", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig8, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_qx.png")
#  fig9 = plot(pdt2["qy"], title = L"q_y", size = (1000,800),  titlefontsize = 30, tickfontsize=25,guidefont=font(24))
#  savefig(fig9, "C:/Users/diana/OneDrive/Desktop/julianeu/alpha"*string(alpha)*"_t"*string(t*10)*"e"*string(e)*"r"*string(r)*"_qy.png")

