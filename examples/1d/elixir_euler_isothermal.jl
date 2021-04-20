module IsothermalEuler

using Trixi

struct IsothermalEulerEquations1D{RealT<:Real} <: Trixi.AbstractEquations{1 #= number of spatial dimensions =#,
                                                3 #= number of primary variables, i.e. scalar =#}
  gamma::RealT
end


varnames(::typeof(cons2cons), ::IsothermalEulerEquations1D) = ("rho", "rho_v1", "rho_e")
varnames(::typeof(cons2prim), ::IsothermalEulerEquations1D) = ("rho", "v1", "p")

function initial_condition_constant(x, t, equations::IsothermalEulerEquations1D)
  rho = 1.0
  rho_v1 = 0.1
  rho_e = 10.0
  return SVector(rho, rho_v1, rho_e)
end

@inline function flux(u, orientation, equations::IsothermalEulerEquations1D)
  rho, rho_v1, rho_e = u
  v1 = rho_v1/rho
  p = (equations.gamma - 1) * (rho_e - 1/2 * rho * v1^2)
  # Ignore orientation since it is always "1" in 1D
  f1 = rho_v1
  f2 = rho_v1 * v1 + p
  #f3 = (rho_e + p) * v1
  return SVector(f1, f2)
end

@inline function flux_chandrashekar(u_ll, u_rr, orientation, equations::IsothermalEulerEquations1D)
  # Unpack left and right state
  rho_ll, rho_v1_ll, rho_e_ll = u_ll
  rho_rr, rho_v1_rr, rho_e_rr = u_rr

  v1_ll = rho_v1_ll/rho_ll
  v1_rr = rho_v1_rr/rho_rr
  p_ll =  (equations.gamma - 1) * (rho_e_ll - 1/2 * rho_ll * (v1_ll^2))
  p_rr =  (equations.gamma - 1) * (rho_e_rr - 1/2 * rho_rr * (v1_rr^2))
  beta_ll = 0.5*rho_ll/p_ll
  beta_rr = 0.5*rho_rr/p_rr
  specific_kin_ll = 0.5*(v1_ll^2)
  specific_kin_rr = 0.5*(v1_rr^2)

  # Compute the necessary mean values
  rho_avg  = 0.5*(rho_ll+rho_rr)
  rho_mean = ln_mean(rho_ll,rho_rr)
  beta_mean = ln_mean(beta_ll,beta_rr)
  beta_avg = 0.5*(beta_ll+beta_rr)
  v1_avg = 0.5*(v1_ll+v1_rr)
  p_mean = 0.5*rho_avg/beta_avg
  velocity_square_avg = specific_kin_ll + specific_kin_rr

  # Calculate fluxes
  # Ignore orientation since it is always "1" in 1D
  f1 = rho_mean * v1_avg
  f2 = f1 * v1_avg + p_mean
  #f3 = f1 * 0.5*(1/(equations.gamma-1)/beta_mean - velocity_square_avg)+f2*v1_avg

  return SVector(f1, f2)
end



end




import .IsothermalEuler
using OrdinaryDiffEq

###############################################################################
# semidiscretization of the isothermal Euler equations

equation = IsothermalEuler.IsothermalEulerEquations1D(1.4)

initial_condition = IsothermalEuler.initial_condition_constant

 surface_flux = IsothermalEuler.flux_chandrashekar
 volume_flux  = IsothermalEuler.flux_chandrashekar
 solver = DGSEM(3, surface_flux, VolumeIntegralFluxDifferencing(volume_flux))

 semi = SemidiscretizationHyperbolic(mesh, equation, initial_condition, solver)
 tspan = (0.0, 0.4)
 ode = semidiscretize(semi, tspan)


###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, maxiters=1e5);
