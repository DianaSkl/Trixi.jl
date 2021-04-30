
struct IsothermalEulerEquations1D{RealT<:Real} <: AbstractIsothermalEulerEquations{1, 2}
  gamma1::RealT
end


varnames(::typeof(cons2cons), ::IsothermalEulerEquations1D) = ("rho", "rho_v1")
varnames(::typeof(cons2prim), ::IsothermalEulerEquations1D) = ("rho", "v1")



function initial_condition_constant(x, t, equations::IsothermalEulerEquations1D)
  rho = 1.0
  rho_v1 = 0.1
  #rho_e = 10.0
  return SVector(rho, rho_v1)
end

function initial_condition_convergence_test(x, t, equations::IsothermalEulerEquations1D)
  c = 2
  A = 0.1
  L = 2
  f = 1/L
  ω = 2 * pi * f
  ini = c + A * sin(ω * (x[1] - t))

  rho = ini
  rho_v1 = ini
  #rho_e = ini^2

  return SVector(rho, rho_v1)
end


@inline function flux(u, orientation, equations::IsothermalEulerEquations1D)
  rho, rho_v1 = u
  v1 = rho_v1/rho
  #p = (equations.gamma1 - 1) * (rho_e - 1/2 * rho * v1^2)
  # Ignore orientation since it is always "1" in 1D
  f1 = rho_v1
  f2 = rho_v1 * v1 + rho
#  f3 = (rho_e + p) * v1
  return SVector(f1, f2)
end

@inline function max_abs_speeds(u, equations::IsothermalEulerEquations1D)
    rho, rho_v1 = u
    v1 = rho_v1 / rho
    #p = (equations.gamma - 1) * (rho_e - 1/2 * rho * v1^2)
    #c = sqrt(equations.gamma * p / rho)
  
    return (abs(v1),)
  end

function flux_lax_friedrichs(u_ll, u_rr, orientation, equations::IsothermalEulerEquations1D)
  # Calculate primitive variables and speed of sound
  rho_ll, rho_v1_ll = u_ll
  rho_rr, rho_v1_rr = u_rr

  v1_ll = rho_v1_ll / rho_ll
  v_mag_ll = abs(v1_ll)
  #p_ll = (equations.gamma - 1) * (rho_e_ll - 1/2 * rho_ll * v_mag_ll^2)
  #c_ll = sqrt(equations.gamma * p_ll / rho_ll)
  v1_rr = rho_v1_rr / rho_rr
  v_mag_rr = abs(v1_rr)
  #p_rr = (equations.gamma - 1) * (rho_e_rr - 1/2 * rho_rr * v_mag_rr^2)
  #c_rr = sqrt(equations.gamma * p_rr / rho_rr)

  # Obtain left and right fluxes
  f_ll = flux(u_ll, orientation, equations)
  f_rr = flux(u_rr, orientation, equations)

  λ_max = max(v_mag_ll, v_mag_rr)
  f1 = 1/2 * (f_ll[1] + f_rr[1]) - 1/2 * λ_max * (rho_rr    - rho_ll)
  f2 = 1/2 * (f_ll[2] + f_rr[2]) - 1/2 * λ_max * (rho_v1_rr - rho_v1_ll)
  #f3 = 1/2 * (f_ll[3] + f_rr[3]) - 1/2 * λ_max * (rho_e_rr  - rho_e_ll)

  return SVector(f1, f2)
end