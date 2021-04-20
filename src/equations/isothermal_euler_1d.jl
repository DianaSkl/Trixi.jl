
struct IsothermalEulerEquations1D{RealT<:Real} <: AbstractCompressibleEulerEquations{1, 3}
  gamma::RealT
end


varnames(::typeof(cons2cons), ::IsothermalEulerEquations1D) = ("rho", "rho_v1", "rho_e")
varnames(::typeof(cons2prim), ::IsothermalEulerEquations1D) = ("rho", "v1", "p")


#
# function initial_condition_constant(x, t, equations::IsothermalEulerEquations1D)
#   rho = 1.0
#   rho_v1 = 0.1
#   rho_e = 10.0
#   return SVector(rho, rho_v1, rho_e)
# end

# function initial_condition_convergence_test(x, t, equations::IsothermalEulerEquations1D)
#   c = 2
#   A = 0.1
#   L = 2
#   f = 1/L
#   ω = 2 * pi * f
#   ini = c + A * sin(ω * (x[1] - t))
#
#   rho = ini
#   rho_v1 = ini
#   rho_e = ini^2
#
#   return SVector(rho, rho_v1, rho_e)
# end

# Calculate 1D flux for a single point
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
