  struct EulerEquations1D{RealT<:Real} <: AbstractCompressibleEulerEquations{1, 4}
    theta::RealT 
  end
  

  varnames(::typeof(cons2prim), ::EulerEquations1D) = ("rho", "v1", "v2", "theta")
  varnames(::typeof(cons2cons), ::EulerEquations1D) = ("rho", "rho_v1", "rho_v2", "rho_e")
  
  

  function initial_condition_constant(x, t, equations::EulerEquations1D)
   
    if (x[1] < 0)
      rho = 3.0 
    else
      rho = 1.0
    end
    
    v1 = 0.0
    v2 = 0.0 
    theta = 1.0

    rho_e = 0.5 * rho * (v1^2 + v2^2) + 3 * rho * theta /2 

    return SVector(rho, rho*v1, rho*v2, rho_e)
  end
  
  @inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::EulerEquations1D)
    v1 = 0.0
    theta = 1.0
    
    c = sqrt(5.0*theta/3.0)

    return abs(v1)+c
  end
  
  @inline function max_abs_speeds(u, equations::EulerEquations1D)
    v1 = 0.0
    theta = 1.0
    
    c = sqrt(5.0*theta/3.0)
  
    return abs(v1) + c
  end
  
  
  # Convert conservative variables to primitive
  # @inline function cons2prim(u, equations::EulerEquations1D)
  #   rho, v1, v2,  theta = u
    
  #   return SVector(rho, v1, v2, theta)
  # end
  
  @inline function cons2prim(u, equations::EulerEquations1D)
    rho, v1, v2, rho_e = u
  
    e = rho_e/rho
    theta = 2*(e - (v1^2 +v2^2)/2)/3

    return SVector(rho, v1, v2, theta)
  end

  
  # Convert conservative variables to entropy
  @inline function cons2entropy(u, equations::EulerEquations1D)
    rho, v1, v2,  theta = u
  
    return SVector(rho, v1, v2, theta)
  end
  


  @inline function density(u, equations::EulerEquations1D)
   rho = u[1]
   return rho
  end
  
  
  @inline function density_pressure(u, equations::EulerEquations1D)
    rho = u[1]
    return rho
   end

  
  
 #Calculate 1D flux for a single point
@inline function flux(u, orientation::Integer, equations::EulerEquations1D)
  rho, v1, v2, rho_e = u

  e = rho_e/rho
  theta = 2*(e - (v1^2 +v2^2)/2)/3

  f1 = rho * v1
  f2 = rho * v1 * v1 + rho * theta
  f3 = rho * v1 * v2
  f4 = v1* ( rho*(v1*v1 + v2*v2)/2 + 5*rho*theta/2)
  return SVector(f1, f2, f3, f4)
end



@inline function source_terms_convergence_test(u, x, t, equations::EulerEquations1D)

    
  return SVector(0.0, 0.0, 0.0, 0.0)

end