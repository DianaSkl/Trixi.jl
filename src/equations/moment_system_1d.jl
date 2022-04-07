struct MomentSystem1D{RealT<:Real} <: AbstractMomentSystem{1, 5} 
    vxr::RealT
    theta_r::RealT
    rho_r::RealT
    tau::RealT
  end

  
varnames(::typeof(cons2prim), ::MomentSystem1D) = ("rho", "vx", "theta")
varnames(::typeof(cons2cons), ::MomentSystem1D) = ("w0", "w0x", "w1", "w0xx", "w1x" )

  
@inline function flux(u, orientation::Integer, equations::MomentSystem1D)
    w0, w0x, w1, w0xx, w1x = u
  
    @unpack vxr, theta_r = equations
    stheta = sqrt(theta_r)
    
    f1  = vxr * w0  + stheta * w0x
    f2  = vxr * w0x + stheta * (w0 - w1 + 2*w0xx)
    f3  = vxr * w1 + stheta * (5.0 * w1x - 2.0 * w0x)/3.0
    f4  = vxr * w0xx + 2.0 * stheta * (w0x - w1x)/3.0 
    f5  = vxr * w1x + stheta * (w1 - 4.0 * w0xx/5.0) 
    
    return SVector(f1, f2, f3, f4, f5)
  end

  
  @inline function cons2prim(prim, equations::MomentSystem1D)
    w0, w0x, w1 = prim
    @unpack vxr, theta_r, rho_r = equations
  
    vx = 1.0 

    dvx = vx - vxr
    rho = w0 * rho_r
    vx = vxr + w0x * sqrt(theta_r) / w0
    theta = theta_r*(1- w1/w0)- dvx^2/3
    
    return SVector(rho, vx, theta)
  end
  
    
  @inline function cons2entropy(u, equations::MomentSystem1D)
    w0, w0x, w1, w0xx, w1x= u
  
    return SVector(w0, w0x, w1, w0xx, w1x)
  end
  
  @inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::MomentSystem1D)
    @unpack vxr, theta_r = equations
    ab1 = abs(vxr) + 2.0 * sqrt(5 * theta_r / 3)
  
    return ab1
  end
  
  @inline function max_abs_speeds(u, equations::MomentSystem1D)
    @unpack vxr, theta_r = equations
    ab1 = abs(vxr) + 2.0 * sqrt(5 * theta_r / 3)
  
    return ab1
  end
  
   
  function initial_condition_convergence_test(x, t, equations::MomentSystem1D)
    
    @unpack vxr, theta_r, rho_r = equations 
    c = 2
    A = 0.1
    L = 2
    f = 1/L
    ω = 2 * pi * f
  
    ini = c + A * sin(ω * (x[1] - t))

    rho = ini
  
    vx = 1.0
  
    dv_x = vx - vxr
    theta = 2 * (ini - 0.5)/3
    dtheta = theta - theta_r

    sigma_xx = 0.1
    q_x = 0.1
  
    w0 = rho/rho_r
    w0x = (rho * dv_x)/(rho_r * sqrt(theta_r))
    w1 = - (dtheta * rho)/(rho_r * theta_r) - rho*dv_x^2/(3 * rho_r * theta_r) 
    w0xx = (sigma_xx + rho * dv_x^2 * 2 /3)/(2 * rho_r * theta_r)
    w1x =  - 2.0 * q_x / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xx * dv_x))/(5.0*rho_r* sqrt(theta_r).^3.0) - (dtheta * dv_x * rho)/ (rho_r * sqrt(theta_r).^3.0) - (rho*dv_x^3 /(5*rho_r*sqrt(theta_r).^3.0))
    
    return SVector(w0, w0x, w1, w0xx, w1x)
  end
  
  @inline function source_terms_convergence_test(u, x, t, equations::MomentSystem1D)
    
    @unpack vxr, theta_r, rho_r, tau = equations 
    c = 2
    A = 0.1
    L = 2
    f = 1/L
    ω = 2 * pi * f
    
    p1 = 2/3
    p2 = 1/3

    stheta = sqrt(theta_r)

    vx = 1.0
    dv_x = vx - vxr

    w0, w0x, w1, w0xx, w1x = u
  
    a1 = -w0 * w0xx + w0x * w0x/3.0
    a2 = -2.0 * w0 * w1x/3.0 + 2.0 * w1 * w0x/3.0 + 4.0 * w0x * w0xx/15.0 
  
    
   
    dw0_x = (A*ω*cos((x[1]-t)*ω))/c
    
    dw0x_x = (A*dv_x*ω*cos((x[1]-t)*ω))/(c*sqrt(theta_r))

    dw1_x = -(A*ω*cos((x[1]-t)*ω)*(p1*(A*sin((x[1]-t)*ω)+c)-theta_r-p2))/(c*theta_r)-(A*p1*ω*cos((x[1]-t)*ω)*(A*sin((x[1]-t)*ω)+c))/(c*theta_r)-(A*dv_x^2*ω*cos((x[1]-t)*ω))/(3*c*theta_r)
    
    dw0xx_x = (A*dv_x^2*p1*ω*cos((x[1]-t)*ω))/(c^2*theta_r)

    dw1x_x =-(A*dv_x*ω*cos((x[1]-t)*ω)*(p1*(A*sin((x[1]-t)*ω)+c)-theta_r-p2))/(c*theta_r^1.5)-(A*dv_x*p1*ω*cos((x[1]-t)*ω)*(A*sin((x[1]-t)*ω)+c))/(c*theta_r^1.5)-(A*dv_x^3*ω*cos((x[1]-t)*ω))/(5*c*theta_r^1.5)
    

    f1  = -dw0_x + vxr * dw0_x  + stheta * dw0x_x
    f2  = -dw0x_x + vxr * dw0x_x + stheta * (dw0_x - dw1_x + 2*dw0xx_x)
    f3  = -dw1_x + vxr * dw1_x + stheta * (5.0 * dw1x_x - 2.0 * dw0x_x)/3.0
    f4  = -dw0xx_x + vxr * dw0xx_x + 2.0 * stheta * (dw0x_x - dw1x_x)/3.0 
    f5  = -dw1x_x + vxr * dw1x_x + stheta * (dw1_x - 4.0 * dw0xx_x/5.0) 


    return(f1,f2,f3,f4,f5)
  end 