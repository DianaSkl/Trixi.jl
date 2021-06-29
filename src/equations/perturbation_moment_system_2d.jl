struct PerturbationMomentSystem2D{RealT<:Real} <: AbstractPerturbationMomentSystem2D{2, 9} 
    vxr::RealT
    vyr::RealT
    theta_r::RealT
  end
  
  varnames(::typeof(cons2prim), ::PerturbationMomentSystem2D) = ("w0", "w0x", "w0y", "w1", "w0xx", "w0yy", "w0xy", "w1x", "w1y")
  varnames(::typeof(cons2cons), ::PerturbationMomentSystem2D) = ("w0", "w0x", "w0y", "w1", "w0xx", "w0yy", "w0xy", "w1x", "w1y")
  
  @inline function flux(u, orientation::Integer, equations::PerturbationMomentSystem2D)
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
    @unpack vxr, vyr, theta_r = equations
  
    if orientation == 1
        
      f1  = vxr * w0  + sqrt(theta_r) * w0x
      f2  = vxr * w0x + sqrt(theta_r) * (w0 - w1) + 2.0 * sqrt(theta_r) * w0xx
      f3  = vxr * w0y + 2.0 * sqrt(theta_r) * w0xy 
      f4  = vxr * w1 + sqrt(theta_r) * (5.0 * w1x - 2.0 * w0x)/3.0
      f5  = vxr * w0xx + 2.0 * sqrt(theta_r)* (w0x - w1x)/3.0 
      f6  = vxr * w0xy + sqrt(theta_r) * (w0y - w1y)/2.0 
      f7  = vxr * w0yy + sqrt(theta_r) * (w1x - w0x)/3.0
      f8  = vxr * w1x + sqrt(theta_r) * (w1 - 4.0 * w0xx/5.0)
      f9  = vxr * w1y - 4.0 * sqrt(theta_r) * w0xy / 5.0
  
    else
  
      f1  = vyr * w0 + sqrt(theta_r) * w0y
      f2  = vyr * w0x + 2.0 * sqrt(theta_r) * w0xy
      f3  = vyr * w0y + sqrt(theta_r) * (w0 - w1 + 2.0 * w0yy)
      f4  = vyr * w1 + sqrt(theta_r) * (5.0 * w1y - 2.0 * w0y)/3.0
      f5  = vyr * w0xx + sqrt(theta_r) * (w1y - w0y)/3.0
      f6  = vyr * w0xy + sqrt(theta_r)  * (w0x - w1x)/2.0
      f7  = vyr * w0yy + sqrt(theta_r) * 2.0 * (w0y - w1y)/3.0
      f8  = vyr * w1x - sqrt(theta_r) * 4.0 * w0xy / 5.0
      f9  = vyr * w1y + sqrt(theta_r) * (w1 - 4.0 * w0yy / 5.0)
    end
    
    return SVector(f1, f2, f3, f4, f5, f6, f7, f8, f9)
  end
  
  @inline function cons2prim(u, equations::PerturbationMomentSystem2D)
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  
    return SVector(w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y)
  end
  
  
  @inline function cons2entropy(u, equations::PerturbationMomentSystem2D)
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  
    return SVector( w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y)
  end
  
  @inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::PerturbationMomentSystem2D)
      return 2.0
  end
  
  @inline function max_abs_speeds(u, equations::PerturbationMomentSystem2D)
    return 2.0, 2.0
  end
  
  #source function for all source terms to be 0:

  # @inline function source_terms_convergence_test(u, x, t, equations::PerturbationMomentSystem2D)
  #   w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  #   tau  = calc_tau()
  
  #   du1 = -w0 * w0xx + w0x * w0x/ 3.0 - w0y * w0y/ 6.0 
  #   du2 = -w0 * w0xy + w0x * w0y/ 4.0 + w0y * w0x/ 4.0
  #   du3 = -w0 * w0yy - w0x * w0x/ 6.0 + w0y * w0y/ 3.0
  #   du4 = 2.0 * w1 * w0x / 3.0 + 4.0 * w0x * w0xx/15.0 + 4.0 * w0y * w0xy/ 15.0 - 2.0 * w0 * w1x / 3.0
  #   du5 = 2.0 * w1 * w0y / 3.0  + 4.0 * w0y * w0yy/ 15.0 + 4.0 * w0x * w0xy/15.0 - 2.0 * w0 * w1y / 3.0
    
  #   #println(SVector(0.0, 0.0, 0.0, 0.0, du1/tau, du2/tau, du3/tau, du4/tau, du5/tau))
  
  #   return SVector(0.0, 0.0, 0.0, 0.0, du1/tau, du2/tau, du3/tau, du4/tau, du5/tau)
  # end


  @inline function source_terms_convergence_test(u, x, t, equations::PerturbationMomentSystem2D)
    
  @unpack vxr, vyr, theta_r = equations

    stheta = sqrt(theta_r)
    tau = calc_tau()
    v = 1
    c = 2
    A = 0.1
    L = 2
    J = 3
    K = 6
    D = 8
    f = 1/L
    g = 5/K
    j = J/5
    d = D/15

    h = 1/tau
    ω = 2 * pi * f

    ini1 = c + A * sin(ω * (x[1] + x[2] - t))
    ini2 = ini1^2

    tmp1 = ω * A * cos(ω * (x[1] + x[2] - t))
    tmp2 = vxr + vyr - v

    du1 = tmp1 * (tmp2 + c * stheta)
    du2 = tmp1 * (tmp2 + c * c * stheta)
    #du3 = tmp1 * (tmp2 + c * c * stheta)
    #du4 = tmp1 * (tmp2 + c * stehta)
    du5 = tmp1 * tmp2 + h * g * ini2
    du6 = tmp1 * tmp2 + h * f * ini2
    #du7 = tmp1 * tmp2 + h * g * ini2
    du8 = tmp1 * (tmp2 - j) - h * d * ini2
    #du9 = tmp1 * (tmp2 - j) - h * d * ini2


    return SVector(du1, du2, du2, du1, du5, du6, du5, du8, du8)
  end 
  
  function initial_condition_convergence_test(x, t, equations::PerturbationMomentSystem2D)
    c = 2
    A = 0.1
    L = 2
    K = 3
    f = 1/L
    g = 1/K
    ω = 2 * pi * f
    ini = c + A * sin(ω * (x[1] + x[2] - t))

    w0 = ini
    w0x = ini
    w0y = ini
    w1 = ini
    w0xx = ini
    w0xy = ini
    w0yy = ini
    w1x = ini
    w1y = ini
  
    return SVector(w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y)
  end
  
  function calc_tau()
    Pr = 2/3
    # gas constant for Xenon
    R = 63.327
    cp = 5/2 * R
    # thermal conductivity for Xenon
    lambda = 0.0055
    
    mu = Pr * lambda/cp
    p = 41.3

    tau = mu/p
    return(tau)

  end
  
  function initial_condition_constant(x, t, equations::PerturbationMomentSystem2D)
  
  
    return SVector(init_w0(x, t, equations::PerturbationMomentSystem2D), 
                   init_w0x(x, t, equations::PerturbationMomentSystem2D), 
                   init_w0y(x, t, equations::PerturbationMomentSystem2D), 
                   init_w1(x, t, equations::PerturbationMomentSystem2D),
                   init_w0xx(x, t, equations::PerturbationMomentSystem2D),
                   init_w0yy(x, t, equations::PerturbationMomentSystem2D),
                   init_w0xy(x, t, equations::PerturbationMomentSystem2D), 
                   init_w1x(x, t, equations::PerturbationMomentSystem2D),
                   init_w1y(x, t, equations::PerturbationMomentSystem2D))
  end
  
  
  @inline function init_w0(x, t, equations::PerturbationMomentSystem2D)
   
    rho_r = 1
    if (x[1] && x[2] < 0)
      drho = 3 - rho_r
    else
      drho = 1 - rho_r
    end
    
   
    w0 = 1.0 + drho / rho_r
    
   return w0
  end
  
  
  @inline function init_w0x(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    drho = 1
    rho_r = 2
    dv_x = 4
    
    w0x = dv_x / sqrt(theta_r) + (drho * dv_x)/(rho_r * sqrt(theta_r))
    
   return w0x
  end
  
  @inline function init_w0y(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    drho = 1 
    rho_r = 2
    dv_y = 4
    
    w0y = dv_y / sqrt(theta_r) + (drho * dv_y)/(rho_r * sqrt(theta_r))
    
   return w0y
  end
  
  @inline function init_w1(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    drho = 1
    rho_r = 2
    dv_y = 4
    dv_x = 4
    dtheta = -14.4
    rho = 3 
    w1 = - (rho * (dv_x *dv_x + dv_y * dv_y) )/(3.0 * (rho_r * theta_r)) - (drho * dtheta)/(rho_r * theta_r) - dtheta / theta_r
    
   return w1
  end
  
  
  @inline function init_w0xx(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    drho = 1
    rho_r = 2
    dv_y = 4
    dv_x = 4
    sigma_xx = 0
    
    w0xx = 0.5 * sigma_xx/(rho_r * theta_r) + (2.0 * dv_x * dv_x - dv_y * dv_y)/(6.0 * theta_r) + (2.0 * drho * dv_x * dv_x - drho * dv_y * dv_y)/(6.0 * rho_r * theta_r)  
    
   return w0xx
  end
  
  @inline function init_w0yy(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    drho = 1
    rho_r = 2
    dv_y = 4
    dv_x = 4
    sigma_yy = 0
    
    w0yy = 0.5 * sigma_yy/(rho_r * theta_r) + (2.0 * dv_y * dv_y - dv_x * dv_x)/(6.0 * theta_r) + (2.0 * drho * dv_y * dv_y - drho * dv_x * dv_x)/(6.0 * rho_r * theta_r)
    
   return w0yy
  end
  
  
  @inline function init_w0xy(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 2
    drho = 1
    sigma_xy = 0
    dv_y = 4
    dv_x = 4
    
    w0xy = 0.5 * sigma_xy/(rho_r * theta_r) + 0.5 * ((rho_r + drho) * (dv_x * dv_y))/(rho_r * theta_r)
    
   return w0xy
  end
  
  
  
  @inline function init_w1x(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    dtheta = -14.4
    rho_r = 2
    sigma_xy = 0
    sigma_xx = 0
    q_x = 138
    dv_y = 4
    dv_x = 4
    rho = 1

    w1x = - 2.0 * q_x / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xx * dv_x + sigma_xy * dv_y))/(5.0*rho_r* sqrt(theta_r).^3.0) - (dtheta * dv_x * rho)/ (rho_r * sqrt(theta_r).^3.0) - rho * (dv_x * dv_y^2 + dv_x^3)/(5.0 * rho_r * sqrt(theta_r).^3.0)
    
   return w1x
  end
  
  
  
  @inline function init_w1y(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    dtheta = -14.4
    rho_r = 2
    sigma_xy = 0
    sigma_yy = 0
    q_y = 138
    dv_x = 4
    dv_y = 4
    rho = 1
  
    w1y = - 2.0 * q_y / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xy * dv_x + sigma_yy * dv_y))/(5.0 * rho_r * sqrt(theta_r).^3.0) - (dtheta * dv_y * rho)/ (rho_r * sqrt(theta_r).^3.0) - rho * (dv_y * dv_x^2 + dv_y^3)/(5.0 * rho_r * sqrt(theta_r).^3.0) 
    
   return w1y
  end
  