struct MomentSystem1D{RealT<:Real} <: AbstractMomentSystem{1, 5} 
    vxr::RealT
    theta_r::RealT
    rho_r::RealT
    tau::RealT
  end

  
varnames(::typeof(cons2prim), ::MomentSystem1D) = ("ρ", "v\u2093", "p")
varnames(::typeof(cons2cons), ::MomentSystem1D) = ("w\u207D\u2070\u207E", "w\u207D\u2070\u207E\u2093", "w\u207D\u00B9\u207E", "w\u207D\u2070\u207E\u2093\u2093", "w\u207D\u00B9\u207E\u2093" )

  
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

  @inline function flux_shima_etal(u_ll, u_rr, orientation::Integer, equations::MomentSystem1D)
    @unpack vxr, theta_r, rho_r = equations
    
    # Unpack left and right state
    rho_ll, vx_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, vx_rr, p_rr = cons2prim(u_rr, equations)
  
    # Average each factor of products in flux
    rho_avg = 1/2 * (rho_ll + rho_rr)
    vx_avg  = 1/2 * ( vx_ll +  vx_rr)
    
    theta_avg = 0.5 * (p_ll/rho_ll + p_rr/rho_rr)
    kin_avg = 1/2 * (vx_ll*vx_rr)
  
    dvx = vx_avg - vxr 
    dvxs = kin_avg 
    drho = rho_avg - rho_r
    dtheta = theta_avg - theta_r 
    sxx = 0.0
    qx = 0.0
  
    w0 = rho_avg/rho_r
    w0x = (rho_avg * dvx)/(rho_r * sqrt(theta_r))
    w1 = -dtheta*rho_avg/(rho_r*theta_r) - (rho_avg*(dvx^2))/(3*rho_r*theta_r)
    w0xx = (sxx+drho*(dvx^2*2/3))/(2*rho_r*theta_r)+(dvx^2*2/3)/(2*theta_r)
    w1x =  -2*(qx+ sxx*dvx)/(5*rho_r*sqrt(theta_r^3)) - (dtheta*dvx*rho_avg)/(rho_r*sqrt(theta_r^3))- (rho_avg*(dvxs*dvx))/(5*sqrt(theta_r^3)*rho_r)
 
    W = SVector(w0, w0x, w1, w0xx, w1x)      
    f1, f2, f3, f4, f5 = flux(W, orientation, equations)
         
    return SVector(f1, f2, f3, f4, f5)
  end
  

  @inline function cons2prim(prim, equations::MomentSystem1D)
    w0, w0x, w1, w0xx, w1x = prim
    @unpack vxr, theta_r, rho_r = equations
   
    rho = w0 * rho_r
    vx = vxr + w0x * sqrt(theta_r) / w0
    theta = theta_r - (w0x^2 * theta_r)/(3 * w0^2) -  (w1 * theta_r)/ w0
    p = rho*theta
    sigmax = - (2 * w0x^2 * theta_r * rho_r)/(3 * w0) + 2*w0xx * theta_r * rho_r
    qx = (w0x^3 * sqrt(theta_r)^3 * rho_r)/ w0^2 - (2*w0x*w0xx*sqrt(theta_r)^3*rho_r)/w0  + (5 * w0x *w1 *sqrt(theta_r)^3*rho_r)/(2*w0) - (5*w1x*sqrt(theta_r)^3*rho_r)/2
    return SVector(rho, vx, p)
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
  
  @inline function density_pressure(u, equations::MomentSystem1D)
    rho = 1.0
    return rho
  end
   
  function initial_condition_convergence_test(x, t, equations::MomentSystem1D)
    
    @unpack vxr, theta_r, rho_r = equations 
    c = 2
    A = 0.1
    L = 2
    f = 1/L
    ω = 2 * pi * f
  
    ini = c + A * sin(ω * (x[1] - t))
    #ini = c + A * sin(ω * (- t))
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
        
   #zeit und ort
   #f1,f2,f3,f4,f5 = source_spatialtime(u, x, t, equations::MomentSystem1D)

   #nur Zeit
   #f1,f2,f3,f4,f5 = source_time(u, x, t, equations::MomentSystem1D)

   #produktionen
   f1,f2,f3,f4,f5 = source_productions(u, x, t, equations::MomentSystem1D)
   

    return(f1,f2,f3,f4,f5)
  end 


  @inline function source_productions(u, x, t, equations::MomentSystem1D)

    @unpack vxr, theta_r, rho_r, tau = equations 
    w0, w0x, w1, w0xx, w1x = u
  
    a1 = -w0 * w0xx + w0x * w0x/3.0
    a2 = -2.0 * w0 * w1x/3.0 + 2.0 * w1 * w0x/3.0 + 4.0 * w0x * w0xx/15.0
  
    #Alternativ wie im Paper:
    # sigma_xx = - (2 * w0x^2 * theta_r * rho_r)/(3 * w0) + 2*w0xx * theta_r * rho_r
    # q_x = (w0x^3 * sqrt(theta_r)^3 * rho_r)/ w0^2 - (2*w0x*w0xx*sqrt(theta_r)^3*rho_r)/w0   + (5 * w0x *w1 *sqrt(theta_r)^3*rho_r)/(2*w0) - (5*w1x*sqrt(theta_r)^3*rho_r)/2
    # vx = vxr + w0x * sqrt(theta_r) / w0
    # dvx = vx - vxr 

    # a1 = - 0.5*sigma_xx/(rho_r * theta_r)
    # a2 = - 4*q_x/(15*rho_r*sqrt(theta_r^3)) + 2*sigma_xx*dvx/(5*rho_r * theta_r) 

    return (0,0,0,a1/tau,a2/tau)
  end 



  @inline function source_spatialtime(u, x, t, equations::MomentSystem1D)
    
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

   return (f1,f2,f3,f4,f5)
  end

  @inline function source_time(u, x, t, equations::MomentSystem1D)

    @unpack vxr, theta_r, rho_r, tau = equations 

    c = 2
    A = 0.1
    L = 2
    f = 1/L
    ω = 2 * pi * f
      
    p1 = 2/3
    p2 = 1/3

    vx = 1.0
    dv_x = vx - vxr

    dw0_t = -(A*ω*cos((-t)*ω))/c
    
    dw0x_t = -(A*dv_x*ω*cos((-t)*ω))/(c*sqrt(theta_r))

    dw1_t =-(-(A*ω*cos((-t)*ω)*(p1*(A*sin((-t)*ω)+c)-theta_r-p2))/(c*theta_r)-(A*p1*ω*cos((-t)*ω)*(A*sin((-t)*ω)+c))/(c*theta_r)-(A*dv_x^2*ω*cos((-t)*ω))/(3*c*theta_r))
    
    dw0xx_t = -(A*dv_x^2*p1*ω*cos((-t)*ω))/(c^2*theta_r)

    dw1x_t =-(-(A*dv_x*ω*cos((-t)*ω)*(p1*(A*sin((-t)*ω)+c)-theta_r-p2))/(c*theta_r^1.5)-(A*dv_x*p1*ω*cos((-t)*ω)*(A*sin((-t)*ω)+c))/(c*theta_r^1.5)-(A*dv_x^3*ω*cos((-t)*ω))/(5*c*theta_r^1.5))
    

    f1  = dw0_t 
    f2  = dw0x_t 
    f3  = dw1_t 
    f4  = dw0xx_t 
    f5  = dw1x_t 


    return (f1,f2,f3,f4,f5)
  end

  function initial_condition_constant(x, t, equations::MomentSystem1D)
  
    return SVector(init_w0(x, t, equations::MomentSystem1D), 
                   init_w0x(x, t, equations::MomentSystem1D), 
                   init_w1(x, t, equations::MomentSystem1D),
                   init_w0xx(x, t, equations::MomentSystem1D),
                   init_w1x(x, t, equations::MomentSystem1D))
  end

  function shocktube_density(x, equations::MomentSystem1D)
    @unpack rho_r = equations 
    if (x[1] < 0)
      drho = 3 - rho_r
    else
      drho = 1 - rho_r
    end
  
    return drho
    end

    function shocktube_temp(x, equations::MomentSystem1D)
      @unpack theta_r = equations 
      if (x[1] < 0)
        dtheta = 1.0 - theta_r 
      else
        dtheta = 1.0 - theta_r
      end
    
      return dtheta
      end
  
  
  @inline function init_w0(x, t, equations::MomentSystem1D)
    @unpack rho_r = equations 
  
    drho = shocktube_density(x, equations)
    w0 = 1 + drho / rho_r
   return w0
  end
  
  
  @inline function init_w0x(x, t, equations::MomentSystem1D)
    
    @unpack theta_r, rho_r = equations 
  
 
    drho = shocktube_density(x, equations)
  
    dv_x = 0
    
    w0x = dv_x / sqrt(theta_r) + (drho * dv_x)/(rho_r * sqrt(theta_r))
  
   return w0x
  end
  

  
  @inline function init_w1(x, t, equations::MomentSystem1D)
    
    @unpack theta_r, rho_r = equations 
    drho = shocktube_density(x, equations)
    dv_y = 0
    dv_x = 0
    dtheta = shocktube_temp(x, equations)
    rho = 1 
    w1 = - (rho * (dv_x *dv_x + dv_y * dv_y) )/(3.0 * (rho_r * theta_r)) - (drho * dtheta)/(rho_r * theta_r) - dtheta / theta_r
    
   return w1
  end
  
  
  @inline function init_w0xx(x, t, equations::MomentSystem1D)
    
    @unpack theta_r, rho_r = equations 
  
    drho = shocktube_density(x, equations)
    dv_y = 0
    dv_x = 0
  
    sigma_xx = 0
    
    w0xx = 0.5 * sigma_xx/(rho_r * theta_r) + (2.0 * dv_x * dv_x - dv_y * dv_y)/(6.0 * theta_r) + (2.0 * drho * dv_x * dv_x - drho * dv_y * dv_y)/(6.0 * rho_r * theta_r)  
    
   return w0xx
  end
    
   
  
  @inline function init_w1x(x, t, equations::MomentSystem1D)
    
    @unpack theta_r, rho_r = equations 
  
    dv_y = 0
    dv_x = 0
    sigma_xy = 0
    sigma_xx = 0
    q_x = 0
  
    dtheta = shocktube_temp(x, equations)
    rho = 1
  
    w1x = - 2.0 * q_x / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xx * dv_x + sigma_xy * dv_y))/(5.0*rho_r* sqrt(theta_r).^3.0) - (dtheta * dv_x * rho)/ (rho_r * sqrt(theta_r).^3.0) - rho * (dv_x * dv_y^2 + dv_x^3)/(5.0 * rho_r * sqrt(theta_r).^3.0)
  
   return w1x
  end
  
  
