
struct MomentSystem1D{RealT<:Real} <: AbstractMomentSystem{1, 5} 
    vxr::RealT
    theta_r::RealT
    rho_r::RealT
    tau::RealT
  end

  
varnames(::typeof(cons2prim), ::MomentSystem1D) = ("ρ", "v\u2093", "p", "σ\u2093\u2093", "q\u2093")
#varnames(::typeof(cons2prim), ::MomentSystem1D) = ("ρ", "v\u2093", "p")
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

  
  @inline function flux_ds(u_ll, u_rr, orientation::Integer, equations::MomentSystem1D)

    @unpack vxr, theta_r, rho_r = equations

    # Unpack left and right state
    # rho_ll, vx_ll, p_ll, sx_ll, qx_ll = cons2prim(u_ll, equations)
    # rho_rr, vx_rr, p_rr, sx_rr, qx_rr = cons2prim(u_rr, equations)
    rho_ll, vx_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, vx_rr, p_rr = cons2prim(u_rr, equations)


    # Average each factor of products in flux
    rho_avg = 1/2 * (rho_ll + rho_rr)
    vx_avg  = 1/2 * (vx_ll +  vx_rr)
    p_avg = 1/2 * (p_ll + p_rr)
    # sx_avg = 1/2 * (sx_ll + sx_rr)
    # qx_avg = 1/2 * (qx_ll + qx_rr)

    #W = SVector(prim2cons((rho_avg, vx_avg, p_avg, sx_avg, qx_avg),equations)) 
    W = SVector(prim2cons((rho_avg, vx_avg, p_avg),equations)) 
    f1, f2, f3, f4, f5 = flux(W, orientation, equations)
         
    return SVector(f1, f2, f3, f4, f5)
  end
  

  @inline function cons2prim(prim, equations::MomentSystem1D)
    w0, w0x, w1, w0xx, w1x = prim
    @unpack vxr, theta_r, rho_r = equations
   
    s3theta = sqrt(theta_r^3)


    rho = w0 * rho_r
    vx = vxr + w0x * sqrt(theta_r) / w0
    dvx = vx - vxr
    theta = theta_r - (w0x^2 * theta_r)/(3 * w0^2) -  (w1 * theta_r)/ w0
    dtheta = theta - theta_r
    p = rho*theta
    sigmax = 2*rho_r*theta_r*w0xx - rho*(2*dvx^2)/3
    qx =  -w1x*rho_r*s3theta*5/2 - sigmax*dvx  - dtheta*dvx*rho*5/2 - 0.5*rho*(dvx^3)
 
    return SVector(rho, vx, p, sigmax, qx)
    #return SVector(rho, vx, p)
  end
  
  @inline function prim2cons(prim, equations::MomentSystem1D)
    #rho, vx, p, sxx, qx = prim
    rho, vx, p = prim
    @unpack vxr, theta_r, rho_r = equations
  
    sxx = qx = 0.0
  
    dvx = vx - vxr 
    theta = p/rho
    dtheta = theta - theta_r 
  

    w0 = rho/rho_r
    w0x = (rho* dvx)/(rho_r * sqrt(theta_r))
    w1 = -dtheta*rho/(rho_r*theta_r) - (rho*(dvx^2))/(3*rho_r*theta_r)
    w0xx = 0.5*sxx/(rho_r*theta_r) + 0.5*(rho*2*(dvx^2)/3)/(rho_r*theta_r)
    w1x =  -2*(qx + sxx*dvx)/(5*rho_r*sqrt(theta_r^3)) - (dtheta*dvx*rho)/(rho_r*sqrt(theta_r^3))- (rho*(dvx^3))/(5*sqrt(theta_r^3)*rho_r)
 
    return SVector(w0, w0x, w1, w0xx, w1x)
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
    a = cons2prim(u, equations)
    rho = a[1]
    p = a[3]
    rho_times_p = rho*p
    return rho_times_p
  end

  @inline function density(u, equations::MomentSystem1D)
    a = cons2prim(u, equations)
    rho = a[1]
    return rho
  end

  @inline function pressures(u, equations::MomentSystem1D)
    a = cons2prim(u, equations)
    p = a[3]
    return p
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
        
   #zeit und ort
   #f1,f2,f3,f4,f5 = source_convergence_spatialtime(u, x, t, equations::MomentSystem1D)

   #nur Zeit
   #f1,f2,f3,f4,f5 = source_time(u, x, t, equations::MomentSystem1D)

   #produktionen
   f1,f2,f3,f4,f5 = source_productions(u, equations::MomentSystem1D)
   

    return(f1,f2,f3,f4,f5)
  end 


  @inline function source_productions(u, equations::MomentSystem1D)

    @unpack tau = equations 
    w0, w0x, w1, w0xx, w1x = u
  
    p1 = -w0 * w0xx + w0x * w0x/3.0
    p2 = -2.0 * w0 * w1x/3.0 + 2.0 * w1 * w0x/3.0 + 4.0 * w0x * w0xx/15.0
  
    return (0,0,0,p1/tau,p2/tau)
  end 



  @inline function source_convergence_spatialtime(u, x, t, equations::MomentSystem1D)
    
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
    dtheta = shocktube_temp(x, equations)
    dv_x = 0
    rho = drho + rho_r 
    w1 = - (rho * (dv_x *dv_x) )/(3.0 * (rho_r * theta_r)) - (drho * dtheta)/(rho_r * theta_r) - dtheta / theta_r
    
   return w1
  end
  
  
  @inline function init_w0xx(x, t, equations::MomentSystem1D)
    
    @unpack theta_r, rho_r = equations 
  
    drho = shocktube_density(x, equations)
    rho = drho + rho_r
    dv_x = 0
  
    sigma_xx = 0.0
    
    w0xx = 0.5 * sigma_xx/(rho_r * theta_r) + (rho* dv_x * dv_x)/(3 * rho_r*theta_r) 
    
   return w0xx
  end
    
   
  
  @inline function init_w1x(x, t, equations::MomentSystem1D)
    
    @unpack theta_r, rho_r = equations 

    dv_x = 0
    sigma_xx = 0.0
    q_x = 0.0
  
    dtheta = shocktube_temp(x, equations)
    drho = shocktube_density(x, equations)
    rho = drho + rho_r
  
    w1x = - 2.0 * q_x / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xx * dv_x))/(5.0*rho_r* sqrt(theta_r).^3.0) - (dtheta * dv_x * rho)/ (rho_r * sqrt(theta_r).^3.0) - rho * (dv_x^3)/(5.0 * rho_r * sqrt(theta_r).^3.0)
  
   return w1x
  end
  
