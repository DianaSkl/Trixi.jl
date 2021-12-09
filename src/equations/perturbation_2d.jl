struct PerturbationMomentSystem2D{RealT<:Real} <: AbstractPerturbationMomentSystem{2, 9} 
    vxr::RealT
    vyr::RealT
    theta_r::RealT
    tau::RealT
  end
  varnames(::typeof(cons2prim), ::PerturbationMomentSystem2D) = ("rho", "vx", "vy", "theta")
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
      f6  = vxr * w0yy + sqrt(theta_r) * (w1x - w0x)/3.0
      f7  = vxr * w0xy + sqrt(theta_r) * (w0y - w1y)/2.0
      f8  = vxr * w1x + sqrt(theta_r) * (w1 - 4.0 * w0xx/5.0)
      f9  = vxr * w1y - 4.0 * sqrt(theta_r) * w0xy / 5.0
  
    else
  
      f1  = vyr * w0 + sqrt(theta_r) * w0y
      f2  = vyr * w0x + 2.0 * sqrt(theta_r) * w0xy
      f3  = vyr * w0y + sqrt(theta_r) * (w0 - w1 + 2.0 * w0yy)
      f4  = vyr * w1 + sqrt(theta_r) * (5.0 * w1y - 2.0 * w0y)/3.0
      f5  = vyr * w0xx + sqrt(theta_r) * (w1y - w0y)/3.0
      f6  = vyr * w0yy + sqrt(theta_r) * 2.0 * (w0y - w1y)/3.0
      f7  = vyr * w0xy + sqrt(theta_r)  * (w0x - w1x)/2.0
      f8  = vyr * w1x - sqrt(theta_r) * 4.0 * w0xy / 5.0
      f9  = vyr * w1y + sqrt(theta_r) * (w1 - 4.0 * w0yy / 5.0)

    end
    
    return SVector(f1, f2, f3, f4, f5, f6, f7, f8, f9)
  end
  
  @inline function cons2prim(u, equations::PerturbationMomentSystem2D)
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
    @unpack vxr, vyr, theta_r = equations
    rho_r = 2.0


    rho = w0 * rho_r
    v_x = vxr + w0x * sqrt(theta_r) / w0
    v_y = vyr + w0y * sqrt(theta_r) / w0
    theta = theta_r - (w0x^2 * theta_r)/(3 * w0^2) - (w0y^2 * theta_r)/(3 * w0^2) - (w1 * theta_r)/ w0
    # sigma_xx = - (2 * w0x^2 * theta_r * rho_r)/(3 * w0) + 2*w0xx * theta_r * rho_r + (w0y^2 * theta_r * rho_r)/(3 * w0)
    # sigma_xy = 2 * w0xy * theta_r * rho_r - (w0x * w0y * theta_r * rho_r)/w0
    # sigma_yy = (w0x^2 * theta_r * rho_r)/(3 * w0) - (2*w0y^2*theta_r*rho_r)/(3*w0)+2*w0yy*theta_r*rho_r
    # q_x = (w0x^3 * sqrt(theta_r)^3 * rho_r)/ w0^2 - (2*w0x*w0xx*sqrt(theta_r)^3*rho_r)/w0 - (2*w0xy*w0y*sqrt(theta_r)^3*rho_r)/w0 + (w0x*w0y^2*sqrt(theta_r)^3 * rho_r)/w0^2 + (5 * w0x *w1 *sqrt(theta_r)^3*rho_r)/(2*w0) - (5*w1x*sqrt(theta_r)^3*rho_r)/2
    # q_y = - (2*w0x * w0xy * sqrt(theta_r)^3*rho_r)/w0 + (w0x^2 * w0y * sqrt(theta_r)^3 * rho_r)/w0^2 + (w0y^3*sqrt(theta_r)^3 * rho_r)/w0^2 - (2*w0y *w0yy* sqrt(theta_r)^3 * rho_r)/w0 + (5*w0y*w1*sqrt(theta_r)^3*rho_r)/2*w0 - (5*w1y*sqrt(theta_r)^3*rho_r)/2

    
    return SVector(rho, v_x, v_y, theta)
  
 
  end


  
  @inline function cons2entropy(u, equations::PerturbationMomentSystem2D)
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  
    return SVector( w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y)
  end
  
  @inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::PerturbationMomentSystem2D)
    @unpack vxr, theta_r = equations

    ab1 = abs(vxr) + 2.0 * sqrt(5 * theta_r / 3)

    return ab1
 
  end
  
  @inline function max_abs_speeds(u, equations::PerturbationMomentSystem2D)
    @unpack vxr, vyr, theta_r = equations

    ab1 = abs(vxr) + 2.0 * sqrt(5 * theta_r / 3)
    ab2 = abs(vyr) + 2.0 * sqrt(5 * theta_r / 3)

    return SVector(ab1, ab2)
    
  end
  
  @inline function density_pressure(u, equations::PerturbationMomentSystem2D)
    rho = 1.0
    return rho
   end

  #source function for all source terms to be 0:

  # @inline function source_terms_convergence_test(u, x, t, equations::PerturbationMomentSystem2D)
  #   w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  #   tau  = calc_tau(equations)
  
  #   du1 = -w0 * w0xx + w0x * w0x/ 3.0 - w0y * w0y/ 6.0 
  #   du2 = -w0 * w0yy - w0x * w0x/ 6.0 + w0y * w0y/ 3.0
  #   du3 = -w0 * w0xy + w0x * w0y/ 4.0 + w0y * w0x/ 4.0
  #   du4 = 2.0 * w1 * w0x / 3.0 + 4.0 * w0x * w0xx/15.0 + 4.0 * w0y * w0xy/ 15.0 - 2.0 * w0 * w1x / 3.0
  #   du5 = 2.0 * w1 * w0y / 3.0 + 4.0 * w0y * w0yy/ 15.0 + 4.0 * w0x * w0xy/15.0 - 2.0 * w0 * w1y / 3.0
  
  
  #   return SVector(0.0, 0.0, 0.0, 0.0, du1/tau, du2/tau, du3/tau, du4/tau, du5/tau)
  # end


  function initial_condition_convergence_test(x, t, equations::PerturbationMomentSystem2D)
    c = 2
    A = 0.1
    L = 2
    f = 1/L
    ω = 2 * pi * f
    ini = c + A * sin(ω * (x[1] + x[2] - t))
 
    @unpack vxr,vyr, theta_r = equations

    rho = ini
    rho_r = 2
    drho = rho-rho_r
    vx = 1
    vy = 1
    dvx = vx-vxr
    dvy = vy-vyr
    theta = 2*(ini-1)/3
    dtheta = theta -  theta_r
    sxx = (ini-2)/100000
    syy = (ini-2)/100000
    sxy = (ini-2)/100000
    qx = (ini-2)/100000
    qy = (ini-2)/100000


    w0 = rho/rho_r
    w0x = (rho * dvx)/(rho_r * sqrt(theta_r))
    w0y = (rho * dvy)/(rho_r * sqrt(theta_r))
    w1 = -dtheta*rho/(rho_r*theta_r) - (rho*(dvx^2 + dvy^2))/(3*rho_r*theta_r)
    w0xx = (sxx+drho*(dvx^2*2/3 - dvy^2*1/3))/(2*rho_r*theta_r)+(dvx^2*2/3 - dvy^2*1/3)/(2*theta_r)
    w0yy = syy/(2*rho_r*theta_r)+ (rho*(dvy^2 * 2/3 - dvx^2*1/3))/(2*theta_r*rho_r)
    w0xy = sxy/(2*rho_r*theta_r) + (dvx*dvy/(2*theta_r)) + drho*dvx*dvy/(2*rho_r*theta_r)
    w1x =  -2*(qx+ sxx*dvx + sxy*dvy )/(5*rho_r*sqrt(theta_r^3)) - (dtheta*dvx*rho)/(rho_r*sqrt(theta_r^3))- (rho*(dvx^3+dvx*dvy^2))/(5*sqrt(theta_r^3)*rho_r)
    w1y =  -2*(qy+ sxy*dvx + syy*dvy )/(5*rho_r*sqrt(theta_r^3)) - (dtheta*dvy*rho)/(rho_r*sqrt(theta_r^3))- (rho*(dvy*dvx^2+dvy^3))/(5*sqrt(theta_r^3)*rho_r)
    return SVector(w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y)
  end
  


  @inline function source_terms_convergence_test(u, x, t, equations::PerturbationMomentSystem2D)
    
  @unpack vxr, vyr, theta_r, tau = equations

    c = 2
    A = 0.1
    L = 2
    f = 1/L
    ω = 2 * pi * f
    z = t
    x1 = x[1]
    x2 = x[2]

    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u

    a1 = -w0 * w0xx + w0x * w0x/ 3.0 - w0y * w0y/ 6.0 
    a2 = -w0 * w0yy - w0x * w0x/ 6.0 + w0y * w0y/ 3.0
    a3 = -w0 * w0xy + w0x * w0y/ 4.0 + w0y * w0x/ 4.0
    a4 = 2.0 * w1 * w0x / 3.0 + 4.0 * w0x * w0xx/15.0 + 4.0 * w0y * w0xy/ 15.0 - 2.0 * w0 * w1x / 3.0
    a5 = 2.0 * w1 * w0y / 3.0 + 4.0 * w0y * w0yy/ 15.0 + 4.0 * w0x * w0xy/15.0 - 2.0 * w0 * w1y / 3.0

    # a1 = 0
    # a2 = 0
    # a3 = 0
    # a4 = 0
    # a5 = 0
 
   du1 = (ω*cos((-z+x2+x1)*ω))/20
   du2 = (20000*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+435003*ω*cos((-z+x2+x1)*ω))/(1000000*sqrt(3))
   du3 = (20000*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+435003*ω*cos((-z+x2+x1)*ω))/(1000000*sqrt(3))
   du4 = -(55000*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+863757*ω*cos((-z+x2+x1)*ω))/1250000
   du5 = (150000*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+2209499*ω*cos((-z+x2+x1)*ω))/25000000 + a1/tau
   du6 = (150000*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+2209499*ω*cos((-z+x2+x1)*ω))/25000000+ a2/tau
   du7 = (450000*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+6628527*ω*cos((-z+x2+x1)*ω))/25000000 + a3/tau
   du8 = (1450000*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+15573677*ω*cos((-z+x2+x1)*ω))/(62500000*sqrt(3)) + a4/tau
   du9 = (1450000*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+15573677*ω*cos((-z+x2+x1)*ω))/(62500000*sqrt(3)) + a5/tau


    return SVector(du1, du2, du3, du4, du5, du6, du7, du8, du9)
  end 
  
 
  function calc_tau(equations::PerturbationMomentSystem2D)
    @unpack tau = equations
    return(tau)

  end


  function shocktube(x, equations::PerturbationMomentSystem2D)
    rho_r = 1

    if (x[1] < -20 || (x[1]<0 && x[1]>=-10) || (x[1]>=10 && x[1]<20))
      drho = 3 - rho_r
    elseif ((x[1] < -10 && x[1] >= -20) || (x[1]>=0 && x[1] <10) || x[1] >= 20)
      drho = 1 - rho_r
    end

    return drho
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
    drho = shocktube(x, equations)
  
   
    w0 = 1 + drho / rho_r
    
   return w0
  end
  
  
  @inline function init_w0x(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 

    rho_r = 1
    drho = shocktube(x, equations)

    dv_x = 0
    
    w0x = dv_x / sqrt(theta_r) + (drho * dv_x)/(rho_r * sqrt(theta_r))

   return w0x
  end
  
  @inline function init_w0y(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    
    w0y = dv_y / sqrt(theta_r) + (drho * dv_y)/(rho_r * sqrt(theta_r))
    
   return w0y
  end
  
  @inline function init_w1(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0
    dtheta = 0
    rho = 1 
    w1 = - (rho * (dv_x *dv_x + dv_y * dv_y) )/(3.0 * (rho_r * theta_r)) - (drho * dtheta)/(rho_r * theta_r) - dtheta / theta_r
    
   return w1
  end
  
  
  @inline function init_w0xx(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 

    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0

    sigma_xx = 0
    
    w0xx = 0.5 * sigma_xx/(rho_r * theta_r) + (2.0 * dv_x * dv_x - dv_y * dv_y)/(6.0 * theta_r) + (2.0 * drho * dv_x * dv_x - drho * dv_y * dv_y)/(6.0 * rho_r * theta_r)  
    
   return w0xx
  end
  
  @inline function init_w0yy(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0
    sigma_yy = 0
    
    w0yy = 0.5 * sigma_yy/(rho_r * theta_r) + (2.0 * dv_y * dv_y - dv_x * dv_x)/(6.0 * theta_r) + (2.0 * drho * dv_y * dv_y - drho * dv_x * dv_x)/(6.0 * rho_r * theta_r)
    
   return w0yy
  end
  
  
  @inline function init_w0xy(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0
    sigma_xy = 0
 
    
    w0xy = 0.5 * sigma_xy/(rho_r * theta_r) + 0.5 * ((rho_r + drho) * (dv_x * dv_y))/(rho_r * theta_r)
    
   return w0xy
  end
  
  
  
  @inline function init_w1x(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0
    sigma_xy = 0
    sigma_xx = 0
    q_x = 0

    dtheta = 0
    rho = 1

    w1x = - 2.0 * q_x / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xx * dv_x + sigma_xy * dv_y))/(5.0*rho_r* sqrt(theta_r).^3.0) - (dtheta * dv_x * rho)/ (rho_r * sqrt(theta_r).^3.0) - rho * (dv_x * dv_y^2 + dv_x^3)/(5.0 * rho_r * sqrt(theta_r).^3.0)
    
   return w1x
  end
  
  
  
  @inline function init_w1y(x, t, equations::PerturbationMomentSystem2D)
    
    @unpack theta_r = equations 
    dtheta = 0
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0
    sigma_xy = 0
    sigma_yy = 0
    q_y = 0

    rho = 1
  
    w1y = - 2.0 * q_y / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xy * dv_x + sigma_yy * dv_y))/(5.0 * rho_r * sqrt(theta_r).^3.0) - (dtheta * dv_y * rho)/ (rho_r * sqrt(theta_r).^3.0) - rho * (dv_y * dv_x^2 + dv_y^3)/(5.0 * rho_r * sqrt(theta_r).^3.0) 
    
   return w1y
  end