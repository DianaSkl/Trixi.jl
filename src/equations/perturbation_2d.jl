struct PerturbationMomentSystem2D{RealT<:Real} <: AbstractPerturbationMomentSystem{2, 9} 
    vxr::RealT
    vyr::RealT
    theta_r::RealT
    rho_r::RealT
    tau::RealT
end

varnames(::typeof(cons2prim), ::PerturbationMomentSystem2D) = ("rho", "v1", "v2", "theta")
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
    @unpack vxr, vyr, theta_r, rho_r = equations


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
    rho, v_x, v_y, theta = cons2prim(u, equations)
    rho_times_p = theta * rho^2 
    return rho_times_p
  end
  
  @inline function density(u, equations::PerturbationMomentSystem2D)
    a = cons2prim(u, equations)
    rho = a[1]
    return rho
  end
   
  
  
  function initial_condition_convergence_test(x, t, equations::PerturbationMomentSystem2D)
    c = 2
    A = 0.1
    L = 2
    f = 1/L
    ω = 2 * pi * f
    ini = c + A * sin(ω * (x[1] + x[2] - t))
  
   
    @unpack vxr, vyr, theta_r, rho_r, tau = equations
  
    gamma = 5/3

    x1 = x[1]
    x2 = x[2]

    rho = ini
    drho = rho-rho_r
    vx = 1.0
    vy = 1.0
    dvx = vx-vxr
    dvy = vy-vyr
    theta = (gamma-1)*(ini-(vx^2 + vy^2)/2)
    dtheta = theta -  theta_r

    r1 = 100000000
    r2 = 100000000
  
    sxx = (ini-2)/r1
    syy = (ini-2)/r1
    sxy = (ini-2)/r1

    qx = (ini-2)/r2
    qy = (ini-2)/r2
    # qx = (15*A*tau*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c)*(A*sin((x2+x1-t)*ω)-0.5*(vy^2+vx^2)+c)*(gamma-1)^2)/4
    # qy = (15*A*tau*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c)*(A*sin((x2+x1-t)*ω)-0.5*(vy^2+vx^2)+c)*(gamma-1)^2)/4
  
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
    
  @unpack vxr, vyr, theta_r, rho_r, tau = equations
  
    c = 2
    A = 0.1
    L = 2
    f = 1/L
    m = 1/L
    ω = 2 * pi * f
    z = t
    x1 = x[1]
    x2 = x[2]
  
    vx = 1.0
    vy = 1.0

    p1 = 2/3
    p2 = 1/3

    r1 = 100000000
    r2 = 100000000
  
    gamma = 5/3
  
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  
    a1 = -w0 * w0xx + w0x * w0x/ 3.0 - w0y * w0y/ 6.0 
    a2 = -w0 * w0yy - w0x * w0x/ 6.0 + w0y * w0y/ 3.0
    a3 = -w0 * w0xy + w0x * w0y/ 4.0 + w0y * w0x/ 4.0
    a4 = 2.0 * w1 * w0x / 3.0 + 4.0 * w0x * w0xx/15.0 + 4.0 * w0y * w0xy/ 15.0 - 2.0 * w0 * w1x / 3.0
    a5 = 2.0 * w1 * w0y / 3.0 + 4.0 * w0y * w0yy/ 15.0 + 4.0 * w0x * w0xy/15.0 - 2.0 * w0 * w1y / 3.0
  

    # dw0_x = (A*ω*cos((x2+x1-t)*ω))/c

    # dw0x_x = (A*(vx-vxr)*ω*cos((x2+x1-t)*ω))/(c*sqrt(theta_r))

    # dw0y_x = (A*(vy-vyr)*ω*cos((x2+x1-t)*ω))/(c*sqrt(theta_r))
 
    # dw1_x = (A*ω*cos((x2+x1-t)*ω)*(1-p1*(A*sin((x2+x1-t)*ω)+c)))/(c*theta_r)-(A*p1*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c))/(c*theta_r)-(A*((vy-vyr)^2+(vx-vxr)^2)*ω*cos((x2+x1-t)*ω))/(3*c*theta_r)
    
    # dw0xx_x = (A*(p1*(vx-vxr)^2-p2*(vy-vyr)^2)*ω*cos((x2+x1-t)*ω)+(A*ω*cos((x2+x1-t)*ω))/r1)/(c^2*theta_r)
    
    # dw0yy_x = (A*(p1*(vy-vyr)^2-p2*(vx-vxr)^2)*ω*cos((x2+x1-t)*ω))/(c^2*theta_r)+(A*ω*cos((x2+x1-t)*ω))/(c^2*r1*theta_r)

    # dw0xy_x = (A*(vx-vxr)*(vy-vyr)*ω*cos((x2+x1-t)*ω))/(2*c*theta_r)+(A*ω*cos((x2+x1-t)*ω))/(c^2*r1*theta_r)

    # dw1x_x = -(-(15*A*tau*ω^2*sin((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c)*(A*sin((x2+x1-t)*ω)-0.5*(vy^2+vx^2)+c)*(gamma-1)^2)/4+(15*A^2*tau*ω^2*cos((x2+x1-t)*ω)^2*(A*sin((x2+x1-t)*ω)-0.5*(vy^2+vx^2)+c)*(gamma-1)^2)/4+(15*A^2*tau*ω^2*cos((x2+x1-t)*ω)^2*(A*sin((x2+x1-t)*ω)+c)*(gamma-1)^2)/4+(A*(vy-vyr)*ω*cos((x2+x1-t)*ω))/r1+(A*(vx-vxr)*ω*cos((x2+x1-t)*ω))/r1)/(5*theta_r^(3/2))-(A*(vx-vxr)*ω*cos((x2+x1-t)*ω)*(p1*(A*sin((x2+x1-t)*ω)+c)-1))/(c*theta_r^(3/2))-(A*p1*(vx-vxr)*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c))/(c*theta_r^(3/2))-(A*((vx-vxr)*(vy-vyr)^2+(vx-vxr)^3)*ω*cos((x2+x1-t)*ω))/(5*c*theta_r^(3/2))
  
    # dw1y_x = -(-(15*A*tau*ω^2*sin((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c)*(A*sin((x2+x1-t)*ω)-0.5*(vy^2+vx^2)+c)*(gamma-1)^2)/4+(15*A^2*tau*ω^2*cos((x2+x1-t)*ω)^2*(A*sin((x2+x1-t)*ω)-0.5*(vy^2+vx^2)+c)*(gamma-1)^2)/4+(15*A^2*tau*ω^2*cos((x2+x1-t)*ω)^2*(A*sin((x2+x1-t)*ω)+c)*(gamma-1)^2)/4+(A*(vy-vyr)*ω*cos((x2+x1-t)*ω))/r1+(A*(vx-vxr)*ω*cos((x2+x1-t)*ω))/r1)/(5*theta_r^(3/2))-(A*(vy-vyr)*ω*cos((x2+x1-t)*ω)*(p1*(A*sin((x2+x1-t)*ω)+c)-1))/(c*theta_r^(3/2))-(A*p1*(vy-vyr)*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c))/(c*theta_r^(3/2))-(A*((vx-vxr)^2*(vy-vyr)+(vy-vyr)^3)*ω*cos((x2+x1-t)*ω))/(5*c*theta_r^(3/2))
   



    dw0_x = (A*ω*cos((x2+x1-t)*ω))/rho_r

    dw0x_x = (A*(vx-vxr)*ω*cos((x2+x1-t)*ω))/(rho_r*sqrt(theta_r))

    dw0y_x = (A*(vy-vyr)*ω*cos((x2+x1-t)*ω))/(rho_r*sqrt(theta_r))
 
    dw1_x = -(A*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c)*(gamma-1))/(rho_r*theta_r)+(A*ω*cos((x2+x1-t)*ω)*(theta_r-(A*sin((x2+x1-t)*ω)-m*(vy^2+vx^2)+c)*(gamma-1)))/(rho_r*theta_r)-(A*((vy-vyr)^2+(vx-vxr)^2)*ω*cos((x2+x1-t)*ω))/(3*rho_r*theta_r)
    
    dw0xx_x = (A*(p1*(vx-vxr)^2-p2*(vy-vyr)^2)*ω*cos((x2+x1-t)*ω)+(A*ω*cos((x2+x1-t)*ω))/r1)/(c*rho_r*theta_r)
    
    dw0yy_x = (A*(p1*(vy-vyr)^2-p2*(vx-vxr)^2)*ω*cos((x2+x1-t)*ω))/(c*rho_r*theta_r)+(A*ω*cos((x2+x1-t)*ω))/(c*r1*rho_r*theta_r)
   
    dw0xy_x = (A*(vx-vxr)*(vy-vyr)*ω*cos((x2+x1-t)*ω))/(c*rho_r*theta_r)+(A*ω*cos((x2+x1-t)*ω))/(c*r1*rho_r*theta_r)
 
    dw1x_x = -(A*(vx-vxr)*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c)*(gamma-1))/(rho_r*theta_r^(3/2))-(A*(vx-vxr)*ω*cos((x2+x1-t)*ω)*((A*sin((x2+x1-t)*ω)-m*(vy^2+vx^2)+c)*(gamma-1)-theta_r))/(rho_r*theta_r^(3/2))-(c*((A*(vy-vyr)*ω*cos((x2+x1-t)*ω))/r1+(A*(vx-vxr)*ω*cos((x2+x1-t)*ω))/r1+(A*ω*cos((x2+x1-t)*ω))/r2))/(5*rho_r*theta_r^(3/2))-(A*((vx-vxr)*(vy-vyr)^2+(vx-vxr)^3)*ω*cos((x2+x1-t)*ω))/(5*rho_r*theta_r^(3/2))
  
    dw1y_x = -(A*(vy-vyr)*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c)*(gamma-1))/(rho_r*theta_r^(3/2))-(A*(vy-vyr)*ω*cos((x2+x1-t)*ω)*((A*sin((x2+x1-t)*ω)-m*(vy^2+vx^2)+c)*(gamma-1)-theta_r))/(rho_r*theta_r^(3/2))-(c*((A*(vy-vyr)*ω*cos((x2+x1-t)*ω))/r1+(A*(vx-vxr)*ω*cos((x2+x1-t)*ω))/r1+(A*ω*cos((x2+x1-t)*ω))/r2))/(5*rho_r*theta_r^(3/2))-(A*((vx-vxr)^2*(vy-vyr)+(vy-vyr)^3)*ω*cos((x2+x1-t)*ω))/(5*rho_r*theta_r^(3/2))
   
   f1 = -dw0_x + vxr * dw0_x  + sqrt(theta_r) * dw0x_x + vyr * dw0_x + sqrt(theta_r) * dw0y_x
   f2 = -dw0x_x + vxr * dw0x_x + sqrt(theta_r) * (dw0_x - dw1_x) + 2.0 * sqrt(theta_r) * dw0xx_x + vyr * dw0x_x + 2.0 * sqrt(theta_r) * dw0xy_x
   f3 = -dw0y_x + vxr * dw0y_x + 2.0 * sqrt(theta_r) * dw0xy_x + vyr * dw0y_x + sqrt(theta_r) * (dw0_x - dw1_x + 2.0 * dw0yy_x)
   f4 = -dw1_x + vxr * dw1_x + sqrt(theta_r) * (5.0 * dw1x_x - 2.0 * dw0x_x)/3.0 +vyr * dw1_x + sqrt(theta_r) * (5.0 * dw1y_x - 2.0 * dw0y_x)/3.0
   f5 = -dw0xx_x + vxr * dw0xx_x + 2.0 * sqrt(theta_r)* (dw0x_x - dw1x_x)/3.0  + vyr * dw0xx_x + sqrt(theta_r) * (dw1y_x - dw0y_x)/3.0 + a1/tau
   f6 = -dw0yy_x + vxr * dw0yy_x + sqrt(theta_r) * (dw1x_x - dw0x_x)/3.0 + vyr * dw0yy_x + sqrt(theta_r) * 2.0 * (dw0y_x - dw1y_x)/3.0 + a2/tau
   f7 = -dw0xy_x + vxr * dw0xy_x + sqrt(theta_r) * (dw0y_x - dw1y_x)/2.0 + vyr * dw0xy_x + sqrt(theta_r)  * (dw0x_x - dw1x_x)/2.0 + a3/tau
   f8 = -dw1x_x + vxr * dw1x_x + sqrt(theta_r) * (dw1_x - 4.0 * dw0xx_x/5.0) + vyr * dw1x_x - sqrt(theta_r) * 4.0 * dw0xy_x / 5.0 + a4/tau
   f9 = -dw1y_x + vxr * dw1y_x - 4.0 * sqrt(theta_r) * dw0xy_x / 5.0 + vyr * dw1y_x + sqrt(theta_r) * (dw1_x - 4.0 * dw0yy_x / 5.0) + a5/tau

  

   du1 = (ω*cos((-z+x2+x1)*ω))/20
   du2 = ((3*r1*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+45*r1*ω*cos((-z+x2+x1)*ω))*gamma-3*r1*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(-15*r1*vxr-30*r1+30)*ω*cos((-z+x2+x1)*ω))/(100*sqrt(3)*r1)
   du3 = ((3*r1*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+45*r1*ω*cos((-z+x2+x1)*ω))*gamma-3*r1*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(-15*r1*vyr-30*r1+30)*ω*cos((-z+x2+x1)*ω))/(100*sqrt(3)*r1)
   du4 = (((2*r1*r2*vyr+2*r1*r2*vxr-7*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(30*r1*r2*vyr+30*r1*r2*vxr-105*r1*r2)*ω*cos((-z+x2+x1)*ω))*gamma+(-2*r1*r2*vyr-2*r1*r2*vxr+7*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(-5*r1*r2*vyr^2+(20-20*r1)*r2*vyr-5*r1*r2*vxr^2+(20-20*r1)*r2*vxr+(100*r1-40)*r2-20*r1)*ω*cos((-z+x2+x1)*ω))/(100*r1*r2)
   du5 = (((2*r1*r2*vyr-4*r1*r2*vxr+2*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(30*r1*r2*vyr-60*r1*r2*vxr+30*r1*r2)*ω*cos((-z+x2+x1)*ω))*gamma+(-2*r1*r2*vyr+4*r1*r2*vxr-2*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(-3*r1*r2*vyr^3+(13*r1*r2-9*r1*r2*vxr)*vyr^2+(12*r1*r2*vxr^2-6*r1*r2*vxr+(11-35*r1)*r2)*vyr+6*r1*r2*vxr^3-20*r1*r2*vxr^2+(73*r1+11)*r2*vxr+(-31*r1-7)*r2+4*r1)*ω*cos((-z+x2+x1)*ω))/(200*r1*r2) + a1/tau
   du6 = -(((4*r1*r2*vyr-2*r1*r2*vxr-2*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(60*r1*r2*vyr-30*r1*r2*vxr-30*r1*r2)*ω*cos((-z+x2+x1)*ω))*gamma+(-4*r1*r2*vyr+2*r1*r2*vxr+2*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(-6*r1*r2*vyr^3+(20*r1*r2-12*r1*r2*vxr)*vyr^2+(9*r1*r2*vxr^2+6*r1*r2*vxr+(-73*r1-11)*r2)*vyr+3*r1*r2*vxr^3-13*r1*r2*vxr^2+(35*r1-11)*r2*vxr+(31*r1+7)*r2-4*r1)*ω*cos((-z+x2+x1)*ω))/(200*r1*r2) + a2/tau
   du7 = -(((3*r1*r2*vyr+3*r1*r2*vxr-6*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(45*r1*r2*vyr+45*r1*r2*vxr-90*r1*r2)*ω*cos((-z+x2+x1)*ω))*gamma+(-3*r1*r2*vyr-3*r1*r2*vxr+6*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(3*r1*r2*vyr^3+(3*r1*r2-12*r1*r2*vxr)*vyr^2+(-12*r1*r2*vxr^2+33*r1*r2*vxr+(-57*r1-3)*r2)*vyr+3*r1*r2*vxr^3+3*r1*r2*vxr^2+(-57*r1-3)*r2*vxr+(93*r1-9)*r2-12*r1)*ω*cos((-z+x2+x1)*ω))/(200*r1*r2) + a3/tau
   du8 = ((((9*r1*r2*vxr-9*r1*r2)*vyr+9*r1*r2*vxr^2-18*r1*r2*vxr+6*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+((135*r1*r2*vxr-135*r1*r2)*vyr+135*r1*r2*vxr^2-270*r1*r2*vxr+90*r1*r2)*ω*cos((-z+x2+x1)*ω))*gamma+((9*r1*r2-9*r1*r2*vxr)*vyr-9*r1*r2*vxr^2+18*r1*r2*vxr-6*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+((9*r1*r2*vxr-9*r1*r2)*vyr^3+(9*r1*r2*vxr^2-36*r1*r2*vxr+(24*r1+18)*r2)*vyr^2+(9*r1*r2*vxr^3-45*r1*r2*vxr^2+(36-84*r1)*r2*vxr+(126*r1-54)*r2-18*r1)*vyr+9*r1*r2*vxr^4-36*r1*r2*vxr^3+(18-96*r1)*r2*vxr^2+((270*r1-54)*r2-18*r1)*vxr+(24-100*r1)*r2+18*r1)*ω*cos((-z+x2+x1)*ω))/(100*sqrt(3)*r1*r2) + a4/tau
   du9 = (((9*r1*r2*vyr^2+(9*r1*r2*vxr-18*r1*r2)*vyr-9*r1*r2*vxr+6*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(135*r1*r2*vyr^2+(135*r1*r2*vxr-270*r1*r2)*vyr-135*r1*r2*vxr+90*r1*r2)*ω*cos((-z+x2+x1)*ω))*gamma+(-9*r1*r2*vyr^2+(18*r1*r2-9*r1*r2*vxr)*vyr+9*r1*r2*vxr-6*r1*r2)*ω*cos((-z+x2+x1)*ω)*sin((-z+x2+x1)*ω)+(9*r1*r2*vyr^4+(9*r1*r2*vxr-36*r1*r2)*vyr^3+(9*r1*r2*vxr^2-45*r1*r2*vxr+(18-96*r1)*r2)*vyr^2+(9*r1*r2*vxr^3-36*r1*r2*vxr^2+(36-84*r1)*r2*vxr+(270*r1-54)*r2-18*r1)*vyr-9*r1*r2*vxr^3+(24*r1+18)*r2*vxr^2+((126*r1-54)*r2-18*r1)*vxr+(24-100*r1)*r2+18*r1)*ω*cos((-z+x2+x1)*ω))/(100*sqrt(3)*r1*r2) + a5/tau
  

  #return SVector(du1, du2, du3, du4, du5, du6, du7, du8, du9)
  return SVector(f1, f2, f3, f4, f5, f6, f7, f8, f9)
  end 
  
    
  # @inline function source_terms_convergence_test(u, x, t, equations::PerturbationMomentSystem2D)
  #   w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  #   @unpack tau = equations
  
  #   du1 = -w0 * w0xx + w0x * w0x/ 3.0 - w0y * w0y/ 6.0 
  #   du2 = -w0 * w0yy - w0x * w0x/ 6.0 + w0y * w0y/ 3.0
  #   du3 = -w0 * w0xy + w0x * w0y/ 4.0 + w0y * w0x/ 4.0
  #   du4 = 2.0 * w1 * w0x / 3.0 + 4.0 * w0x * w0xx/15.0 + 4.0 * w0y * w0xy/ 15.0 - 2.0 * w0 * w1x / 3.0
  #   du5 = 2.0 * w1 * w0y / 3.0 + 4.0 * w0y * w0yy/ 15.0 + 4.0 * w0x * w0xy/15.0 - 2.0 * w0 * w1y / 3.0

  #   return SVector(0.0, 0.0, 0.0, 0.0, du1/tau, du2/tau, du3/tau, du4/tau, du5/tau)
  # end
  
  
  
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