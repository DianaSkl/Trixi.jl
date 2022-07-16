# This file was used in the following master thesis:
# - Diana Sklema (2022)
#   "Untersuchung eines gestoerten Momentensystems in der kompressiblen Stroeomungsmechanik"
#   University of Cologne, advisors: Gregor Gassner, Michael Schlottke-Lakemper


struct MomentSystem2D{RealT<:Real} <: AbstractMomentSystem{2, 9} 
    vxr::RealT
    vyr::RealT
    theta_r::RealT
    rho_r::RealT
    tau::RealT
end

varnames(::typeof(cons2prim), ::MomentSystem2D) = ("ρ", "vx", "σxx", "p", "vy", "σyy", "qx", "qy","σxy")
varnames(::typeof(cons2cons), ::MomentSystem2D) = ("w0", "w0x", "w0y", "w1", "w0xx", "w0yy", "w0xy", "w1x", "w1y")
  

# Calculate 2D flux 
@inline function flux(u, orientation::Integer, equations::MomentSystem2D)
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
  
"""
flux_ds(u_ll, u_rr, orientation, equations::MomentSystem2D)
Kinetic energy preserving two-point flux for the 2d perturbated moment system
"""  
@inline function flux_ds(u_ll, u_rr, orientation::Integer, equations::MomentSystem2D)

  @unpack vxr, vyr, theta_r, rho_r = equations

  # Unpack left and right state
  a_ll = cons2prim(u_ll, equations)
  a_rr = cons2prim(u_rr, equations)
  rho_ll = a_ll[1]
  vx_ll = a_ll[2]
  vy_ll = a_ll[5]
  p_ll = a_ll[4]
  rho_rr = a_rr[1]
  vx_rr = a_rr[2]
  vy_rr = a_rr[5]
  p_rr = a_rr[4]
  

  # Average each factor of products in flux
  rho_avg = 1/2 * (rho_ll + rho_rr)
  vx_avg  = 1/2 * (vx_ll +  vx_rr)
  vy_avg  = 1/2 * (vy_ll +  vy_rr)
  p_avg = 1/2 * (p_ll + p_rr)

  W = SVector(prim2cons((rho_avg, vx_avg, vy_avg, p_avg ),equations))      
  f1, f2, f3, f4, f5, f6, f7, f8, f9 = flux(W, orientation, equations)
       
  return SVector(f1, f2, f3, f4, f5, f6, f7, f8, f9)
end

  
@inline function cons2prim(u, equations::MomentSystem2D)
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
    @unpack vxr, vyr, theta_r, rho_r = equations

    rho = w0 * rho_r
    vx = vxr + w0x * sqrt(theta_r) / w0
    dvx = vx -vxr
    vy = vyr + w0y * sqrt(theta_r) / w0
    dvy = vy - vyr 
    theta = theta_r - (w0x^2 * theta_r)/(3 * w0^2) - (w0y^2 * theta_r)/(3 * w0^2) - (w1 * theta_r)/ w0
    dtheta = theta -theta_r
    p = rho*theta
    s3theta = sqrt(theta_r^3)

    sigmaxx = 2*rho_r*theta_r*w0xx - rho*(2*dvx^2 - dvy^2)/3
    sigmayy = 2*rho_r*theta_r*w0yy - rho*(2*dvy^2 - dvx^2)/3
    sigmaxy = 2*rho_r*theta_r*w0xy - rho*dvx*dvy
    qx = -w1x*rho_r*s3theta*5/2 - sigmaxx*dvx - sigmaxy*dvy - dtheta*dvx*rho*5/2 - 2*rho*(dvx^3 + dvx*dvy^2)
    qy = -w1y*rho_r*s3theta*5/2 - sigmayy*dvy - sigmaxy*dvx - dtheta*dvy*rho*5/2 - 2*rho*(dvy^3 + dvy*dvx^2)
    
    return SVector(rho, vx, sigmaxx, p, vy, sigmayy, qx, qy, sigmaxy)
end
  
# Convert primitive to conservative variables
@inline function prim2cons(prim, equations::MomentSystem2D)
  rho, vx, vy, p = prim
  @unpack vxr, vyr, theta_r, rho_r = equations

  sxx = sxy = syy = qx = qy = 0.0

  dvx = vx - vxr 
  dvy = vy - vyr 
  theta = p/rho
  dtheta = theta - theta_r 

  w0 = rho/rho_r
  w0x = (rho * dvx)/(rho_r * sqrt(theta_r))
  w0y = (rho * dvy)/(rho_r * sqrt(theta_r))
  w1 = -dtheta*rho/(rho_r*theta_r) - (rho*(dvx^2 + dvy^2))/(3*rho_r*theta_r)
  w0xx = 0.5*sxx/(rho_r*theta_r) + 0.5*(rho*(2*(dvx^2)/3- (dvy^2)/3))/(rho_r*theta_r)
  w0yy = 0.5*syy/(rho_r*theta_r) + 0.5*(rho*(2*(dvy^2)/3- (dvx^2)/3))/(rho_r*theta_r)
  w0xy = 0.5*sxy/(rho_r*theta_r) + 0.5*(rho*dvx*dvy)/(rho_r*theta_r)
  w1x = -2*(qx + sxx*dvx +sxy*dvy)/(5*rho_r*sqrt(theta_r^3)) - (dtheta*dvx*rho)/(rho_r*sqrt(theta_r^3))- (rho*(dvx^3 + dvx*dvy^2))/(5*sqrt(theta_r^3)*rho_r)
  w1y = -2*(qy + sxy*dvx *syy*dvy)/(5*rho_r*sqrt(theta_r^3)) - (dtheta*dvy*rho)/(rho_r*sqrt(theta_r^3))- (rho*(dvx^2*dvy + dvy^3))/(5*sqrt(theta_r^3)*rho_r)

  return SVector(w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y)
end
  
  # Convert conservative variables to primitive
  @inline function cons2entropy(u, equations::MomentSystem2D)
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  
    return SVector( w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y)
  end
  
  # Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
  # maximum velocity magnitude plus the maximum speed of sound
  @inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::MomentSystem2D)
    @unpack vxr, theta_r = equations
  
    ab1 = abs(vxr) + 2.0 * sqrt(5 * theta_r / 3)
  
    return ab1
  
  end
  
  @inline function max_abs_speeds(u, equations::MomentSystem2D)
    @unpack vxr, vyr, theta_r = equations
  
    ab1 = abs(vxr) + 2.0 * sqrt(5 * theta_r / 3)
    ab2 = abs(vyr) + 2.0 * sqrt(5 * theta_r / 3)
  
    return SVector(ab1, ab2)
    
  end
  
  # possible indicator ρp for the SC-Methode
  @inline function density_pressure(u, equations::MomentSystem2D)
    a = cons2prim(u, equations)
    rho = a[1]
    p = a[4]
    rho_times_p = rho*p
    return rho_times_p
  end
  
  # possible indicator ρ for the SC-Methode
  @inline function density(u, equations::MomentSystem2D)
    a = cons2prim(u, equations)
    rho = a[1]
    return rho
  end
     
  """
  initial_condition_convergence_test(x, t, equations::MomentSystem2D)
  periodical initial condition for the method of manufactured solution
  """
  function initial_condition_convergence_test(x, t, equations::MomentSystem2D)
    c = 2
    A = 0.1
    L = 2
    f = 1/L
    ω = 2 * pi * f
    ini = c + A * sin(ω * (x[1] + x[2] - t))
  
   
    @unpack vxr, vyr, theta_r, rho_r, tau = equations
  
    gamma = 5/3

    rho = ini
    drho = rho-rho_r
    vx = 1.0
    vy = 1.0
    dvx = vx-vxr
    dvy = vy-vyr
    theta = (gamma-1)*(ini-(vx^2 + vy^2)/2)
    dtheta = theta -  theta_r

    sxx = 0.1
    syy = 0.1
    sxy = 0.1

    qx = 0.1
    qy = 0.1

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
  
  
  
  @inline function source_terms_convergence_test(u, x, t, equations::MomentSystem2D)
    
  @unpack vxr, vyr, theta_r, rho_r, tau = equations
  
   #convtest zeit&raum
   #f1,f2,f3,f4,f5,f6,f7,f8,f9 = source_convergence_spatialtime(u, x, t, equations::MomentSystem2D)
    
   #produktionen
   f1,f2,f3,f4,f5,f6,f7,f8,f9 = source_productions(u, equations::MomentSystem2D)

   
   return (f1,f2,f3,f4,f5,f6,f7,f8,f9)

 end 
  
    
  
  @inline function source_productions(u, equations::MomentSystem2D)

    @unpack tau = equations 
    w0, w0x, w0y, w1, w0xx, w0yy, w0xy, w1x, w1y = u
  

     p1 = -w0 * w0xx + w0x * w0x/ 3.0 - w0y * w0y/ 6.0 
     p2 = -w0 * w0yy - w0x * w0x/ 6.0 + w0y * w0y/ 3.0
     p3 = -w0 * w0xy + w0x * w0y/ 4.0 + w0y * w0x/ 4.0
     p4 = 2.0 * w1 * w0x / 3.0 + 4.0 * w0x * w0xx/15.0 + 4.0 * w0y * w0xy/ 15.0 - 2.0 * w0 * w1x / 3.0
     p5 = 2.0 * w1 * w0y / 3.0 + 4.0 * w0y * w0yy/ 15.0 + 4.0 * w0x * w0xy/15.0 - 2.0 * w0 * w1y / 3.0

    return (0,0,0,0,p1/tau,p2/tau,p3/tau,p4/tau,p5/tau)
  end 



  @inline function source_convergence_spatialtime(u,x,t,equations::MomentSystem2D)
    @unpack vxr, vyr, theta_r, rho_r, tau = equations
    c = 2
    A = 0.1
    L = 2
    f = 1/L
    m = 1/L
    ω = 2 * pi * f
    x1 = x[1]
    x2 = x[2]
  
    vx = 1.0
    vy = 1.0

    p1 = 2/3
    p2 = 1/3
  

    dw0_x = (A*ω*cos((x2+x1-t)*ω))/c

    dw0x_x = (A*(vx-vxr)*ω*cos((x2+x1-t)*ω))/(c*sqrt(theta_r))

    dw0y_x = (A*(vy-vyr)*ω*cos((x2+x1-t)*ω))/(c*sqrt(theta_r))

    dw1_x = (A*ω*cos((x2+x1-t)*ω)*(1-p1*(A*sin((x2+x1-t)*ω)+c)))/(c*theta_r)-(A*p1*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c))/(c*theta_r)-(A*((vy-vyr)^2+(vx-vxr)^2)*ω*cos((x2+x1-t)*ω))/(3*c*theta_r)

    dw0xx_x = (A*(p1*(vx-vxr)^2-p2*(vy-vyr)^2)*ω*cos((x2+x1-t)*ω))/(c^2*theta_r)

    dw0yy_x = (A*(p1*(vy-vyr)^2-p2*(vx-vxr)^2)*ω*cos((x2+x1-t)*ω))/(c^2*theta_r)

    dw0xy_x = (A*(vx-vxr)*(vy-vyr)*ω*cos((x2+x1-t)*ω))/(2*c*theta_r)

    dw1x_x = -(A*(vx-vxr)*ω*cos((x2+x1-t)*ω)*(p1*(A*sin((x2+x1-t)*ω)+c)-1))/(c*theta_r^(3/2))-(A*p1*(vx-vxr)*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c))/(c*theta_r^(3/2))-(A*((vx-vxr)*(vy-vyr)^2+(vx-vxr)^3)*ω*cos((x2+x1-t)*ω))/(5*c*theta_r^(3/2))
    
    dw1y_x = -(A*(vy-vyr)*ω*cos((x2+x1-t)*ω)*(p1*(A*sin((x2+x1-t)*ω)+c)-1))/(c*theta_r^(3/2))-(A*p1*(vy-vyr)*ω*cos((x2+x1-t)*ω)*(A*sin((x2+x1-t)*ω)+c))/(c*theta_r^(3/2))-(A*((vx-vxr)^2*(vy-vyr)+(vy-vyr)^3)*ω*cos((x2+x1-t)*ω))/(5*c*theta_r^(3/2))
    
    
   f1 = -dw0_x + vxr * dw0_x  + sqrt(theta_r) * dw0x_x + vyr * dw0_x + sqrt(theta_r) * dw0y_x
   f2 = -dw0x_x + vxr * dw0x_x + sqrt(theta_r) * (dw0_x - dw1_x) + 2.0 * sqrt(theta_r) * dw0xx_x + vyr * dw0x_x + 2.0 * sqrt(theta_r) * dw0xy_x
   f3 = -dw0y_x + vxr * dw0y_x + 2.0 * sqrt(theta_r) * dw0xy_x + vyr * dw0y_x + sqrt(theta_r) * (dw0_x - dw1_x + 2.0 * dw0yy_x)
   f4 = -dw1_x + vxr * dw1_x + sqrt(theta_r) * (5.0 * dw1x_x - 2.0 * dw0x_x)/3.0 +vyr * dw1_x + sqrt(theta_r) * (5.0 * dw1y_x - 2.0 * dw0y_x)/3.0
   f5 = -dw0xx_x + vxr * dw0xx_x + 2.0 * sqrt(theta_r)* (dw0x_x - dw1x_x)/3.0  + vyr * dw0xx_x + sqrt(theta_r) * (dw1y_x - dw0y_x)/3.0 
   f6 = -dw0yy_x + vxr * dw0yy_x + sqrt(theta_r) * (dw1x_x - dw0x_x)/3.0 + vyr * dw0yy_x + sqrt(theta_r) * 2.0 * (dw0y_x - dw1y_x)/3.0 
   f7 = -dw0xy_x + vxr * dw0xy_x + sqrt(theta_r) * (dw0y_x - dw1y_x)/2.0 + vyr * dw0xy_x + sqrt(theta_r)  * (dw0x_x - dw1x_x)/2.0 
   f8 = -dw1x_x + vxr * dw1x_x + sqrt(theta_r) * (dw1_x - 4.0 * dw0xx_x/5.0) + vyr * dw1x_x - sqrt(theta_r) * 4.0 * dw0xy_x / 5.0 
   f9 = -dw1y_x + vxr * dw1y_x - 4.0 * sqrt(theta_r) * dw0xy_x / 5.0 + vyr * dw1y_x + sqrt(theta_r) * (dw1_x - 4.0 * dw0yy_x / 5.0) 

    return (f1,f2,f3,f4,f5,f6,f7,f8,f9)
 end  
  
  function shocktube(x, equations::MomentSystem2D)
    rho_r = 1
  
    if (x[1] < -20 || (x[1]<0 && x[1]>=-10) || (x[1]>=10 && x[1]<20))
      drho = 3 - rho_r
    elseif ((x[1] < -10 && x[1] >= -20) || (x[1]>=0 && x[1] <10) || x[1] >= 20)
      drho = 1 - rho_r
    end
  
    return drho
    end
  
  function initial_condition_constant(x, t, equations::MomentSystem2D)
  
    return SVector(init_w0(x, t, equations::MomentSystem2D), 
                   init_w0x(x, t, equations::MomentSystem2D), 
                   init_w0y(x, t, equations::MomentSystem2D), 
                   init_w1(x, t, equations::MomentSystem2D),
                   init_w0xx(x, t, equations::MomentSystem2D),
                   init_w0yy(x, t, equations::MomentSystem2D),
                   init_w0xy(x, t, equations::MomentSystem2D), 
                   init_w1x(x, t, equations::MomentSystem2D),
                   init_w1y(x, t, equations::MomentSystem2D))
  end
  
  
  @inline function init_w0(x, t, equations::MomentSystem2D)
   
    rho_r = 1
    drho = shocktube(x, equations)
  
   
    w0 = 1 + drho / rho_r
    
   return w0
  end
  
  
  @inline function init_w0x(x, t, equations::MomentSystem2D)
    
    @unpack theta_r = equations 
  
    rho_r = 1
    drho = shocktube(x, equations)
  
    dv_x = 0
  
    w0x = dv_x / sqrt(theta_r) + (drho * dv_x)/(rho_r * sqrt(theta_r))
  
   return w0x
  end
  
  @inline function init_w0y(x, t, equations::MomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    
    w0y = dv_y / sqrt(theta_r) + (drho * dv_y)/(rho_r * sqrt(theta_r))
    
   return w0y
  end
  
  @inline function init_w1(x, t, equations::MomentSystem2D)
    
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
  
  
  @inline function init_w0xx(x, t, equations::MomentSystem2D)
    
    @unpack theta_r = equations 
  
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0
  
    sigma_xx = 0
   
    w0xx = 0.5 * sigma_xx/(rho_r * theta_r) + (2.0 * dv_x * dv_x - dv_y * dv_y)/(6.0 * theta_r) + (2.0 * drho * dv_x * dv_x - drho * dv_y * dv_y)/(6.0 * rho_r * theta_r)  
    
   return w0xx
  end
  
  @inline function init_w0yy(x, t, equations::MomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0
    sigma_yy = 0
    
    w0yy = 0.5 * sigma_yy/(rho_r * theta_r) + (2.0 * dv_y * dv_y - dv_x * dv_x)/(6.0 * theta_r) + (2.0 * drho * dv_y * dv_y - drho * dv_x * dv_x)/(6.0 * rho_r * theta_r)
    
   return w0yy
  end
  
  
  @inline function init_w0xy(x, t, equations::MomentSystem2D)
    
    @unpack theta_r = equations 
    rho_r = 1
    drho = shocktube(x, equations)
    dv_y = 0
    dv_x = 0
    sigma_xy = 0
  
    
    w0xy = 0.5 * sigma_xy/(rho_r * theta_r) + 0.5 * ((rho_r + drho) * (dv_x * dv_y))/(rho_r * theta_r)
    
   return w0xy
  end
  
  
  
  @inline function init_w1x(x, t, equations::MomentSystem2D)
    
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
  
  
  
  @inline function init_w1y(x, t, equations::MomentSystem2D)
    
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