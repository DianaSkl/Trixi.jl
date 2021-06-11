drho = 1
rho_r = 2
dv_x = 4
dv_y = 4
theta_r = 4
sigma_xx = 0
sigma_xy = 0
sigma_yy = 0
dtheta = -14.4
rho = drho + rho_r
dv2 = 9
q_x = 138
q_y = 138


w0 = 1.0 + drho / rho_r
w0x = dv_x / sqrt(theta_r) + (drho * dv_x)/(rho_r * sqrt(theta_r))
w0y = dv_y / sqrt(theta_r) + (drho * dv_y)/(rho_r * sqrt(theta_r))
w1 = - (rho * (dv_x *dv_x + dv_y * dv_y) )/(3.0 * (rho_r * theta_r)) - (drho * dtheta)/(rho_r * theta_r) - dtheta / theta_r

w0xx = 0.5 * sigma_xx/(rho_r * theta_r) + (2.0 * dv_x * dv_x - dv_y * dv_y)/(6.0 * theta_r) + (2.0 * drho * dv_x * dv_x - drho * dv_y * dv_y)/(6.0 * rho_r * theta_r)  
w0yy = 0.5 * sigma_yy/(rho_r * theta_r) + (2.0 * dv_y * dv_y - dv_x * dv_x)/(6.0 * theta_r) + (2.0 * drho * dv_y * dv_y - drho * dv_x * dv_x)/(6.0 * rho_r * theta_r)
w0xy = 0.5 * sigma_xy/(rho_r * theta_r) + 0.5 * ((rho_r + drho) * (dv_x * dv_y))/(rho_r * theta_r) 
w1x = - 2.0 * q_x / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xx * dv_x + sigma_xy * dv_y))/(5.0*rho_r* sqrt(theta_r).^3.0) - (dtheta * dv_x)/ sqrt(theta_r).^3.0 - (dv_x * dv2)/(5.0*sqrt(theta_r).^3.0) - (drho*dtheta*dv_x)/(rho_r*sqrt(theta_r).^3.0)- (drho*dv_x*dv2)/(5.0*rho_r*sqrt(theta_r).^3.0) 
w1y = - 2.0 * q_y / (5.0* rho_r * sqrt(theta_r).^3.0) - (2.0 * (sigma_xy * dv_x + sigma_yy * dv_y))/(5.0*rho_r* sqrt(theta_r).^3.0) - (dtheta * dv_y)/ sqrt(theta_r).^3.0 - (dv_y * dv2)/(5.0*sqrt(theta_r).^3.0) -(drho*dtheta*dv_y)/(rho_r*sqrt(theta_r).^3.0)- (drho*dv_y*dv2)/(5.0*rho_r*sqrt(theta_r).^3.0)


du1 = -w0 * w0xx + w0x * w0x/ 3.0 - w0y * w0y/ 6.0 
du2 = -w0 * w0xy + w0x * w0y/ 4.0 + w0y * w0x/ 4.0
du3 = -w0 * w0yy - w0x * w0x/ 6.0 + w0y * w0y/ 3.0
du4 = 2.0 * w1 * w0x / 3.0 + 4.0 * w0x * w0xx/15.0 + 4.0 * w0y * w0xy/ 15.0 - 2.0 * w0 * w1x / 3.0
du5 = 2.0 * w1 * w0y / 3.0  + 4.0 * w0y * w0yy/ 15.0 + 4.0 * w0x * w0xy/15.0 - 2.0 * w0 * w1y / 3.0
  