using LaTeXStrings
x = LinRange(-1,1,257)
y = LinRange(-1,1,257)
slope = 15
B = zeros(length(y))
rho = zeros(length(y))
vx = zeros(length(y))
vy = zeros(length(x))
theta = zeros(length(y))
tmp = zeros(length(y))
diffi = zeros(length(y))
dtheta = zeros(length(y)) 
global bla = 0.0



for i = 1:length(y)
    B[i] = tanh(slope * y[i] + 7.5) - tanh(slope * y[i] - 7.5)
    rho[i] = (0.5 + 0.75 * B[i])
    
    vx[i] = 0.5 * (B[i] - 1)
    vy[i] = 0.1 * sin(2 * pi * x[i])
   
    theta[i] =  1/rho[i]

    dtheta[i] = -(0.75*(15*sech(15*y[i]+7.5)^2-15*sech(15*y[i]-7.5)^2))/(0.75*(tanh(15*y[i]+7.5)-tanh(15*y[i]-7.5))+0.5)^2

        # a = 0.05
        # z1 = 0.5
        # z2 = 1.5
        # s=0.2
        # tmp[i] = tanh((y[i]-z1)/a)- tanh((y[i]-z2)/a)
        # rho[i] = 1 + 0.5*tmp[i]
        # vx[i] = tmp[i] -1
        # vy[i] = 0.01*sin(2 * pi * y[i])*(exp(-(y[i]-z1)^2/s^2)+ exp(-(y[i]-z2)^2/s^2))
        # local p = 10
        # theta[i] = p/rho[i]
        # diffi[i] = 0.01*((exp(-(-1.5 + y[i])^2)/25)*(-1.5 + y[i]) - (exp(-(-0.5 + y[i])^2)/25)*(-0.5 + y[i]))*sin(2*pi*x[i])/50
    global bla = bla + theta[i]
end
    

p1 = plot(y, rho, title="ρ", legend = false, xlabel = "y", titlefontsize = 20, xlabelfontsize=12, tickfontsize=12)
p2 = plot(y, vx, title = "v\u2093", legend = false, xlabel = "y", titlefontsize = 20, xlabelfontsize=12, tickfontsize=12)
p3 = plot(x, vy, title = "vy", legend = false, xlabel = "x", titlefontsize = 20, xlabelfontsize=12, tickfontsize=12)

plot(p1,p2,p3, layout = (1,3), guidefontsize = 6, legendfontsize= 5, titlefontsize = 15, tickfontsize=9, linewidth = 2, guidefont=font(15), size=(400,250))

