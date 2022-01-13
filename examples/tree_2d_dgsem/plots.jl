x = LinRange(0,1,257)
y = LinRange(0,1,257)
slope = 15
B = zeros(length(y))
rho = zeros(length(y))
vx = zeros(length(y))
vy = zeros(length(x))
theta = zeros(length(y))
tmp = zeros(length(y))
diffi = zeros(length(y))
global bla = 0.0

for j = 1:length(x)

    a = 0.05
    z1 = 0.5
    z2 = 1.5
    s=0.2
    vy[j] = 0.01*sin(2 * pi * x[j])*(exp(-(x[j]-z1)^2/s^2)+ exp(-(x[j]-z2)^2/s^2))
    
end

for i = 1:length(y)
    # B[i] = tanh(slope * y[i] + 7.5) - tanh(slope * y[i] - 7.5)
    # rho[i] = 0.5 + 0.75 * B[i]
    
    # vx[i] = 0.5 * (B[i] - 1)
    # vy[i] = 0.1 * sin(2 * pi * x[i])
    # tmp = 1/(5/3-1.0) + 0.5*(rho[i]*vx[i]^2+rho[i]*vy[i]^2)
    # #theta[i] = (5/3-1.0)*(tmp/rho[i]- 0.5*(rho[i]*vx[i]^2 + rho[i]*vy[i]^2))
    # theta[i] = 5*1.0/(3*rho[i])

        a = 0.05
        z1 = 0.5
        z2 = 1.5
        s=0.2
        tmp[i] = tanh((y[i]-z1)/a)- tanh((y[i]-z2)/a)
        rho[i] = 1 + 0.5*tmp[i]
        vx[i] = tmp[i] -1
        local p = 10
        theta[i] = p/rho[i]
        diffi[i] = 0.01*((exp(-(-1.5 + y[i])^2)/25)*(-1.5 + y[i]) - (exp(-(-0.5 + y[i])^2)/25)*(-0.5 + y[i]))*sin(2*pi*x[i])/50
        global bla = bla + diffi[i]
end
    

p1 = plot(y, rho, title="rho")
p2 = plot(y, vx, title = "vx")
p3 = plot(x, vy, title ="vy")
p4 = plot(y, theta, title ="theta")

plot(p1,p2,p3,p4)


plot(y, diffi, title = "diff")
bla/257