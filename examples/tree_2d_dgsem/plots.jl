x = LinRange(-1,1,257)
y = LinRange(-1,1,257)
slope = 15
B = zeros(length(y))
rho = zeros(length(y))
vx = zeros(length(y))
vy = zeros(length(y))
theta = zeros(length(y))
p = 1.0
global bla = 0.0

for i = 1:length(y)
    B[i] = tanh(slope * y[i] + 7.5) - tanh(slope * y[i] - 7.5)
    rho[i] = 0.5 + 0.75 * B[i]
    
    vx[i] = 0.5 * (B[i] - 1)
    vy[i] = 0.1 * sin(2 * pi * x[i])
    theta[i] = p/(rho[i]*8.314)
    global bla = bla + theta[i]
end

p1 = plot(x, rho, title="rho")
p2 = plot(x, vx, title = "vx")
p3 = plot(x, vy, title ="vy")
p4 = plot(x, theta, title ="theta")

plot(p1,p2,p3,p4)