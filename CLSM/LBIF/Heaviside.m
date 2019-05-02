function H = Heaviside(u,epsilon)
H = 0.5*(1+2/pi*atan(u./epsilon));