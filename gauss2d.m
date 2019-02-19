function z = gauss2d(x,y,mx,my,sx,sy)

X = ((x-mx).^2)/(2*sx^2);
Y = ((y-my).^2)/(2*sy^2);
z = exp(-(X+Y));