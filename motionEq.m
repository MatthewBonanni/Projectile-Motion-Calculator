function vdot = motionEq(t, v, constants)

b = constants.b;
c = constants.c;
m = constants.mass;
g = constants.g;

vdot = zeros(2,1);
vdot(1) = (-b*v(1) - c * sqrt(v(1)^2 + v(2)^2) * v(1)) / m;
vdot(2) = (-m*g - b*v(1) - c * sqrt(v(1)^2 + v(2)^2) * v(2)) / m;

end