
function xdot=duffing(t,x)

global delta alpha gamma beta omega

xdot(1)=-delta*x(1)+alpha*x(2)-beta*x(2)^3+gamma*cos(omega*t);
xdot(2)=x(1);

xdot=xdot';

end
