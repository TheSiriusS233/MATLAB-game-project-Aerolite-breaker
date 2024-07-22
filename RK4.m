%% function to implement Runge-Kutta (RK) 4th order method
% inputs: ti = current time (s)
%         xi = previous state vector estimate (4x1): position and velocity
%          h = sample time (s)
%  w1,w2,w3,w4,c2,c3,c4,a21,a31,a32,a41,a42,a43 = "Classic" RK4 parameters
% outputs: x_new = new state vector estimate (4x1): position and velocity
function x_new=RK4(ti,xi,h,w1,w2,w3,w4,c2,c3,c4,...
a21,a31,a32,a41,a42,a43)
k1=h*f(ti,xi); % output of f(t,x) is 4x1, so k1 is 4x1
k2=h*f(ti+c2*h,xi+a21*k1);
k3=h*f(ti+c3*h,xi+a31*k1+a32*k2);
k4=h*f(ti+c4*h,xi+a41*k1+a42*k2+a43*k3);
x_new=xi+w1*k1+w2*k2+w3*k3+w4*k4;
end

