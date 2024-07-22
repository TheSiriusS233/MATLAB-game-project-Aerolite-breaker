%% function containing the differential equation dx/dt=f(t,x) to solve
% inputs:     t = time (s)
%             x = state vector (4x1) consists of mass position and velocity
function dxdt=f(t,x)
global mu;
r=norm(x(1:2),2); % magnitude of position vector
dxdt=[zeros(2,2) eye(2);
-(mu/r^3)*eye(2) zeros(2,2)]*x; % dxdt is a 4x4 vector
end