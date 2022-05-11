function T = raytracing(f, dfdx, dfdy, dfdz, T0, v, step = 0.1, maxIter = 100, maxRef = 10)
% raytracing(f, T0, v) projects a ray from the origin point T0 in the 
% direction v and finds the points where the ray hits the plane, given by
% the function f, and returns them
% step is the amount the ray moves in every iteration
% maxIter is the maximum number of iterations
% maxRef is the maximum number of allowed reflections

% express the parameter z from the plane function
fz = @(x, y) -f(x, y, 0)/f(0, 0, 1);

g = @(t) f(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t);
gdot = @(t) v(1)*dfdx(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t) + v(2)*dfdy(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t) + v(3)*dfdz(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t);


% plot the plane
x = [10 -10 -10 10];
y = [10 10 -10 -10];
[x y] = meshgrid(x, y);
z = fz(x, y);
surf(x, y, z);
hold on

% initialize the vector of hit points as empty
T = [];

t = step;

% multiply the direction vector with the step size
stepv = step*v;

% initialize the previous and current observed point
prevT = T0;
currT = T0 + stepv;

% initialize the function sign in the previous and current point
prevSign = sign(f(T0(1), T0(2), T0(3)));
currSign = sign(f(currT(1), currT(2), currT(3)));

% set the number of iterations to 0 
iter = 0;

% while the sign remains unchanged and the number of iterations is smaller
% then the maximum number of iterations...
while (prevSign == currSign && iter <= maxIter)
   % plot the current point
   plot3(currT(1), currT(2), currT(3), '.r', 'markersize', 10);
  
   % increment the number of iterations
   iter++;
  
   t = t + step;
  
   % set the previous point to the current point
   prevT = currT;
   
   % generate a new point
   currT = currT + stepv;
   
   % set the previous sign to the current sign
   prevSign = currSign;
   
   % generate a new sign
   currSign = sign(f(currT(1), currT(2), currT(3)));
  
endwhile;

% if a hit point is found, use Newton's iteration to get a better approximation
if (iter <= maxIter)

  % define a starting approximation
  startingApproximation = (t + (t - step)) / 2;
  
  % use Newton's method to get a better approximation of the parameter
  u = newton(g, gdot, startingApproximation);
  
  % determine the better hit point
  U = T0 + v*u
  
  % add the hit point to the list of hit points
  T = [T; U];
  
  % plot the hit point
  plot3(U(1), U(2), U(3), '.m', 'markersize', 30);
endif          
endfunction

function [X, n] = newton(F, JF, X0, tol = 1e-10, maxit = 100)
%X = newton(F, JF, X0, tol, maxit) solves the (nonlinear) 
%system F(X) = 0 using the Newton's iteration with initial
%guess X0. (JF is the Jacobi matrix of F.)

for n = 1:maxit
	%Execute one step of Newton's iteration...
	X = X0 - feval(JF, X0)\feval(F, X0);
	%... and check if the new approximation is within prescribed tolerance.
	if(norm(X - X0) < tol)
		break;
	end
	X0 = X;
end

%A warning in case the last approximation is not within specified tolerance.
if(n == maxit)
	warning("no convergence after maxit iterations")
end
endfunction