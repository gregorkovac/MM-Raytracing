function T = raytracing(f, T0, v, step, maxRef)
% raytracing(f, T0, v) projects a ray from the origin point T0 in the 
% direction v and finds the points where the ray hits the plane, given by
% the function f.
% step is the ???
% maxRef is the maximum number of reflections

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