function color = raytracing(f, f1, f2, dfdx, dfdy, dfdz, df1dx, df1dy, df1dz, df2dx, df2dy, df2dz, T0, v, lightOrigin, step, maxIter, testRef, def_colors, colorOption1, colorOption2, colorOption3, skyColor)
% raytracing(f, T0, v) projects a ray from the origin point T0 in the 
% direction v and finds the point where the ray hits the planes, given by
% the functions f1, f2 and f3, and returns a 3x1 vector color, representing the color of the
% point intersection of the ray and an object.
% dfdx, dfdy, dfdz, df1dx, df1dy, df1dz, df2dx, df2dy, df2dz are partial
% derivatives of f1, f2 and f3 by x, y and z.
% lightOrigin is the position of the light source.
% step is the amount the ray moves in every iteration.
% maxIter is the maximum number of iterations.
% testRef is used for checking if we have already sent a ray from the hit
% point to the light source (when calling the function it should be 0).
% def_colors is a 3x3 matrix, where rows represent the default color of a
% corresponding function.
% colorOption1, colorOption2 and colorOption3 are used to determine how we
% color the hit point. 0 means default color, 1 means random color, 2 means
% checkerboard, 3 means gradient, 4 means reflective surface 
% and everything else means black.

% initialize color to black
color = [0; 0; 0];

% initialize default color
def_color = def_colors(:, 1);

% compute functions and derivatives needed for Newton's method
g = @(t) f(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t);
gdot = @(t) v(1)*dfdx(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t) + v(2)*dfdy(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t) + v(3)*dfdz(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t);

g1 = @(t) f1(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t);
g1dot = @(t) v(1)*df1dx(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t) + v(2)*df1dy(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t) + v(3)*df1dz(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t);

g2 = @(t) f2(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t);
g2dot = @(t) v(1)*df2dx(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t) + v(2)*df2dy(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t) + v(3)*df2dz(T0(1) + v(1)*t, T0(2) + v(2)*t, T0(3) + v(3)*t);

% initialize the vector of hit points as empty
T = [];

% initialize the parameter t
t = step;

% multiply the direction vector with the step size
stepv = step*v;

% initialize the previous and current observed point
prevT = T0;
currT = T0 + stepv;

% initialize the function sign in the previous and current point
prevSign = sign(f(T0(1), T0(2), T0(3)));
currSign = sign(f(currT(1), currT(2), currT(3)));

prevSign1 = sign(f1(T0(1), T0(2), T0(3)));
currSign1 = sign(f1(currT(1), currT(2), currT(3)));

prevSign2 = sign(f2(T0(1), T0(2), T0(3)));
currSign2 = sign(f2(currT(1), currT(2), currT(3)));

% set the number of iterations to 0 
iter = 0;

% initialize ret to an empty array. ret is used for finding objects between
% the hit point and the light source
ret = [];

% while the sign remains unchanged and the number of iterations is smaller
% then the maximum number of iterations...
while (iter <= maxIter)
   % increment the number of iterations
   iter = iter + 1;
  
   t = t + step;
  
   % set the previous point to the current point
   prevT = currT;
   
   % generate a new point
   currT = currT + stepv;
   
   % set the previous sign to the current sign
   prevSign = currSign;
   prevSign1 = currSign1;
   prevSign2 = currSign2;
   
   % generate a new sign
   currSign = sign(f(currT(1), currT(2), currT(3)));
   currSign1 = sign(f1(currT(1), currT(2), currT(3)));
   currSign2 = sign(f2(currT(1), currT(2), currT(3)));
   
   % check for intersection with f
   if (currSign ~= prevSign)
    % define a starting approximation
    startingApproximation = (t + (t - step)) / 2;
  
    % use Newton's method to get a better approximation of the parameter
    u = newton(g, gdot, startingApproximation, 1e-10, 100);
  
    % determine the better hit point
    U = T0 + v*u;
  
    % add the hit point to the list of hit points
    T = [T; U];

    % if we haven't already done a second reflection, send another ray from
    % the hit point to the light source
    if (testRef == 0)
        ret = raytracing(@(x,y,z) 0, f1, f2, @(x,y,z) 0, @(x,y,z) 0, @(x,y,z) 0, df1dx, df1dy, df1dz, df2dx, df2dy, df2dz, U, lightOrigin - U, lightOrigin, step, maxIter, 1, def_colors, colorOption1, colorOption2, colorOption3, skyColor);
    end
  
    % compute the angle between the normal and lightsource and the
    % reflection vector
    [cos_reflAngle, reflVec] = reflectionAngle(v, U, lightOrigin, dfdx, dfdy, dfdz);
    
    % choose color for coloring the point   
    switch colorOption1
        case 0
            def_color = def_colors(:, 1);
        case 1
            def_color = rand(3, 1);
        case 2
            if (mod(floor(U(1).*2),2) == 0 && mod(floor(U(2).*2),2) == 0)
                def_color =   [61; 53; 54]./255;
            elseif (mod(floor(U(1).*2),2) ~= 0 && mod(floor(U(2).*2),2) == 0)
                def_color =  [255; 255; 255]./255;
            elseif (mod(floor(U(1).*2),2) == 0 && mod(floor(U(2).*2),2) ~= 0)
                def_color =  [186; 177; 177]./255;
            else 
                def_color = [77; 73; 74]./255;
            end
        case 3
            def_color = [cos(U(1)); U(2); U(3)];
        case 4
            def_color = raytracing(@(x,y,z) 0, f1, f2, @(x,y,z) 0, @(x,y,z) 0, @(x,y,z) 0, df1dx, df1dy, df1dz, df2dx, df2dy, df2dz, U, reflVec, lightOrigin, step, maxIter, 1, def_color, colorOption1, colorOption2, colorOption3, skyColor);
        otherwise
            def_color = [0; 0; 0];
    end
  break;
  end
 
  % check for intersection with f1
  if (currSign1 ~= prevSign1)
    % define a starting approximation
    startingApproximation = (t + (t - step)) / 2;
  
    % use Newton's method to get a better approximation of the parameter
    u = newton(g1, g1dot, startingApproximation, 1e-10, 100);
  
    % determine the better hit point
    U = T0 + v*u;
  
    % add the hit point to the list of hit points
    T = [T; U];

    % if we haven't already done a second reflection, send another ray from
    % the hit point to the light source
    if (testRef == 0)
        ret = raytracing(f, @(x,y,z) 0, f2, dfdx, dfdy, dfdz, @(x,y,z) 0, @(x,y,z) 0, @(x,y,z) 0, df2dx, df2dy, df2dz, U, lightOrigin - U, lightOrigin, step, maxIter, 1, def_colors, colorOption1, colorOption2, colorOption3, skyColor);
    end
   
    % compute the angle between the normal and lightsource and the
    % reflection vector
    [cos_reflAngle, reflVec] = reflectionAngle(v, U, lightOrigin, df1dx, df1dy, df1dz);

    % choose color for coloring the point   
    switch colorOption2
        case 0
            def_color = def_colors(2, :);
        case 1
            def_color = rand(3, 1);
        case 2
            if (mod(floor(U(1).*2),2) == 0 && mod(floor(U(2).*2),2) == 0)
                def_color =   [61; 53; 54]./255;
            elseif (mod(floor(U(1).*2),2) ~= 0 && mod(floor(U(2).*2),2) == 0)
                def_color =  [255; 255; 255]./255;
            elseif (mod(floor(U(1).*2),2) == 0 && mod(floor(U(2).*2),2) ~= 0)
                def_color =  [186; 177; 177]./255;
            else 
                def_color = [77; 73; 74]./255;
            end
        case 3
            def_color = [cos(U(1)); U(2); U(3)];
        case 4
            def_color = raytracing(f, @(x,y,z) 0, f2, dfdx, dfdy, dfdz, @(x,y,z) 0, @(x,y,z) 0, @(x,y,z) 0, df2dx, df2dy, df2dz, U, reflVec, lightOrigin, step, maxIter, 1, def_color, colorOption1, colorOption2, colorOption3, skyColor);
        otherwise
            def_color = [0; 0; 0];
    end

    break;
  end

  % check for intersection with f1
  if (currSign2 ~= prevSign2)
    % define a starting approximation
    startingApproximation = (t + (t - step)) / 2;
  
    % use Newton's method to get a better approximation of the parameter
    u = newton(g2, g2dot, startingApproximation, 1e-10, 100);
  
    % determine the better hit point
    U = T0 + v*u;
  
    % add the hit point to the list of hit points
    T = [T; U];

    if (testRef == 0)
        ret = raytracing(f, f1, @(x, y, z) 0, dfdx, dfdy, dfdz, df1dx, df1dy, df1dz, @(x, y, z) 0, @(x, y, z) 0, @(x, y, z) 0, U, lightOrigin - U, lightOrigin, step, maxIter, 1, def_colors, colorOption1, colorOption2, colorOption3, skyColor);
    end
   
    % compute the angle between the normal and lightsource and the
    % reflection vector
    [cos_reflAngle, reflVec] = reflectionAngle(v, U, lightOrigin, df2dx, df2dy, df2dz);
    
    % choose color for coloring the point   
    switch colorOption3
        case 0
            def_color = def_colors(3, :);
        case 1
            def_color = rand(3, 1);
        case 2
            if (mod(floor(U(1).*2),2) == 0 && mod(floor(U(2).*2),2) == 0)
                def_color =   [61; 53; 54]./255;
            elseif (mod(floor(U(1).*2),2) ~= 0 && mod(floor(U(2).*2),2) == 0)
                def_color =  [255; 255; 255]./255;
            elseif (mod(floor(U(1).*2),2) == 0 && mod(floor(U(2).*2),2) ~= 0)
                def_color =  [186; 177; 177]./255;
            else 
                def_color = [77; 73; 74]./255;
            end
        case 3
            def_color = [cos(U(1)); U(2); U(3)];
        case 4
            def_color = raytracing(f, f1, @(x,y,z) 0, dfdx, dfdy, dfdz, df1dx, df1dy, df1dz, @(x,y,z) 0, @(x,y,z) 0, @(x,y,z) 0, U, reflVec, lightOrigin, step, maxIter, 1, def_color, colorOption1, colorOption2, colorOption3, skyColor);
        otherwise
            def_color = [0; 0; 0];
    end

    break;
    end
   
end

if (iter > maxIter)
    cos_reflAngle = 0;
    % optional: for reflective surfaces
    % if there are more iterations than allowed, set color to 
    % default color of background
    % color = [0.2; 0.8; 1];
    % return;
end
if (isempty(T) && testRef == 0)
    % in this case we haven't found and intersection
    color = skyColor;
elseif (isempty(ret) || (ret(1,1) == 0 && ret(2, 1) == 0 && ret(3, 1) == 0))
    % normal coloring
    color = (2.*cos_reflAngle).*def_color.*(iter./200); 
else
    % coloring if an object is between the intersection and the light
    % source
    color = (cos_reflAngle).*def_color.*(iter./200);
end
end

% newton function from lab sessions
function [X, n] = newton(F, JF, X0, tol, maxit)
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
end

function [cosf, r] = reflectionAngle(v, T, lightOrigin, dfdx, dfdy, dfdz)  
  % fi = reflectionAngle(f, T) returnes the angle 
  % needed for determining the color in raytracing and a reflection vector.  
  % Depeneding on the chosen option it can calculate angle:
  %    1. between normal vector and light source 
  %    2. between reflection vector of the ray and light source 

  % n is a normal to the plane from point T
  n = [feval(dfdx, T(1), T(2), T(3));  feval(dfdy, T(1), T(2), T(3)); feval(dfdz, T(1), T(2), T(3))];
  
  % normalizing n
  n = n./norm(n); 

  % computing the reflection vector
  vn = 2*dot(v, n);
  p = vn .* n;
  r = v - p;
  
  % computing the vector between the light origin and the intersection of
  % ray from the camera with the object
  lvec = lightOrigin - T;
    
  % first option
  cosf = dot(n, lvec) ./ (norm(n) .* norm(lvec));
  
  % second option 
  % cosf = dot(r, lvec) ./ (norm(r) .* norm(lvec));
end
