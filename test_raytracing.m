% REFLECTIVE SPHERE
%{
% equation of sphere and corresponding partial derivatives
f = @(x, y, z) (x + 0).^2 + (y - 0.5).^2 + (z + 1).^2 - 0.5.^2;
dfdx = @(x, y, z) 2.*(x + 0);
dfdy = @(x, y, z) 2.*(y - 0.5);
dfdz = @(x, y, z) 2.*(z + 1);

% equation of plane and corresponding partial derivatives
f1 = @(x, y, z) 0.*x + -1.*y + 1.*z + 2;
df1dx = @(x, y, z) 0
df1dy = @(x, y, z) -1;
df1dz = @(x, y, z) 1;

% in this case, third object doesnt exsist
f2 = @(x, y, z) 0;
df2dx = @(x, y, z) 0;
df2dy = @(x, y, z) 0;
df2dz = @(x, y, z) 0;
%}

% TORUS SPHERE PLANE
% equation of sphere and corresponding partial derivatives
f = @(x, y, z) (x - 0).^2 + (y - 0.5).^2 + (z + 1).^2 - 0.2.^2;
dfdx = @(x, y, z) 2.*x;
dfdy = @(x, y, z) 2.*(y - 0.5);
dfdz = @(x, y, z) 2.*(z + 1);

% equation of sphere and corresponding partial derivatives
f1 = @(x, y, z) -2.*x - 2.*y + 5.*z + 7;
df1dx = @(x, y, z) -2;
df1dy = @(x, y, z) -2;
df1dz = @(x, y, z) 5;

% equation of torus and corresponding partial derivatives
f2 = @(x, y, z) (0.4 - ((x - 0).^2 + (y - 0.5).^2).^0.5).^2 + (z+0.7).^2 - 0.01;
df2dx = @(x, y, z) (x - 0) .* (2 - 0.8./((x - 0).^2 + (y - 0.5).^2).^0.5);
df2dy = @(x, y, z) (y - 0.5) .* (2 - 0.8./((x - 0).^2 + (y - 0.5).^2).^0.5);
df2dz = @(x, y, z) 2 .* (z+0.7);
%}

% additional objects
%{
% equation of cylinder and corresponding partial derivatives
f2 = @(x, y, z) (x + 0.2).^2 + (z+1).^2 - 0.1;
% computing partial derivatives to the parametrization of the sphere
df2dx = @(x, y, z) 2.*(x + 0.2);
df2dy = @(x, y, z) 0;
df2dz = @(x, y, z) 2.*(z + 1);

% equation of sphere and corresponding partial derivatives
%f2 = @(x, y, z) (x - 0.1).^2 + (y - 0.7).^2 + (z + 0.6).^2 - 0.1.^2;
%df2dx = @(x, y, z) 2.*(x - 0.1);
%df2dy = @(x, y, z) 2.*(y - 0.7);
%df2dz = @(x, y, z) 2.*(z + 0.6);
%}

image_width = 200;
image_height = 150;

%initializing white image
image = zeros(image_height, image_width, 3); 

%initializing origin of light and camera
lightOrigin = [0; 0; 0];
origin = [0; 0.5; 0];

origin_height = 2;
origin_width = origin_height.*image_width/image_height; 

%initializing vectors needed for computing direction vector
horizontal = [origin_width; 0; 0];
vertical = [0; origin_height; 0];
lower_left_corner = origin - horizontal./2 - vertical./2 - [0; 0; 1];

%getting color of every point of the matrix representing the image
for i=1:image_height
    for j=1:image_width
    u = j / image_width;
    v = i / image_height;
    
    % direction vector from the camera
    direction_vector = lower_left_corner + u.*horizontal + v.*vertical - origin;

    % default colors of objects, we want to draw on the image
    def_colors = [0.8 0 0.1; 0.5 0.5 0.5; 0.1 0.9 0.2]';
    
    % raytracing with camera
    image(i, j, :) = raytracing(f, f1, f2, dfdx, dfdy, dfdz, df1dx, df1dy, df1dz, df2dx, df2dy, df2dz, origin, direction_vector, lightOrigin, 0.01, 500, 0, def_colors, 0, 0, 0, [0.2; 0.8; 1]);
    end
end

% optional: softening the image
%for i=2:image_height-1
%    for j=2:image_width-1
%        image(i, j, :) = 0.5.* image(i,j, :) + (image(i + 1, j, :) + image(i + 1, j + 1, :) + image(i + 1,j - 1, :) + image(i - 1, j, :) + image(i - 1, j + 1, :) + image(i - 1, j - 1, :) + image(i, j - 1, :) + image(i, j + 1))./16;
%    end
%end

% showing image
imshow(image)