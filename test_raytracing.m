% initializing sphere with the center in (0, 0, -1) and with radius 0.5
f = @(x, y, z) (x - 0).^2 + (y - 0).^2 + (z + 1).^2 - 0.5.^2;
% computing partial derivatives to the parametrization of the sphere
dfdx = @(x, y, z) 2.*x;
dfdy = @(x, y, z) 2.*y;
dfdz = @(x, y, z) 2.*(z + 1);


image_width = 400;
image_height = 300;
%initializing white image
image = ones(image_height, image_width, 3); 

%initializing origin of rayes
origin = [0; 0; 0];
origin_height = 2;
origin_width = origin_height.*image_width/image_height; 

%initializing vectors needed for computing direction vector
horizontal = [origin_width; 0; 0]
vertical = [0; origin_height; 0]
lower_left_corner = origin - horizontal./2 - vertical./2 - [0; 0; 1]

%getting color of every point of the matrix representing the image
for i=1:image_height
    for j=1:image_width
    u = j / image_width;
    v = i / image_height;

    direction_vector = lower_left_corner + u.*horizontal + v.*vertical - origin;
 
    image(i, j, :) = raytracing(f, dfdx, dfdy, dfdz, origin, direction_vector, 0.01, 100, 10);
    end
end
imshow(image)