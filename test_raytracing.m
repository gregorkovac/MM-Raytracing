f = @(x, y, z) (x - 0).^2 + (y - 0.5).^2 + (z + 1).^2 - 0.3.^2;
% computing partial derivatives to the parametrization of the sphere
dfdx = @(x, y, z) 2.*x;
dfdy = @(x, y, z) 2.*(y - 0.5);
dfdz = @(x, y, z) 2.*(z + 1);

f1 = @(x, y, z) z + 2;
df1dx = @(x, y, z) 0;
df1dy = @(x, y, z) 0;
df1dz = @(x, y, z) 1;

%{
f1 = @(x, y, z) z + 1;
df1dx = @(x, y, z) 0;
df1dy = @(x, y, z) 0;
df1dz = @(x, y, z) 1;
%}
image_width = 800;
image_height = 600;
%initializing white image
image = zeros(image_height, image_width, 3); 

%initializing origin of rayes
lightOrigin = [0.5; 1; 0];
origin = [0; 0.5; 0];
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
    
    def_colors = [1 0 0; 0.9 0.9 0.9]';
    image(i, j, :) = raytracing(def_colors, f, f1, dfdx, dfdy, dfdz, df1dx, df1dy, df1dz, origin, direction_vector, lightOrigin, 0.01, 200, 10, 0);

    end
end

for i=2:image_height-1
    for j=2:image_width-1
        image(i, j, :) = 0.5.* image(i,j, :) + (image(i + 1, j, :) + image(i + 1, j + 1, :) + image(i + 1,j - 1, :) + image(i - 1, j, :) + image(i - 1, j + 1, :) + image(i - 1, j - 1, :) + image(i, j - 1, :) + image(i, j + 1))./16;
    end
end

imshow(image)