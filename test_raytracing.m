% initializing sphere with the center in (0, 0, -1) and with radius 0.5
%{
f = @(x, y, z) 1 + x + y + z;
dfdx = @(x, y, z) 1;
dfdy = @(x, y, z) 1;
dfdz = @(x, y, z) 1;
%}
f = @(x, y, z) (x - 0).^2 + (y - 0).^2 + (z + 1).^2 - 0.5.^2;
% computing partial derivatives to the parametrization of the sphere
dfdx = @(x, y, z) 2.*x;
dfdy = @(x, y, z) 2.*y;
dfdz = @(x, y, z) 2.*(z + 1);
%}

image_width = 400;
image_height = 300;
%initializing white image
image = ones(image_height, image_width, 3); 

% second photo
image2 = zeros(image_height, image_width, 3);

%initializing origin of rayes
origin = [0.3; 0.3; 0.3];
origin_height = 2;
origin_width = origin_height.*image_width/image_height; 

%initalizing position of the camera
origin_camera = origin + [0.2; 0.2; 0.2];

%initializing vectors needed for computing direction vector
horizontal = [origin_width; 0; 0]
vertical = [0; origin_height; 0]
lower_left_corner = origin - horizontal./2 - vertical./2 - [0; 0; 1]

%imshow(image2)
%getting color of every point of the matrix representing the image
for i=1:image_height
    for j=1:image_width
    u = j / image_width;
    v = i / image_height;

    direction_vector = lower_left_corner + u.*horizontal + v.*vertical - origin;
    %direction_vector = image_point - origin_camera;
   
    
    [color, intersection] = raytracing(f, dfdx, dfdy, dfdz, origin, direction_vector, 0.01, 100, 10);
    %image2(i, j, :) = raytracing(f, dfdx, dfdy, dfdz, origin_camera, direction_vector, 0.01, 100, 10);
    if(~isempty(intersection))% & intersection(1, :) > 0 & intersection(2, :) > 0)
        %intersection = round(intersection + [1 1]);
       
        x_idx = image_width .* intersection(1, 1) ./ ((origin_camera(1, :) + 0.01.*image_width./2)) + image_width ./ 2;
        y_idx = image_height .* intersection(2, 1) ./ ((origin_camera(2, :) + 0.01.*image_height./2)) + image_height ./ 2;
        
        if(round(abs(x_idx)) ~= 0 & round(abs(y_idx)) ~= 0)
            if(image2(floor(x_idx), floor(y_idx), :) == [0;0;0])
                image2(floor(x_idx), floor(y_idx), :) = color;
            elseif(image2(floor(x_idx), ceil(y_idx), :) == [0;0;0])
                image2(floor(x_idx), ceil(y_idx), :) = color;
            elseif(image2(ceil(x_idx), floor(y_idx), :) == [0;0;0])
                image2(ceil(x_idx), floor(y_idx), :) = color;
            else image2(ceil((x_idx)), ceil((y_idx)), :) = color;
            end
        end
 %       image2(abs(round(intersection(1, :).*100)), abs(round(intersection(2, :).*100)), :) = color;
    end    
   % image2(inse, j, :) = (cos_reflAngle).*[1; 1; 1];  
    end
end
imshow(image2)