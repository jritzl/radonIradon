function output= radontransform(image_input,step_theta,num_samples)

 if isstruct(image_input)
        fieldNames = fieldnames(image_input);
        image = image_input.(fieldNames{1});
    elseif ischar(image_input) || isstring(image_input)
        [~, ~, ext] = fileparts(image_input);
        if strcmp(ext, '.mat')
            % Load the MAT file and extract the image data
            loaded_data = load(image_input);
            fieldNames = fieldnames(loaded_data);
            fieldNames = fieldNames(~startsWith(fieldNames, '__')); % Ignore MATLAB's internal variables
            image = loaded_data.(fieldNames{1}); % Use the first non-internal field
        else
            % Read the image file (jpg or png)
            image = imread(image_input);
        end
    else
        % If the input is already a matrix, use it directly
        image = image_input;
    end

    % Convert to 2D grayscale if not already
    if ndims(image) > 2
        image = rgb2gray(image);
    elseif ismatrix(image) == false
        image = squeeze(image);
    end

%step2
    % Extracting image dimensions and computing the diagonal length.
    [imgHeight, imgWidth] = size(image);
    maxDistance = sqrt(imgHeight^2 + imgWidth^2) / 2;
%step4
    % Generating vectors for theta and ts values based on the image dimensions.
    theta_values = (0:step_theta:180-step_theta); 
    t_values = linspace(-maxDistance, maxDistance, num_samples);
    gridX = linspace(-imgWidth / 2, imgWidth / 2, imgWidth+1);
    gridY = linspace(-imgHeight / 2, imgHeight / 2, imgHeight+1);

    % Initializing the matrix to store Radon transform results.
    output = zeros(length(t_values), length(theta_values));

    % Iteratively calculating Radon transform for each angle and t.
    for theta = 1:length(theta_values)
        current_angle = theta_values(theta);
        for t = 1:length(t_values)
            current_ts = t_values(t)+0.001; 
%step 5
            % Calculating intersection points for each ray and angle.
            y_intersects = (current_ts - gridX * cosd(current_angle)) / sind(current_angle);
            x_intersects = (current_ts - gridY * sind(current_angle)) / cosd(current_angle);

            % Combining and sorting the intersection points.
            all_intersections =[gridX' y_intersects' ; x_intersects' gridY'];
%step 7
           % Sort the intersection points by rows
            all_intersections = sortrows(all_intersections);
            % Find the unique rows (duplicate removal)
            unique_indices = [true; any(diff(all_intersections, 1, 1) ~= 0, 2)];
            % Keep only the unique rows
            all_intersections = all_intersections(unique_indices, :);
%step 6
            % Filtering intersections that fall outside the image.
           invalid_indices = all_intersections(:,1) < -imgWidth / 2 | all_intersections(:,1) > imgWidth / 2 | all_intersections(:,2) < -imgHeight / 2 | all_intersections(:,2) > imgHeight / 2;
           valid_intersections = all_intersections;
           valid_intersections(invalid_indices, :) = [];

%%I could not find a way that changing the order of step 7 and 6
            
            % Handling edge case with a single or zero intersection.
            if size(valid_intersections, 1) == 1 || size(valid_intersections, 1) == 0
                output(t, theta) = 0;
                continue
            end
%step 8
            % Calculate the differences between consecutive points for x and y coordinates
            dx = diff(valid_intersections(:,1));
            dy = diff(valid_intersections(:,2));

            % Calculate the distance (length) between consecutive points
            path_lengths = sqrt(dx.^2 + dy.^2);

            % Calculate the midpoints between consecutive points
              mid_x_points = (valid_intersections(1:end-1, 1) + valid_intersections(2:end, 1)) / 2;
              mid_y_points = (valid_intersections(1:end-1, 2) + valid_intersections(2:end, 2)) / 2;


            % Identifying the pixels intersected by the ray

            [Ygrid, Ymid] = meshgrid(gridY, mid_y_points);
            row_data = sum(Ygrid >= Ymid, 2);

            [Xgrid, Xmid] = meshgrid(gridX, mid_x_points);
            column_data = sum(Xgrid < Xmid, 2);
   
            % Gathering pixel values from the image.
            pixel_values = zeros(length(column_data), 1);
            for pixel_index = 1:length(column_data)
                if row_data(pixel_index) > 0 && row_data(pixel_index) <= imgHeight && column_data(pixel_index) > 0 && column_data(pixel_index) <= imgWidth
                    pixel_values(pixel_index) = image(row_data(pixel_index), column_data(pixel_index));
             
                end
                
            end

            % Calculating the projection value as a weighted sum.
            radon_val = sum(path_lengths .* pixel_values, 1);
            output(t, theta) = radon_val;
        end

    end


end


