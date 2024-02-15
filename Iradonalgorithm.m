function Image = Iradonalgorithm(Row, Column, projection_data, filter_type, window_size)

    switch lower(filter_type)
         case 'none'
            projectiondata =projection_data; % Triangular filter
        otherwise

            num_beams = size(projection_data, 1); % Number of beams
            num_angles = size(projection_data, 2); % Number of angles
        
            % Validate and adjust the window size if necessary
            window_size = min(window_size, num_beams);
        
            % Initialize the filtered projection data matrix
            filtered_projection = zeros(size(projection_data));
        
            % Step 1: Design the filter
            n = linspace(-window_size/2, window_size/2, window_size);
        
            switch lower(filter_type)
                case 'triangle'
                    filter = abs(n / (window_size / 2)); % Triangular filter
                case 'ramlak'
                    filter = abs(n);
                case 'ramlak-hamming'
                    % Define basic Ram-Lak filter
                    filter = [0, 1:(window_size/2), (window_size/2-1):-1:1];
                    % Apply Hamming window to Ram-Lak filter
                    w_hamming = window(@hamming, window_size);
                    filter = filter .* fftshift(w_hamming');
                    filter = filter / max(filter); % Normalize the filter
                case 'ramlak-barthann'
                    % Define basic Ram-Lak filter
                    filter = [0, 1:(window_size/2), (window_size/2-1):-1:1];
                    % Apply  window to Ram-Lak filter
                    w_barthann = window(@barthannwin, window_size);
                    filter = filter .* fftshift(w_barthann');
                    filter = filter / max(filter); % Normalize the filter
                case 'ramlak-rectwin'
                    % Define basic Ram-Lak filter
                    filter = [0, 1:(window_size/2), (window_size/2-1):-1:1];
                    % Apply  window to Ram-Lak filter
                    w_barthann = window(@rectwin, window_size);
                    filter = filter .* fftshift(w_barthann');
                    filter = filter / max(filter); % Normalize the filter

                case 'ramlak-blackman'
                    % Define basic Ram-Lak filter
                    filter = [0, 1:(window_size/2), (window_size/2-1):-1:1];
                    % Apply  window to Ram-Lak filter
                    w_barthann = window(@blackman, window_size);
                    filter = filter .* fftshift(w_barthann');
                    filter = filter / max(filter); % Normalize the filter


                case 'ramlak-hann'
                    % Define basic Ram-Lak filter
                    filter = [0, 1:(window_size/2), (window_size/2-1):-1:1];
                    % Apply  window to Ram-Lak filter
                    w_barthann = window(@hann, window_size);
                    filter = filter .* fftshift(w_barthann');
                    filter = filter / max(filter); % Normalize the filter

                case 'ramlak-nuttallwin'
                    % Define basic Ram-Lak filter
                    filter = [0, 1:(window_size/2), (window_size/2-1):-1:1];
                    % Apply  window to Ram-Lak filter
                    w_barthann = window(@nuttallwin, window_size);
                    filter = filter .* fftshift(w_barthann');
                    filter = filter / max(filter); % Normalize the filter


                case 'ramlak-gausswin'
                    % Define basic Ram-Lak filter
                    filter = [0, 1:(window_size/2), (window_size/2-1):-1:1];
                    % Apply  window to Ram-Lak filter
                    w_barthann = window(@gausswin, window_size);
                    filter = filter .* fftshift(w_barthann');
                    filter = filter / max(filter); % Normalize the filter   
                    
                case 'none'
                      filter = abs(n / (window_size / 2));
                otherwise
                    error('Unknown filter type');
            end
        
            % Pad the filter to the size of the number of beams if necessary
            if length(filter) < num_beams
                padding = (num_beams - length(filter)) / 2;
                filter = [zeros(1, floor(padding)), filter, zeros(1, ceil(padding))];
            end
        
            filter = fftshift(filter); % Center the filter
        
            % Apply the filter for each angle
            for angle = 1:num_angles
                % Step 2: Obtain the DFT of the projection data
                DFT_data = fft(projection_data(:, angle));
        
                % Step 3: Apply the filter to the DFT
                filtered_DFT = DFT_data .* filter';
        
                % Step 4: Take the inverse DFT
                filtered_projection(:, angle) = ifft(filtered_DFT, 'symmetric');
            end
                 projectiondata =filtered_projection;
      end

    [num_samples, angle] = size(projectiondata);
    step_theta=180/angle;

%step2
    % Extracting image dimensions and computing the diagonal length.
    imgHeight = Row;
    imgWidth = Column;
    Image = zeros(imgHeight, imgWidth);

    maxDistance = sqrt(imgHeight^2 + imgWidth^2) / 2;
%step4
    % Generating vectors for theta and ts values based on the image dimensions.
    theta_values = (0:step_theta:180-step_theta); 
    t_values = linspace(-maxDistance, maxDistance, num_samples);
    gridX = linspace(-imgWidth / 2, imgWidth / 2, imgWidth+1);
    gridY = linspace(-imgHeight / 2, imgHeight / 2, imgHeight+1);



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
            
            
            for i = 1:length(column_data)
                if row_data(i) > 0 && row_data(i) <= imgHeight && column_data(i) > 0 && column_data(i) <= imgWidth
                     Image(row_data(i), column_data(i)) = Image(row_data(i), column_data(i)) + path_lengths(i) * projectiondata(t, theta);
                 end
            end
                
          end

      end


end
