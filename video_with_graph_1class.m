close all
clear all

base_folder_local = '\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Poland\2025-01-22_Nanoparticles+T7\Drop 14\drop14'; % folder with images
base_folder = '\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Poland\2025-01-22_Nanoparticles+T7\Drop 14\exp14\labels'; % folder with txt detections
time_folder = '\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Poland\2025-01-22_Nanoparticles+T7\Drop 14\drop14';

nb_slices = 24;
nb_time_points = 17;
motion_range = 5;

size_pic = 640;
width = 480;

% Save movie
v = VideoWriter('bacteria_movie_T7+np_lysis_highMOI.mp4');
v.FrameRate = 12;
open(v);

%%%%%%%%%%%%

lst_files = dir(fullfile(base_folder, '*txt')); % Gets a list of all txt files in folder
lst_files_im = dir(fullfile(base_folder_local, '*.jpeg')); % Gets a list of all jpg files in folder
lst_files_time = dir(fullfile(time_folder, '*.jpeg')); % Gets a list of all jpeg files in folder

% Natural sort order of list_files
[~, Index] = natsort({lst_files.name});
lst_files2 = lst_files(Index);

[~, Index] = natsort({lst_files_im.name});
lst_files_im2 = lst_files_im(Index);

[~, Index] = natsort({lst_files_time.name});
lst_files_time2 = lst_files_time(Index);

number = zeros(length(lst_files2), 1); % total number of cells
class_in = zeros(length(lst_files2), 1); % class out-of-plane = 1 in class array
nb_cells = zeros(nb_time_points, 1); % array to hold the number of cells
time_tag = zeros(nb_time_points, 1); % array to hold the time tags

% Get time tags for each time point
FileInfoini = fullfile(time_folder, 'foo0_0.jpeg'); % Function to extract time point of the first image
time_t0 = GetFileTime(FileInfoini); % First image's timestamp
time_t0 = GetFileTime(FileInfoini);
time_tag0 = (time_t0.Write(4)*3600 + time_t0.Write(5)*60 + time_t0.Write(6)) + (time_t0.Write(3)-1)*86400; % time-point of the first image in seconds
for i = 1:nb_time_points % Extract time-points of all subsequent foo_i_0 images 
    FileInfo = fullfile(time_folder, sprintf('foo%d_0.jpeg', i-1));
    time_t = GetFileTime(FileInfo);
    time_t= GetFileTime(FileInfo);
    time_tag(i)= (time_t.Write(4)*3600 + time_t.Write(5)*60 + time_t.Write(6)) + (time_t.Write(3)-1)*86400 - time_tag0;
end

% Initialize a figure with two side-by-side subplots
figure;

set(gcf, 'Position', [100, 100, 1200, 400]); % Set the figure position and size

for t = 1:nb_time_points % change to number of data points
    t % display time point
    waitbar(t/nb_time_points);
    clf; % clear figure for each time point

    coord_center = []; % store all coordinates for centers of bounding boxes for all slices

    for k = 1:nb_slices % loop over each slice

        area = [];

        full_file_name = fullfile(strcat(base_folder, '\foo', num2str(t-1), '_', num2str(k-1), '.txt'));
        full_file_name_im = fullfile(base_folder_local, lst_files_im2((t-1)*nb_slices+k).name); % corresponding name for image

        data = importdata(full_file_name);

        number(k) = size(data, 1);
        class_in(k) = number(k);

        if number(k) ~= 0
            nb = number(k);
            bboxes = round([data(:, 3).*width, data(:, 2).*size_pic, (data(:, 3)+data(:, 5)).*width, (data(:, 2)+data(:, 4)).*size_pic]);
            scores = data(:, 6);
            class = data(:, 1);

            % Calculate area of all bounding boxes     
            for i = 1:nb
                area = [area; (bboxes(i, 4)-bboxes(i, 2)) * (bboxes(i, 3)-bboxes(i, 1))];
            end

            for i = 1:nb
                if class(i) == 0
                    coord_center = [coord_center; [bboxes(i, 2)+(bboxes(i, 4)-bboxes(i, 2))/2, bboxes(i, 1)+(bboxes(i, 3)-bboxes(i, 1))/2, class(i), k, area(i), scores(i)]];
                end        
            end
        end
    end

    % Find the slice with most number of in-plane cells
    [~, b] = max(class_in);

    % Subplot 1: Display the image corresponding to the slice with most detections
    subplot(1, 2, 1);
    full_file_name_im_max = fullfile(base_folder_local, lst_files_im2((t-1)*nb_slices + b).name); % Image with max detections
    im = imread(full_file_name_im_max);
    
    % Check if the image is RGB (3 channels), then convert to grayscale
    if size(im, 3) == 3
        im = rgb2gray(im);
    end
    
    imagesc(flip(im, 1)); % Display the flipped grayscale image
    colormap gray; % Apply grayscale colormap
    set(gca, 'Ydir', 'normal'); 
    hold on;
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);

    dist = [];
    % Process coord_center array
    for l = 1:size(coord_center, 1)
        if coord_center(l, 4) == b
            coord_center = [coord_center; coord_center(l, :)];
        end
    end

    coord_center = flipud(coord_center);
    coord_center = unique(coord_center, 'rows', 'stable');

    for i = 1:size(coord_center, 1)
        for j = 1:size(coord_center, 1)
            dist(i, j) = sqrt((coord_center(i, 1)-coord_center(j, 1))^2 + (coord_center(i, 2)-coord_center(j, 2))^2);

            if dist(i, j) < motion_range && i ~= j && coord_center(i, 4) ~= coord_center(j, 4)
                coord_center(j, :) = NaN;
            end
        end
    end

    if ~isempty(coord_center)
        coord_1 = coord_center(:, 1);
        coord_1(isnan(coord_1)) = [];

        nb_cells(t) = length(coord_1);

        for i = 1:size(coord_center, 1)
            hold on; 
            plot(coord_center(i, 1), width-coord_center(i, 2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'black'); % plot detected cell in red
        end
        xlim([0 size_pic]);
        ylim([0 width]);
    end

    % Display elapsed time and number of cells on the frame
    elapsed_time = time_tag(t)/60; % Time in minutes
    num_cells = nb_cells(t); % Number of cells
    text(0.1*size_pic, 0.9*width, sprintf('Time: %.1f min', elapsed_time), 'Color', 'black', 'FontSize', 14); % Display time in minutes
    text(0.1*size_pic, 0.85*width, sprintf('Cells: %d', num_cells), 'Color', 'black', 'FontSize', 14); % Display number of cells

    % Subplot 2: Real-time plot of the number of cells vs time
    subplot(1, 2, 2);
    plot(time_tag(1:t)/60, nb_cells(1:t), '-', 'Color', 'green', 'MarkerFaceColor', 'green');
    xlabel('Time (minutes)');
    ylabel('Number of cells');
    title('Real-time Cell Count');
    xlim([0 max(time_tag)/60]);
    ylim([0 max(nb_cells) + 5]); % Add some padding for better visualization
    grid on;

    % Capture and write frame to video
    next_frame = getframe(gcf);
    writeVideo(v, next_frame);

    max_class(t) = b; % slice where number of in-plane cells is max
end

close(v); % close video

figure(2);
plot(nb_cells);

save('nb_cells.mat'); % saves number of cell array for all time points
