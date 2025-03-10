% Initialize
clear all;

base_folder_local = 'F:\DATA ANALYSIS\2024-08-19-PA+SA+p278\Drop 15\Masked'; % Folder with images
base_folder = 'F:\DATA ANALYSIS\2024-08-19-PA+SA+p278\Drop 15\exp5\labels'; % Folder with txt detections
time_folder = 'F:\DATA ANALYSIS\2024-08-19-PA+SA+p278\Drop 15\drop15';

nb_slices = 24; % Number of slices acquired per time point
nb_time_points = 350; % Total number of time points
motion_range = 5; % Determines exclusion radius in pixel number
size_pic = 640;
width = 480;
score_threshold = 0.4;

% Initialize VideoWriter
v = VideoWriter('bacteria_movie_PASA_graph.mp4');
v.FrameRate = 12;
open(v);

% Get list of files
lst_files = dir(fullfile(base_folder, '*.txt')); % Gets a list of all txt files in folder
lst_files_im = dir(fullfile(base_folder_local, '*.jpg')); % Gets a list of all jpg files in folder

% Natural sort order of list_files
[~, Index] = natsort({lst_files.name});
lst_files2 = lst_files(Index);

[~, Index] = natsort({lst_files_im.name});
lst_files_im2 = lst_files_im(Index);

number = zeros(length(lst_files2), 1); % Total number of cells
class_out = zeros(length(lst_files2), 1); % Class out-of-plane
class_in = zeros(length(lst_files2), 1); % Class in-plane
doublets = zeros(round(length(lst_files2) / nb_slices), 1);
nb_cells = zeros(round(length(lst_files2) / nb_slices), 1);
time_tag = zeros(nb_time_points, 1); % Array to hold the time tags
cell_identity = [];

% Get time tags for each time point
FileInfoini = fullfile(time_folder, 'foo0_0.jpeg'); % Function to extract time point of the first image
time_t0 = GetFileTime(FileInfoini); % First image's timestamp
time_tag0 = (time_t0.Write(4)*3600 + time_t0.Write(5)*60 + time_t0.Write(6)) + (time_t0.Write(3)-1)*86400; % time-point of the first image in seconds
for i = 1:nb_time_points % Extract time-points of all subsequent foo_i_0 images 
    FileInfo = fullfile(time_folder, sprintf('foo%d_0.jpeg', i-1));
    time_t = GetFileTime(FileInfo);
    time_tag(i)= (time_t.Write(4)*3600 + time_t.Write(5)*60 + time_t.Write(6)) + (time_t.Write(3)-1)*86400 - time_tag0;
end

% Initialize figure
figure('Position', [100, 100, 1200, 400]);

% Loop through each time point
for t = 1:nb_time_points
    t % Display time point
    waitbar(t/nb_time_points);
    clf; % Clear figure for each time point

    coord_center = []; % Store all coordinates for centers of bounding boxes for all slices

    for k = 1:nb_slices
        area = [];
        full_file_name = fullfile(base_folder, sprintf('foo%d_%d.txt', t-1, k-1));
        full_file_name_im = fullfile(base_folder_local, lst_files_im2((t-1)*nb_slices+k).name);

        data = importdata(full_file_name);
        number(k) = size(data, 1);
        class_in(k) = number(k);

        if number(k) ~= 0
            nb = number(k);
            bboxes = round([data(:,3).*width, data(:,2).*size_pic, (data(:,3)+data(:,5)).*width, (data(:,2)+data(:,4)).*size_pic]);
            scores = data(:,6);
            class = data(:,1);

            for i = 1:nb
                area = [area; [(bboxes(i,4)-bboxes(i,2)) * (bboxes(i,3)-bboxes(i,1))]];
                CENTER_X = bboxes(i,2) + (bboxes(i,4) - bboxes(i,2)) / 2;
                CENTER_Y = bboxes(i,1) + (bboxes(i,3) - bboxes(i,1)) / 2;
                coord_center = [coord_center; [CENTER_X CENTER_Y class(i) k area(i) scores(i)]];
            end
        end
    end

    coord_center_2 = coord_center; % Keep a copy of coord_center / coord_center_2 keeps only final red dots coords

    [a, b] = max(class_in); % Find the slice with the most number of in-plane cells

    % Find bacteria that are close by and within motion_range
    dist = [];
    count = 0;

    for l = 1:size(coord_center, 1)
        if coord_center(l, 4) == b
            coord_center = [coord_center; coord_center(l, :)];
        end
    end
    coord_center = flipud(coord_center);
    coord_center = unique(coord_center, 'rows', 'stable');

    for i = 1:size(coord_center_2, 1)
        for j = 1:size(coord_center_2, 1)
            dist(i, j) = sqrt((coord_center_2(i, 1) - coord_center_2(j, 1))^2 + (coord_center_2(i, 2) - coord_center_2(j, 2))^2);
            if dist(i, j) < motion_range && i ~= j && coord_center_2(i, 4) ~= coord_center_2(j, 4)
                coord_center_2(j, :) = NaN;
            end
        end
    end

    coord_center_2(any(isnan(coord_center_2), 2), :) = [];

    nb_cells_total(t) = size(coord_center_2, 1); % Total number of cells in time point

    if ~isempty(coord_center_2)
        dual_detect_class = zeros(nb_cells_total(t), 2);
        for i = 1:size(coord_center_2, 1)
            for j = 1:size(coord_center, 1)
                dist = sqrt((coord_center_2(i, 1) - coord_center(j, 1))^2 + (coord_center_2(i, 2) - coord_center(j, 2))^2);
                if dist < motion_range && i ~= j && coord_center(i, 4) ~= coord_center(j, 4)
                    if coord_center(j, 3) == 0 && coord_center(j, 6) > score_threshold
                        dual_detect_class(i, 1) = dual_detect_class(i, 1) + 1;
                    elseif coord_center(j, 3) == 1 && coord_center(j, 6) > score_threshold
                        dual_detect_class(i, 2) = dual_detect_class(i, 2) + 1;
                    end
                end
            end

            cell_identity_local = dual_detect_class(:, 1) ./ dual_detect_class(:, 2);
            cell_identity = [cell_identity; cell_identity_local];
            cell_identity(cell_identity == Inf) = 5;

            nb_pseudo(t) = sum(cell_identity_local >= 1);
            nb_staph(t) = sum(cell_identity_local <= 1);
        end
    end

    % Load the image with the most detections and convert to grayscale
    [a, b] = max(class_in); % Slice with the most number of in-plane cells
    full_file_name_im = fullfile(base_folder_local, lst_files_im2((t-1)*nb_slices+b).name);
    best_image = imread(full_file_name_im);
    
    % Convert to grayscale if the image is RGB or already grayscale
    if size(best_image, 3) == 3
        best_image_gray = rgb2gray(best_image); % Convert RGB to grayscale
    else
        best_image_gray = best_image; % Image is already grayscale
    end
    
    % Left subplot: Display the image
    subplot(1, 2, 1);
    imagesc(flip(best_image_gray, 1));
    colormap gray; % Set colormap to grayscale
    set(gca, 'YDir', 'normal');
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    xlim([0 size_pic]);
    ylim([0 width]);
    hold on;

    % Plot the detected cells
    for i = 1:size(coord_center_2, 1)
        if coord_center_2(i, 3) == 0 % PA cells
            plot(coord_center_2(i, 1), width - coord_center_2(i, 2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'black'); % Green
                elseif coord_center_2(i, 3) == 1 % SA cells
            plot(coord_center_2(i, 1), width - coord_center_2(i, 2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black'); % Red
        end
    end

    % Add timestamp and cell counts as text
    elapsed_time = time_tag(t)/60; % Time in minutes
    text(0.05 * size_pic, 0.95 * width, sprintf('Time: %.1f min', elapsed_time), 'Color', 'white', 'FontSize', 15); % Display time in minutes
    text(0.05 * size_pic, 0.87 * width, sprintf('PA Cells: %d', nb_pseudo(t)), 'Color', 'white', 'FontSize', 15); % Display number of Pseudomonas cells
    text(0.05 * size_pic, 0.80 * width, sprintf('SA Cells: %d', nb_staph(t)), 'Color', 'white', 'FontSize', 15); % Display number of Staph cells

    % Right subplot: Plot the cell counts
    subplot(1, 2, 2);
    plot(time_tag(1:t) / 60, nb_cells_total(1:t), 'k', 'DisplayName', 'Total');
    hold on;
    plot(time_tag(1:t) / 60, nb_pseudo(1:t), 'g', 'DisplayName', 'Pseudomonas');
    plot(time_tag(1:t) / 60, nb_staph(1:t), 'r', 'DisplayName', 'Staphylococcus');
    legend('Total', 'Pseudomonas', 'Staphylococcus');
    xlabel('Time (min)');
    ylabel('Count');
    title('Cell Counts Over Time');
    grid on;
    xlim([0 max(time_tag) / 60]);
    
    % Capture the frame and write it to the video
    next_frame = getframe(gcf);
    writeVideo(v, next_frame);
end

% Close the video file
close(v);
disp('Video saved successfully!');

% Save results to .mat file
% save(fullfile(output_folder, 'results.mat'), 'nb_cells_total', 'nb_pseudo', 'nb_staph', 'time_tags');
