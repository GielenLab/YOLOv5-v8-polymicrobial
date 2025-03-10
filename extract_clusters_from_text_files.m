close all
clear all

base_folder_local='\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Poland\2025-01-22_Nanoparticles+T7\Drop 17\drop17';% folder with images

base_folder='\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Poland\2025-01-22_Nanoparticles+T7\Drop 17\exp17\labels'; % folder with txt detections

time_folder='\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Poland\2025-01-22_Nanoparticles+T7\Drop 17\drop17';

% Function to save nb_cells and time_tagged in an Excel file in the output folder
output_folder = '\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Poland\2025-01-22_Nanoparticles+T7\Drop 17'; % Define the output folder

nb_slices=24;

nb_time_points=17;

motion_range=5; %12.5

size_pic=640;
width=480;

%save movie
v=VideoWriter('bacteria_movie.mp4');
v.FrameRate=5;
open(v);
%%%%%%%%%%%%

lst_files = dir(fullfile(base_folder,'*txt')); % Gets a list of all npz files in folder

lst_files_im = dir(fullfile(base_folder_local,'*.jpeg')); % Gets a list of all tif files in folder

lst_files_time = dir(fullfile(time_folder,'*.jpeg')); % Gets a list of all tif files in folder

%natural sort order of list_files
[~, Index] = natsort({lst_files.name});
lst_files2   = lst_files(Index);

 [~, Index] = natsort({lst_files_im.name});
 lst_files_im2   = lst_files_im(Index);
 
  %[~, Index] = natsort({lst_files_time.name});
 lst_files_time2   = lst_files_time(Index);


number=zeros(length(lst_files2),1); % total number of cells
class_in=zeros(length(lst_files2),1); % class out-of-plane = 1 in class array

nb_cells=zeros(length(lst_files2)/nb_slices,1);

for t=1:nb_time_points %change to number of data points

t %display time point
    waitbar(t/nb_time_points);
    clf; %clear figure for each time point


    coord_center=[]; % store all cordinates for centers of bounding boxes for all slices

    for k = 1: nb_slices   % loop over each slice
      
        area=[];
        
        %meta=split(lst_files2(k).name,'_');
        full_file_name = fullfile(strcat(base_folder,'\foo',num2str(t-1),'_',num2str(k-1),'.txt'));   %lst_files2(k).name); % file name for detections
        full_file_name_im = fullfile(base_folder_local, lst_files_im2((t-1)*nb_slices+k).name); % corresponding name for image
       
              %%%FIGURE WITH RED DOTS AND FIRST SLICE
          figure(1);
          if k==14
           im=imread(full_file_name_im);
           imagesc(flip(im,1)); 
           set(gca,'Ydir','normal'); hold on;
           set(gca,'XTickLabel', []);
           set(gca,'YTickLabel', []);
          end
          
        %[nb, bboxes, class,scores]=read_npy(full_file_name);
        data=importdata(full_file_name);


        number(k)=size(data,1);
        class_in(k)=number(k);
        
        if number(k)~=0
            nb=number(k);
            %bboxes=round([data(:,2).*size_pic,data(:,3).*width,(data(:,2)+data(:,4)).*size_pic,(data(:,3)+data(:,5)).*width]);  
            bboxes=round([data(:,3).*width,data(:,2).*size_pic,(data(:,3)+data(:,5)).*width,(data(:,2)+data(:,4)).*size_pic]);  
            
            scores=data(:,6);
            class=data(:,1);
        
       %calculate area of all bounding boxes     
            for i=1:nb
                area=[area;  (bboxes(i,4)-bboxes(i,2)) * (bboxes(i,3)-bboxes(i,1))];
            end
            %%%%%%%%%%%%%%%%%
                
          
        
            
            for i=1:nb % nb is total number of beads detected
                
                %rectangle('Position',[bboxes(i,2) bboxes(i,1) bboxes(i,4)-bboxes(i,2) bboxes(i,3)-bboxes(i,1)]);
                
                %here choose whether to keep only in-plane cells
                if class(i)==0  %&& ((area(i)>150 && scores(i)>0.8) || area(i)>400) %|| class(i)==1% only keep in-plane cells = class 0   -- here triage based on area of bounding box 
                    
    
                    %coord_center : X - Y - CLASS - NB_SLICE - AREA -SCORE
                    coord_center=[coord_center;[bboxes(i,2)+(bboxes(i,4)-bboxes(i,2))/2 bboxes(i,1)+(bboxes(i,3)-bboxes(i,1))/2  class(i) k area(i) scores(i)]];
               
                end        
                
           
            end

        end
           
    
    end

%%%% get time tags
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the slice with most number of in-plane cells

    [a,b]=max(class_in);
     
    %now find bacteria that are close by and within motion_range
    dist=[];
    dist3=[];
    

count=0;

for l=1:size(coord_center,1) % put the coord of the max plane at the end of coord_center to process it first
    if coord_center(l,4)==b
        coord_center=[coord_center; coord_center(l,:)];
     end
end
coord_center=flipud(coord_center); %flip array
coord_center=unique(coord_center,'rows','stable'); %remove duplicates while keeping order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:size(coord_center,1)
    
        %find correspinding bin for center of bounding box and put 1 if one
        %cell and more if more
                   
        for j=1:size(coord_center,1)
        
            %check if many cells close on same plane

            dist(i,j)=sqrt((coord_center(i,1)-coord_center(j,1))^2+ (coord_center(i,2)-coord_center(j,2))^2); %calculate all distances
               
            % dist3(i,j)=sqrt((coord_center(i,1)-coord_center(j,1))^2+ (coord_center(i,2)-coord_center(j,2))^2 + (coord_center(i,4)-coord_center(j,4))^2); %calculate all distances
    
% %             %plot all cells on plane with max number of cells first
            if coord_center(i,4)==b
                % plot(coord_center(i,1),250-coord_center(i,2),'o','MarkerSize',8,'MarkerFaceColor',[1  0 0]); hold on; %plot detected cell in red
                 %scatter3(coord_center(i,1),250-coord_center(i,2),b*0.5,'o','MarkerFaceColor',[1  0 0]);hold on; %in 3D
                %plot(coord_center(i,1),b*0.5,'o','MarkerSize',8,'MarkerFaceColor',[1  0 0]); hold on; %plot detected cell in red
                 
                % xlim([0 250]);
                % ylim([0 250]);
                % zlim([0 10]);
%                  %DRAW EXCLUSION CIRCLE
%                 c=[coord_center(i,1) 250-coord_center(i,2)]; 
%                 pos = [c-motion_range 2*motion_range 2*motion_range];
%                 rectangle('Position',pos,'Curvature',[1 1]);
           
            end


        if dist(i,j)<motion_range &&  i~=j && coord_center(i,4)~=coord_center(j,4) %&&  i~=j% if two points are close and not on same plane set RED coord_center(i,4)~=coord_center(j,4) 
         
      
                  % plot(coord_center(i,1),250-coord_center(i,2),'o','MarkerSize',8,'MarkerFaceColor',[1  0 0]); hold on; %plot detected cell in red
                   
                  % scatter3(coord_center(i,1),250-coord_center(i,2),b*0.5,'o','MarkerFaceColor',[1  0 0]);hold on;%in 3D
                  % plot(coord_center(i,1),coord_center(i,4)*0.5,'o','MarkerSize',8,'MarkerFaceColor',[1  0 0]); hold on; %plot detected cell in red
                
                   %remove point close by from coord_center
                  
                   coord_center(j,:)=NaN; %coord_center(j,:)=NaN;

                 
         end
       
        end

    end

if ~isempty(coord_center)

    coord_1=coord_center(:,1);
    coord_1(isnan(coord_1))=[];

    nb_cells(t)=length(coord_1);
    
    %plot all remaining detected cells as red dots
    for i=1:size(coord_center,1)
        hold on; 
        plot(coord_center(i,1),width-coord_center(i,2),'o','MarkerSize',8,'MarkerFaceColor',[1  0 0]); hold on; %plot detected cell in red
    end
    xlim([0 size_pic]);
    ylim([0 width]);

end

    text(0.8*size_pic,0.8*width,num2str(nb_cells(t)),'FontSize',18);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    next_frame=getframe(gcf);
    writeVideo(v,next_frame);
    
% 
    max_class(t)=b; %slice where number of in--plane cells is max
     
       
    
     

end

close(v); %close video

figure(2);
plot(nb_cells);

%plot(movmean(max_class,10));

% Convert time_tagged to minutes
time_minutes = time_tag / 60;  % Convert time_tagged to minutes

% Ensure time_minutes and nb_cells are of the same length
len_time = length(time_minutes);
len_cells = length(nb_cells);

% Adjust lengths if necessary
if len_time > len_cells
    time_minutes = time_minutes(1:len_cells);  % Truncate time_minutes if it's longer
elseif len_cells > len_time
    nb_cells = nb_cells(1:len_time);  % Truncate nb_cells if it's longer
end

% Combine nb_cells and time_minutes into one matrix
data_to_save = [time_minutes', nb_cells];  % Transpose time_minutes and combine with nb_cells

% Define output folder and create it if necessary
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Define the full path for the Excel file
filename = fullfile(output_folder, 'bacteria_data.xlsx');

% Write headers and data to Excel
headers = {'Time (minutes)', 'Number of Cells'};
xlswrite(filename, headers, 'Sheet1', 'A1');  % Write headers
xlswrite(filename, data_to_save, 'Sheet1', 'A2');  % Write data starting from the second row

disp(['Data saved to ', filename]);

%%%%%%%%%%%%%%%%% HOW TO READ NUMPIES EXAMPLE %%%%%%%%%%%%%%%%%%%%%
function [nz,bboxes,class,scores] =read_npy(file_name)
   
    a=unzip(file_name); 

    b=readNPY(a{1}); %read numpies 
    scores=readNPY(a{2});
    b=squeeze(b); %emove first dimension which is not used
    nz=nnz(b)/4 ;%number of non zero elements,divide by 4 coordinates = number of bacteria in total
    class=readNPY(a{3});
    bboxes=b;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%

