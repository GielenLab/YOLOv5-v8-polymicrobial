%lose all
clear all

base_folder_local='\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Robyn\PW Lysis\2024-09-10_lowMOI_PA\Drop 15';% folder with images

base_folder='\\mcrlsinas01.ex.ac.uk\lab-gielen\Anuj\Robyn\PW Lysis\2024-09-10_lowMOI_PA\Drop 15'; % folder with txt detections

time_folder='F:\DATA ANALYSIS\2024-08-11_PA+SA_Growth\Drop 26\drop26';

output_folder= 'F:\DATA ANALYSIS\2024-08-11_PA+SA_Growth\Drop 26\New Model'

nb_slices=24; % number of slices acquired per time point

nb_time_points=60 % total number of time points

motion_range=5; %12.5 % determines exclusion radius in pixel number

size_pic=640;
width=480;

score_threshold=0.4;

 save movie
 v=VideoWriter('bacteria_movie_PASA.mp4');
 v.FrameRate=5;
 open(v);
 %%%%%%%%%%%%

lst_files = dir(fullfile(base_folder,'*.txt')); % Gets a list of all npz files in folder

lst_files_im = dir(fullfile(base_folder_local,'*.jpg')); % Gets a list of all tif files in folder

%lst_files_time = dir(fullfile(time_folder,'*.jpeg')); % Gets a list of all tif files in folder

%lst_files_time = dir(fullfile(time_folder,'*.jpeg')); % Gets a list of all tif files in folder
%natural sort order of list_files
[~, Index] = natsort({lst_files.name});
lst_files2   = lst_files(Index);

 [~, Index] = natsort({lst_files_im.name});
 lst_files_im2   = lst_files_im(Index);





number=zeros(length(lst_files2),1); % total number of cells
class_out=zeros(length(lst_files2),1); % class out-of-plane = 1 in class array
class_in=zeros(length(lst_files2),1); % class out-of-plane = 1 in class array

doublets=zeros(round(length(lst_files2)/nb_slices),1);
nb_cells=zeros(round(length(lst_files2)/nb_slices),1);

cell_identity=[];

for t=1: nb_time_points %change to number of data points

t %display time point
    waitbar(t/nb_time_points);
    clf; %clear figure for each time point


    coord_center=[]; % store all cordinates for centers of bounding boxes for all slices

    for k = 1:nb_slices   % loop over each slice
      
        area=[];
        
        %meta=split(lst_files2(k).name,'_');
        full_file_name = fullfile(strcat(base_folder,'\foo',num2str(t-1),'_',num2str(k-1),'.txt'));   %lst_files2(k).name); % file name for detections
       
        full_file_name_im = fullfile(base_folder_local, lst_files_im2((t-1)*nb_slices+k).name); % corresponding name for image
       
    
        %[nb, bboxes, class,scores]=read_npy(full_file_name);
        data=importdata(full_file_name);

        number(k)=size(data,1);
        class_in(k)=number(k);

        if number(k)~=0
            nb=number(k);
            
            bboxes=round([data(:,3).*width,data(:,2).*size_pic,(data(:,3)+data(:,5)).*width,(data(:,2)+data(:,4)).*size_pic]);  
            
            scores=data(:,6);
            class=data(:,1);
        
        
       %calculate area of all bounding boxes     
            for i=1:nb
                area=[area; [ (bboxes(i,4)-bboxes(i,2)) * (bboxes(i,3)-bboxes(i,1))]];
            end
            %%%%%%%%%%%%%%%%%
                
                %%%FIGURE WITH RED DOTS AND FIRST SLICE
           %figure(1);
           %if k==10 %plot only the k slice (best one)
            %a=imread(full_file_name_im);
            %imagesc(flip(a,1)); 
            %set(gca,'Ydir','normal'); hold on;
            %set(gca,'XTickLabel', []);
            %set(gca,'YTickLabel', []);
          %xlim([0 size_pic]);
          %ylim([0 width]);

          % end
        
            
            for i=1:nb % nb is total number of beads detected
                               
                %here choose whether to keep only in-plane cells
                
                    %coord_center : X - Y - CLASS - NB_SLICE - AREA -SCORE
                    CENTER_X=bboxes(i,2)+(bboxes(i,4)-bboxes(i,2))/2;
                    CENTER_Y=bboxes(i,1)+(bboxes(i,3)-bboxes(i,1))/2;
                    coord_center=[coord_center;[CENTER_X CENTER_Y  class(i) k area(i) scores(i)]];
               
                    %plot dot for all detections
                    %plot(CENTER_X,width-CENTER_Y,'o','MarkerSize',8,'MarkerFaceColor',[0 class(i) (1-class(i))]); hold on;
          
            end
       
        end

    end


coord_center_2=coord_center; %keep a copy of coord_center / coord_center_2 keeps only final red dots coords

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

    for i=1:size(coord_center_2,1)
    
        %find correspinding bin for center of bounding box and put 1 if one
        %cell and more if more
                   
        for j=1:size(coord_center_2,1)
        
            %check if many cells close on same plane

            dist(i,j)=sqrt((coord_center_2(i,1)-coord_center_2(j,1))^2+ (coord_center_2(i,2)-coord_center_2(j,2))^2); %calculate all distances
               
            if dist(i,j)<motion_range &&  i~=j && coord_center_2(i,4)~=coord_center_2(j,4) %&&  i~=j% if two points are close and not on same plane set RED coord_center(i,4)~=coord_center(j,4) 
         
                   
                  % plot(coord_center(i,1),250-coord_center(i,2),'o','MarkerSize',8,'MarkerFaceColor',[1  0 0]); hold on; %plot detected cell in red
                           
                  
                   coord_center_2(j,:)=NaN; %remove point close by from coord_center

            end
       
        end

    end

coord_center_2(~any(~isnan(coord_center_2), 2),:)=[]; %remove NaN from list

nb_cells_total(t)=size(coord_center_2,1); %total number of cells in time point


if ~isempty(coord_center_2)
   

    dual_detect_class=zeros(nb_cells_total(t),2);
    %plot all remaining detected cells as red dots
    for i=1:size(coord_center_2,1)
    
        %build a table with [number Pseusomonas ; number Staph ; average score Pseudomonas; average score Staph]
        for j=1:size(coord_center,1)

         dist=sqrt((coord_center_2(i,1)-coord_center(j,1))^2+ (coord_center_2(i,2)-coord_center(j,2))^2); %calculate all distances
             
                if dist<motion_range &&  i~=j && coord_center(i,4)~=coord_center(j,4) %&&  i~=j% if two points are close and not on same plane set RED coord_center(i,4)~=coord_center(j,4) 
                  
                        if coord_center(j,3)==0 && coord_center(j,6)>score_threshold %IF class Pseudomonas and threshold by score to eliminate low conf detections
                            dual_detect_class(i,1)=dual_detect_class(i,1)+1;                        
                        elseif coord_center(j,3)==1 &&  coord_center(j,6)>score_threshold %IF class Staph
                            dual_detect_class(i,2)=dual_detect_class(i,2)+1;                        

                        end
                        
                        %coord_center(j,:)=NaN; %remove point close by from coord_center
                               

                end
            
        end


         % plot(coord_center_2(i,1),size_pic-coord_center_2(i,2),'o','MarkerSize',8,'MarkerFaceColor',[1  0 0]); hold on; %plot detected cell in red
    
            cell_identity_local=dual_detect_class(:,1)./dual_detect_class(:,2);

            cell_identity=[cell_identity; cell_identity_local ];%ratio nb pseudomonas, nb staph
            cell_identity(cell_identity==Inf)=5; % Inf when 0 Staph, replace by 5

            nb_pseudo(t)=sum(cell_identity_local>=1); %total number Pseudomonas in this time point
            nb_staph(t)=sum(cell_identity_local<=1); %total number Staph in this time point

    end
    
    
end

    %text(0.8*size_pic,0.8*size_pic,num2str(nb_cells_total(t)),'FontSize',18);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%     next_frame=getframe(gcf);
%     writeVideo(v,next_frame);
    
% 
    max_class(t)=b; %slice where number of in--plane cells is max
     
        

end

% close(v); %close video

figure(2);
plot(nb_cells_total);
hold on;
plot(nb_pseudo);
hold on;
plot(nb_staph);

legend('Total','Pseudomonas','Staph');
FileInfoini = [time_folder '\foo0_0.jpeg'];
time_t0 = GetFileTime(FileInfoini);
time_tag0 = (time_t0.Write(4)*3600 + time_t0.Write(5)*60 + time_t0.Write(6)) + (time_t0.Write(3)-1)*86400; % time-point of the first image in seconds

for i = 1 : nb_time_points
    FileInfo = [time_folder strcat('\foo',int2str(i-1),'_0.jpeg')];
    time_t= GetFileTime(FileInfo);
    time_tag(i)= (time_t.Write(4)*3600 + time_t.Write(5)*60 + time_t.Write(6)) + (time_t.Write(3)-1)*86400 - time_tag0; % time-point of the i-th image in seconds
end

 plot(movmean(max_class,10));
figure(3);
plot(time_tag./60,nb_cells_total);
hold on;
plot(time_tag./60,nb_pseudo);
hold on;
plot(time_tag./60,nb_staph);

% figure(4);
% plot(movmean(nb_cells_total,10));
% hold on;
% plot(movmean(nb_pseudo,10));
% hold on;
% plot(movmean(nb_staph,10));

% Define the file name for the Excel file
output_excel_file = fullfile(output_folder, 'bacteria_analysis_results.xlsx');

% Prepare data for export by transposing the arrays
nb_pseudo_transposed = nb_pseudo';
nb_staph_transposed = nb_staph';
time_tag_transposed = time_tag';

% Combine data into a table for better readability in Excel
T = table(time_tag_transposed/60, nb_pseudo_transposed, nb_staph_transposed, ...
    'VariableNames', {'Time_Tag', 'Nb_Pseudo', 'Nb_Staph'});

% Write the table to the Excel file
writetable(T, output_excel_file);

% Save the time-tagged data
time_tagged_file = fullfile(output_folder, 'time_tagged_data.mat');
save(time_tagged_file, 'time_tag_transposed', 'nb_pseudo_transposed', 'nb_staph_transposed');

disp(['Analysis results saved to ', output_excel_file]);
disp(['Time-tagged data saved to ', time_tagged_file]);
%plot(time_tag./60,nb_cells_total)
save('nb_cells.mat'); %saves number of cell array for all time points

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


