%% 1 importing videos from a folder of videos named sequentially

files = dir('/Users/hilina_baby/Documents/MATLAB/*.m4v');
file = files';

delete video_frames2/*.jpg

video = VideoReader('/Users/hilina_baby/Documents/MATLAB/trial_2.m4v'); % reads the video,give the file name

%% 2 count # frames
first_NumFr = video.NumberOfFrames; % Tells you the number of frames that have been taken from the video 
fps = 30; % Frames per second: velocity of the camara 
%Fs = 1200; % Sampling frequency: (datapoints or pixels) / (width); 

%% 3 Find all frames in video and write them to a folder video_frames2/frame#.jpg
% Extract frames from video

for img = 1: first_NumFr
    FrameN = strcat('frame', num2str(img),'.jpg'); % Order the frames in a horizontal vector and give them the file format jpg
           n = read(video,img); %get the values of each frame NumFr
    imwrite(n,['video_frames2/frame',sprintf('%03d',img),'.jpg']);
end
% deleting frames in between
for i = 1:first_NumFr
    file_to_delete = ['video_frames2/' 'frame', sprintf('%03d',i),'.jpg'];
    if(mod(i,10) == 0)
         disp('Nothing');    
    else
        delete(file_to_delete);   
    end
end
fls = dir('video_frames2');
NumFr = length(fls)-2;
A_temp=dir('./video_frames2');
% renaming frame numbers
for i=1:NumFr
    v1= ['video_frames2/',A_temp(i+2).name];
    v2= ['video_frames2/frame',num2str(i),'.jpg']
    
    movefile(v1,v2);
end
%% 4 extract the protractor lines in vector form and measure the angles
%find the lines

list_excel=zeros(NumFr,1);

for im = 1:NumFr   % For loop which processes all the frames, one by one
  %  subject_frame = imread(['video_frames2/frame', num2str(im),'.jpg']);% CREATE this file to store frames first. reads the frame
     %subject_frame = imread(['video_frames2/frame', num2str(im),'.jpg']);
     %cropped_frame = imcrop( subject_frame,[10 50 800 700]); 
     gray_frame = rgb2gray(subject_frame);  % Turns sup2 into image into a binary image and define the thershold value to minimize the interclass variance of the black and white 
imshow(gray_frame)
     BW_frame= edge(gray_frame,'Prewitt',0.03,'nothinning');
     [H, theta,rho]= hough(BW_frame,'ThetaResolution',0.5); %add specifications to the hough function
     peaks=houghpeaks(H,7);
     lines = houghlines(BW_frame,theta,rho,peaks,'FillGap',5,'MinLength',7);
     
     %figure, imshow(BW_frame), 
     hold on
     max_len = 0;
     length_store=zeros(length(lines),1);
     length_field= zeros(length(lines),2);
             for k = 1:length(lines)
             xy = [lines(k).point1; lines(k).point2];
             plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    
             % Plot beginnings and ends of lines
             plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
             plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

             % Determine the endpoints of the longest line segment
             len = norm(lines(k).point1 - lines(k).point2);
                  if ( len > max_len)
                     max_len = len;
                     xy_long = xy;
                  end
             length_store(k) = norm(lines(k).point1 - lines(k).point2);
             length_field(k,:)= [length_store(k),k];
             end
             
             length_field;
             sorted_field= sortrows(length_field);
             
             
            long_line_field_1=sorted_field(k,2); %longest line 
            long_line_field_2=sorted_field(k-1,2); %second longest line
            
            x=lines(long_line_field_1).point1; 
            y=lines(long_line_field_1).point2;
            
            a=lines(long_line_field_2).point1;
            b=lines(long_line_field_2).point2;
            
            %gamma_between = (atan((y(2)-y(1))/(x(2)-x(1))) - atan((b(2)-b(1))/(a(2)-a(1)))) * 180/pi; %measure angle between
    
u1=y(1)-x(1);
u2=y(2)-x(2);
v1=b(1)-a(1);
v2=b(2)-a(2);

u=[u1 u2 0];
v=[v1 v2 0];

ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v))  
  
    for lil_loop=2:50;
                if ThetaInDegrees< 10;
                long_line_field_2=sorted_field(k-lil_loop,2); %use the next longest line
                a=lines(long_line_field_2).point1;
                b=lines(long_line_field_2).point2;
                v1=b(1)-a(1); v2=b(2)-a(2); v=[v1 v2 0];
               
                ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v)) %NEW, CORRECT THETA
                end
    end

%% 5 store values from a loop

list_excel(im)= ThetaInDegrees;
end
list_excel
figure;  
t_boundary= (1/fps)*NumFr*10;
t_fr= 0: 1/fps*10: (NumFr-1)/fps*10;
y_angle= [list_excel(:,1)'];
plot(t_fr, y_angle, 'g*')
title('time vs angle plot');
xlabel('time [s]');
ylabel('arm angle [theta]');

%output time
time_output=[t_fr'];

%export
arm_angle_file = 'armangledata.xlsx';
%xlswrite(arm_angle_file,list_excel,'Arm Angles')

%subplots!!




