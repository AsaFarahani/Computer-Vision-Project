% Asa Borzabadi Farahani
clc
close all
clear all
%% PART A
% You need to correct the address for the image 
img = imread('skyscrapers.jpg');

% You can change the value of sigma by chanaging the caue of this variable.
sigma = 4;

% Dimentions of the array
size_array = 6*sigma+1;

% Generate a 1D Gaussian filter 
Gaussian = fspecial('gaussian',[size_array+4 1],sigma);

% Plot the 1D Gaussian 
figure;set(gcf,'color','w');
m =((size(Gaussian,1)-1)/2) % Centered around 0
subplot(1,2,1);plot(-m:m,Gaussian,'-or');
title(['1D gaussian - simga = ' num2str(sigma)]);

% Generate the first derivative of the Gaussian
D_Gaussian = conv(Gaussian, [-1/2 0  1/2], 'valid');

% Generate the second derivative of the Gaussian
D2_Gaussian = conv(D_Gaussian, [-1/2 0  1/2], 'valid');

% Plot the second derivative of the 1D Gaussian
m = ((size(D2_Gaussian,1)-1)/2);   % Centered around 0
subplot(1,2,2);plot(-m:m,D2_Gaussian, '-ob');
title(['Second derivative 1D gaussian - simga = ' num2str(sigma)]);

% Create the 2D filter array (size = (2sigma + 1)* 2sigma +1))
Filter_A = repmat(D2_Gaussian,length(size_array),size_array);

% Plot the 2D array filter 
figure;set(gcf,'color','w');
m = (size(Filter_A,1)-1)/2; % centered around 0
subplot(1,2,1);mesh(-m:m,-m:m,Filter_A);ylabel('x');xlabel('y');colorbar;
subplot(1,2,2);imagesc(-m:m,-m:m,Filter_A);ylabel('x');xlabel('y');colorbar;
%% PART B
% Modifing the previous filter
sigma_y = sigma*2;
Gaussian_y = fspecial('gaussian',[1 size_array],sigma_y);
for i = 1:1: size_array
    Filter_B(i,:) = Filter_A(i,:).*Gaussian_y;
end

% Plot the Modified filter 
m = (size(Filter_B,1)-1)/2; % centered around 0 
figure;set(gcf,'color','w');
subplot(1,2,1);mesh(-m:m,-m:m,Filter_B);ylabel('x');xlabel('y');;colorbar;
subplot(1,2,2);imagesc(-m:m,-m:m,Filter_B);ylabel('x');xlabel('y');colorbar;
%% PART C
% Plotting the rotated filters
theta = 0;
[m, out0] = rotation_filter(Filter_B, theta);
figure;set(gcf,'color','w');
subplot(4,2,1);mesh(-m:m,-m:m,out0);colorbar;title(['theta ' num2str((theta)*180/pi)]);
subplot(4,2,2);imagesc(-m:m,-m:m,out0);colorbar;title(['theta ' num2str((theta)*180/pi)]);
theta = pi/4;
[m, out45] = rotation_filter(Filter_B, theta);
subplot(4,2,3);mesh(-m:m,-m:m,out45);colorbar;title(['theta ' num2str((theta)*180/pi)]);
subplot(4,2,4);imagesc(-m:m,-m:m,out45);colorbar;title(['theta ' num2str((theta)*180/pi)]);
theta = pi/2;
[m, out90] = rotation_filter(Filter_B, theta);
subplot(4,2,5);mesh(-m:m,-m:m,out90);colorbar;title(['theta ' num2str((theta)*180/pi)]);
subplot(4,2,6);imagesc(-m:m,-m:m,out90);colorbar;title(['theta ' num2str((theta)*180/pi)]);
theta = 3*pi/4;
[m, out135] = rotation_filter(Filter_B, theta);
subplot(4,2,7);mesh(-m:m,-m:m,out135);colorbar;title(['theta ' num2str((theta)*180/pi)]);
subplot(4,2,8);imagesc(-m:m,-m:m,out135);colorbar;title(['theta ' num2str((theta)*180/pi)]);
%% PART D
% Laod the image

figure;set(gcf,'color','w');
% Plot the colored image 
subplot(1,2,1);imshow(img); title('Original picture');
% Conver the image to gray scale 
imgG = rgb2gray(img);
% Plot the gray scaled image 
subplot(1,2,2);imshow(imgG); title('Black and white picture');
% Apply the retated filters to the image that we have
figure;set(gcf,'color','w');
F0_imgG = conv2((imgG),out0,'same');
subplot(2,2,1);imshow(F0_imgG,[]);title('Filtered image - 0');
F45_imgG = conv2((imgG),out45,'same');
subplot(2,2,2);imshow(F45_imgG,[]);title('Filtered image - 45');
F90_imgG = conv2((imgG),out90,'same');
subplot(2,2,3);imshow(F90_imgG,[]);title('Filtered image - 90');
F135_imgG = conv2((imgG),out135,'same');
subplot(2,2,4);imshow(F135_imgG,[]);title('Filtered image - 135');

% Find the Zero_rossings (separately for each filtering version)
% Here, we used shifting to find the zero scorrings 
figure;set(gcf,'color','w');
temp0 = zeros(size(F0_imgG));
temp0(:,:,1) = sign(circshift(F0_imgG,-1,1)) ~= sign(F0_imgG);
temp0(:,:,2) = sign(circshift(F0_imgG,-1,2)) ~= sign(F0_imgG);
t0 = temp0(:,:,1) + temp0(:,:,2);
subplot(2,2,1);imshow(t0); title('Zero-crossings obtained by a 0-degree oriented filter');
temp45 = zeros(size(F45_imgG));
temp45(:,:,1) = sign(circshift(F45_imgG,-1,1)) ~= sign(F45_imgG);
temp45(:,:,2) = sign(circshift(F45_imgG,-1,2)) ~= sign(F45_imgG);
t45 = temp45(:,:,1) + temp45(:,:,2);
subplot(2,2,2);imshow(t45);title('Zero-crossings obtained by a 45-degree oriented filter');
temp90 = zeros(size(F90_imgG));
temp90(:,:,1) = sign(circshift(F90_imgG,-1,1)) ~= sign(F90_imgG);
temp90(:,:,2) = sign(circshift(F90_imgG,-1,2)) ~= sign(F90_imgG);
t90 = temp90(:,:,1) + temp90(:,:,2);
subplot(2,2,3);imshow(t90);title('Zero-crossings obtained by a 90-degree oriented filter');
temp135 = zeros(size(F135_imgG));
temp135(:,:,1) = sign(circshift(F135_imgG,-1,1)) ~= sign(F135_imgG);
temp135(:,:,2) = sign(circshift(F135_imgG,-1,2)) ~= sign(F135_imgG);
t135 = temp135(:,:,1) + temp135(:,:,2);
subplot(2,2,4);imshow(t135);title('Zero-crossings obtained by a 135-degree oriented filter');

% Merge all the Zero-crossings obained from different filters
tempall = t0+t45+t90+t135;
figure;set(gcf,'color','w');
imshow(tempall);title('Zero-crossings');

%% PART E - Laplacian of a Gaussian filter with ? = 6 and width 40
% h = fspecial('log',hsize,sigma) returns a rotationally symmetric Laplacian 
% of Gaussian filter of size hsize with standard deviation sigma.

% Define the laplacian filter
hsize = 40;
sigma = 6;
filter_e = fspecial('log',hsize,sigma);
% Plot the laplacian filter
figure;set(gcf,'color','w');
subplot(2,1,2),surf(filter_e); 
subplot(2,1,1),imagesc(filter_e);title(' The lacplacian filter used in part E');

% Find the filtered image
Filtered_img_partE = conv2((imgG),filter_e,'same');
% Plot the filtered image
figure;set(gcf,'color','w');
subplot(2,1,1),imshow(Filtered_img_partE); title(' The filtered image - Part E');

% Find zero crossing for the filtered image
temp_e = zeros(size(Filtered_img_partE));
temp_e(:,:,1) = sign(circshift(Filtered_img_partE,-1,1)) ~= sign(Filtered_img_partE);
temp_e(:,:,2) = sign(circshift(Filtered_img_partE,-1,2)) ~= sign(Filtered_img_partE);
zerocrossings_laplacian = temp_e(:,:,1) + temp_e(:,:,2);
subplot(2,1,2);imshow(zerocrossings_laplacian); title('Zero-crossings obtained by using the Laplacian filter');
%% PART C - Rotation Function
function [m, out] = rotation_filter(Fy_D2G, theta)
    m = (size(Fy_D2G,1) - 1)/2;
    [x,y] = meshgrid(-m:m,-m:m);     % non-rotated coordinate system, contains (0,0)
    u = cos(theta)*x - sin(theta)*y; % rotated coordinate system
    v = sin(theta)*x + cos(theta)*y; % rotated coordinate system
    udecimal = int64(round(u)+m+1);
    vdecimal = int64(round(v)+m+1);
    out = nan(size(Fy_D2G));
    for i = 1:1:size(Fy_D2G,1)
        for j = 1:1:size(Fy_D2G,2)
            if (udecimal(i,j) >= 1 && vdecimal(i,j) >= 1 && vdecimal(i,j) <= size(Fy_D2G,2) && udecimal(i,j) <= size(Fy_D2G,1))
                out(i,j) = Fy_D2G(udecimal(i,j),vdecimal(i,j));
            end

        end
    end
    % A more complexed method - interpolation of the NAN values    
    for i =1:1:size(Fy_D2G,2)
       for j = 1:1:size(Fy_D2G,2)
           if (udecimal(i,j) < 1)
                udecimal(i,j)  = 1;
           end
           if (vdecimal(i,j) < 1)
                vdecimal(i,j)  = 1;
           end
           if (udecimal(i,j) > size(Fy_D2G,2))
                udecimal(i,j)  = size(Fy_D2G,2);
           end
           if (vdecimal(i,j) > size(Fy_D2G,2))
                vdecimal(i,j) = size(Fy_D2G,2);
           end
       end
    end
    for i = 1:1:size(Fy_D2G,1)
        for j = 1:1:size(Fy_D2G,2)
            if (udecimal(i,j) >= 1 && vdecimal(i,j) >= 1 && vdecimal(i,j) <= size(Fy_D2G,2) && udecimal(i,j) <= size(Fy_D2G,1))
                out(i,j) = Fy_D2G(udecimal(i,j),vdecimal(i,j));
            end

        end
    end
end