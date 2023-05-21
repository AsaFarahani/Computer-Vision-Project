% Asa Borzabadi Farahani
clc
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laod the image
% You need to correct the address for the image 
img = imread('skyscrapers.jpg');

% Define thresholds
dr = 1
dt = 1 % degree 
Cth = 10;

% Convert degree to Radian
dt = dt*pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to grayscale
im = imresize(rgb2gray(img),0.5);

% Creates an array edges where each row is (x, y, cos theta, sin theta)   
Iedges = edge(im,'canny'); 
[~,grad_dir]=imgradient(im);
grad_dir = - grad_dir;

%  Now find all the edge locations, and add their orientations (cos theta,sin theta). 
%  row, col is  y,x
[row, col] = find(Iedges);
% Each edge is a 4-tuple:   (x, y, cos theta, sin theta)   
edges = [col, row, zeros(length(row),1), zeros(length(row),1) ];
for k = 1:length(row)
     edges(k,3) = cos(grad_dir(row(k),col(k))/180.0*pi);
     edges(k,4) = sin(grad_dir(row(k),col(k))/180.0*pi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RANSAC
% Apply randomness to the reported edge list 
edges = edges(randperm(size(edges, 1)), :);

% Number of reported edges
numbers = size(edges,1);

% Calculate r and theta(t) for all the edges / here the fifth column
% contains the r values and the sixth one contains the theta value
% r = xcos(t) + ysin(t)
edges(:,5) = edges(:,1).*edges(:,3) + edges(:,2).*edges(:,4);
% t = atan(sint/cost) 
for i = 1:1:size(edges,1)
    if ((edges(i,4)>0 & edges(i,3)>0) | (edges(i,4)<0 & edges(i,3)>0))
        edges(i,6) = atan(edges(i,4)./edges(i,3));
    end
    if ((edges(i,4)>0 & edges(i,3)<0))
        edges(i,6) = atan(edges(i,4)./edges(i,3))+pi;
    end
    if ((edges(i,4)<0 & edges(i,3)<0))
        edges(i,6) = +atan(edges(i,4)./edges(i,3))- pi;
    end
end

% Create a temp array containing the edges 
temp = edges;

counter_mainloop = size(edges,1);
counter_anything_left_or_not = size(edges,1);
Number_Significant_lines = 0;

for N = 1:1:counter_mainloop
    % Choose a random edge 
    selectedR_edge = temp(1,:);
    % Find r and theta(t) of the selected sample
    % Find theta of the chosen edge
    if ((selectedR_edge(4)>0 & selectedR_edge(3)>0) | (selectedR_edge(4)<0 & selectedR_edge(3)>0))
        tselected = atan(selectedR_edge(4)./selectedR_edge(3));
    end
    if ((selectedR_edge(4)>0 & selectedR_edge(3)<0))
        tselected = atan(selectedR_edge(4)./selectedR_edge(3)) + pi;
    end
    if ((selectedR_edge(4)<0 & selectedR_edge(3)<0))
        tselected = +atan(selectedR_edge(4)./selectedR_edge(3))- pi;
    end
    % Find r of the chosen edge
    rselected = selectedR_edge(1)*selectedR_edge(3) + selectedR_edge(2)*selectedR_edge(4);

    % Check the remaing edges to see weather they are inline with the selected edge 
    edges_inline(:,1) = ((abs(temp(:,5) - rselected) <= dr ) & (abs(temp(:,6) - tselected) <= dt));

    % Number of inline edges with the selected sample
    C(N,1) = sum(edges_inline);

    % Find the significant inline edges
    if (C(N,1)>= Cth)
        Number_Significant_lines = Number_Significant_lines +1;
        % Which edges are inline? identify them and save the mean information 
        % of all lines in a cluster in mean_values_edge_inline
        ci = 0;
        for i = 1:1:size(temp,1)
            if(edges_inline(i,1) == 1)
                ci = ci+1;
                values_edge_inline(ci,:) = temp(i,:);
                % Delete the selected sample and the samaples which are inline with it
                % from the obtained edges
                temp(i,:) = nan;
            end
        end
        % Save the significant lines
        significant_edges(Number_Significant_lines,:) = mean(values_edge_inline);
        counter_anything_left_or_not = counter_anything_left_or_not - ci;
    end
    % If C < Cmin, the line is nor powerful enough
    if (C(N,1)< Cth)
        ssize = size(temp,1);
        %temp(ssize+1,:) = temp(1,:);
        temp(1,:) = nan; % Discart the sample and continue
        counter_anything_left_or_not = counter_anything_left_or_not - 1;
    end
    value_on_break = N;
    % This is going to set the number of iterations for the main loop
    if counter_anything_left_or_not == 0
        N
        break;
    end 
    temp(isnan(temp(:,1)),:) = [];

    clear edges_inline
    clear values_edge_inline
end

% Create a blanck image and show the desired lines 
img_temp = zeros(size(im));
for n = 1:1:Number_Significant_lines
    for i = 1 :1:size(im,1)
        for j = 1 :1:size(im,2)
            if (abs(significant_edges(n,5) - j*(significant_edges(n,3)) - i*(significant_edges(n,4))) <=1)
                img_temp(i,j) = 1;
            end
        end
    end
end

% Display the results 
figure;set(gcf,'color','w');
subplot(1,4,1);imshow(im); title('Gray scaled image');
subplot(1,4,2);imshow(Iedges); title('canny edge detector');
subplot(1,4,3);imshow(img_temp,[]); title('RANSAC');
subplot(1,4,4);imshow(Iedges);
% Superimpose a contour map 
hold on;
contour((img_temp));title('RANSAC + detected edges');
hold off;
