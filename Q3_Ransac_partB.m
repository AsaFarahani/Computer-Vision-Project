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
dt = 1 %degree 
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

    % Check the remaing edges to see weather they are inline with the
    % selected edge 
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
        Sdata{Number_Significant_lines} = values_edge_inline;
        significant_edges(Number_Significant_lines,:) = mean(values_edge_inline);
        counter_anything_left_or_not = counter_anything_left_or_not - ci;
    end
    if (C(N,1)< Cth)
        ssize = size(temp,1);
        %temp(ssize+1,:) = temp(1,:);
        temp(1,:) = nan; % Discart the sample and continue
        counter_anything_left_or_not = counter_anything_left_or_not - 1;
    end
    %if (C(N,1)< Cth)
        %for i = 1:1:size(temp,1)
            %if(edges_inline(i,1) == 1)
              %  temp(i,:) = nan;
            %end
        %end
    %end
    value_on_break = N;
    if counter_anything_left_or_not == 0
        break;
    end 
    temp(isnan(temp(:,1)),:) = [];

    clear edges_inline
    clear values_edge_inline
end
figure;set(gcf,'color','w');
subplot(1,2,1);imshow(im)

% Total Least squares
for m = 1:1:size(Sdata,2)
    %{
        % Orthogonal linear regression method in 2D for model: y = a + bx   
        %
        % Input parameters:
        %  - XData: input data block -- x: axis
        %  - YData: input data block -- y: axis
        %  - vizualization: figure ('yes','no')
        % Return parameters:
        %  - Err: error - sum of orthogonal distances 
        %  - P: vector of model parameters [b-slope, a-offset
    %}
    YData = Sdata{m}(:,2);
    XData = Sdata{m}(:,1);
    kx=length(XData);
    ky=length(YData);
    n=size(YData,2);
    sy=sum(YData)./ky;
    sx=sum(XData)./kx;
    sxy=sum(XData.*YData);
    sy2=sum(YData.^2);
    sx2=sum(XData.^2);
    B=0.5.*(((sy2-ky.*sy.^2)-(sx2-kx.*sx.^2))./(ky.*sx.*sy-sxy));
    b1=-B+(B.^2+1).^0.5;
    b2=-B-(B.^2+1).^0.5;
    a1=sy-b1.*sx;
    a2=sy-b2.*sx;
    R=corrcoef(XData,YData);
    if R(1,2) > 0 
        P=[b1 a1];
        Yhat = XData.*b1 + a1;
        Xhat = ((YData-a1)./b1);
    end
    if R(1,2) < 0
        P=[b2 a2];
        Yhat = XData.*b2 + a2;
        Xhat = ((YData-a2)./b2);
    end   
    alpha = atan(abs((Yhat-YData)./(Xhat-XData)));
    d=abs(Xhat-XData).*sin(alpha);

    Err=sum(d.^2);
    hold on
    subplot(1,2,2);plot(XData,YData,'blue*'); 
    hold on;
    subplot(1,2,2);plot(XData,Yhat,'black');
    hold on
end
hold off

figure;set(gcf,'color','w');
imshow(im)

% Total Least squares
for m = 1:1:size(Sdata,2)
    %{
    % Orthogonal linear regression method in 2D for model: y = a + bx   
    %
    % Input parameters:
    %  - XData: input data block -- x: axis
    %  - YData: input data block -- y: axis
    %  - vizualization: figure ('yes','no')
    % Return parameters:
    %  - Err: error - sum of orthogonal distances 
    %  - P: vector of model parameters [b-slope, a-offset
    %}
    YData = Sdata{m}(:,2);
    XData = Sdata{m}(:,1);
    kx=length(XData);
    ky=length(YData);
    n=size(YData,2);
    sy=sum(YData)./ky;
    sx=sum(XData)./kx;
    sxy=sum(XData.*YData);
    sy2=sum(YData.^2);
    sx2=sum(XData.^2);
    B=0.5.*(((sy2-ky.*sy.^2)-(sx2-kx.*sx.^2))./(ky.*sx.*sy-sxy));
    b1=-B+(B.^2+1).^0.5;
    b2=-B-(B.^2+1).^0.5;
    a1=sy-b1.*sx;
    a2=sy-b2.*sx;
    R=corrcoef(XData,YData);
    if R(1,2) > 0 
        P=[b1 a1];
        Yhat = XData.*b1 + a1;
        Xhat = ((YData-a1)./b1);
    end
    if R(1,2) < 0
        P=[b2 a2];
        Yhat = XData.*b2 + a2;
        Xhat = ((YData-a2)./b2);
    end   
    alpha = atan(abs((Yhat-YData)./(Xhat-XData)));
    d=abs(Xhat-XData).*sin(alpha);

    Err=sum(d.^2);
    hold on
    plot(XData,YData,'blue*'); 
    hold on;
    plot(XData,Yhat,'black');
    hold on
end
hold off