function [lineStack_h_propotion_filter, centerStack_h_propotion_filter, num_PointsStack] = ...
    LineStackFilter2(Selected_corner_org, Selected_corner_left,...
    lineStack_h,centerStack_h,k_h,...
    CaliImg, radius)
%%% Paras of line which links left point and right point
Line_SubImg = zeros(3,1);
Line_SubImg_k = (Selected_corner_left(2,1) - Selected_corner_org(2,1)) /...
    (Selected_corner_left(1,1) - Selected_corner_org(1,1));
Line_SubImg_b = -Line_SubImg_k*Selected_corner_org(1,1)+Selected_corner_org(2,1);
%{
    figure;imshow(CenterSubImg);hold on;
    x_temp = round(Selected_corner_left(1,1)):1:round(Selected_corner_org(1,1));
    y_temp = Line_SubImg_k*x_temp + Line_SubImg_b;
    plot(x_temp,y_temp,'g-','LineWidth',1);    
%}
%%selected 2D point
%%====================== important variables =====
lineStack_h_propotion_filter = zeros(3,k_h);
centerStack_h_propotion_filter = zeros(2,k_h);
num_PointsStack = 0;
%%==================================
figure;imshow(CaliImg);hold on;
for i = 1:k_h % 当前这里面的变量还是以像素为单位的吧
     %%%筛选，原始linefeather的斜率如果比子孔径图像中的初始linefeather的斜率Line_SubImg_k，夹角大于5度，就不要了。
    if abs(lineStack_h(1:2,i)'*[Line_SubImg_k; -1])/(norm(lineStack_h(1:2,i))*norm([Line_SubImg_k; -1]))...
            >= cos(5*pi/180) % 是8度，5度，还是几度，这里是一个手调参数
            % 两条线的方向向量的夹角(即内积除以模值，即夹角的余弦)，也即两条线的法向量的夹角
        num_PointsStack = num_PointsStack + 1;
        lineStack_h_propotion_filter(:,num_PointsStack) = lineStack_h(:, i);
        centerStack_h_propotion_filter(:,num_PointsStack) = centerStack_h(:, i);
    end
end
color = 'g-';
Plot_LineFeather_h(radius, lineStack_h_propotion_filter, num_PointsStack, centerStack_h_propotion_filter,color);
hold off;
disp('Wait');close all;
end