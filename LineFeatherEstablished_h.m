function [Point1_On3Dline, Point2_On3Dline, L_3Dline] = LineFeatherEstablished_h(...
    Selected_corner_org, Selected_corner_left,...
        lineStack_h,centerStack_h,k_h,...
        CaliImg, pixelHeight, pixelWidth, radius, l_dis, pixelPitch)
%%====================================================================
%%===== calculating 3D line:筛选器 filtering
%%====================================================================
tic %%计时开始
d_img = 14.01; % Coarse distance between centers of adjacent macro images

filter1 = false;
filter2 = true;
nofilter = false;
%%%%%
if filter1 == true
    [lineStack_h_propotion_filter, centerStack_h_propotion_filter, num_PointsStack] = ...
        LineStackFilter1(Selected_corner_org, Selected_corner_left,...
        lineStack_h,centerStack_h,k_h,...
        l_dis, CaliImg, d_img);
end
%%%%%
if filter2 == true
    [lineStack_h_propotion_filter, centerStack_h_propotion_filter, num_PointsStack] = ...
        LineStackFilter2(Selected_corner_org, Selected_corner_left,...
        lineStack_h,centerStack_h,k_h,...
        CaliImg, radius);
end
%%%%%
if nofilter == true
    lineStack_h_propotion_filter = lineStack_h;
    centerStack_h_propotion_filter = centerStack_h;
    num_PointsStack = k_h;
end
%%====================================================================
%%===== calculating 3D line:initial solution:scheme1,2d line to 3d line
%%====================================================================
scheme1 = true;
if scheme1 == true
    figure;imshow(uint8(CaliImg));hold on;

%%% 使用lineStack_h中压栈的linefeather的参数，计算空间3D line
[Point1_On3Dline, Point2_On3Dline, L_3Dline] ...
    =Calculate3Dline(lineStack_h_propotion_filter, centerStack_h_propotion_filter, num_PointsStack,l_dis,pixelPitch);
%%% 根据两个空间点和 centerStack_h计算在每个宏像素中的重投影的line feather的方程
[Point_left, Point_right]= ExportLine3DLocalTerminals(...
    Point1_On3Dline, Point2_On3Dline, Selected_corner_org, Selected_corner_left, d_img, true);
[lineStack_h_filtered, k_h_filtered, lineStack_reproject]=...
    reprojection_linefeather(Point_left, Point_right, l_dis, pixelPitch, ...
    centerStack_h, lineStack_h, k_h, radius);
%%%%%% 画图
%Plot_3DLineXYProjection(Point1_On3Dline, Point2_On3Dline, Selected_corner_org, Selected_corner_left, d_img);
color = 'y-';Plot_LineFeather_h(radius, lineStack_reproject, k_h, centerStack_h,color);
hold off;
disp('Wait');close all;
%%%%%%
Point_left = [2313;2260;-95];
Point_right = [3332;2187;-115];
%{
   %figure;imshow(uint8(CaliImg));hold on;
%%% 使用lineStack_h中压栈的linefeather的参数，计算空间3D line
[Point1_On3Dline, Point2_On3Dline, L_3Dline] ...
    =Calculate3Dline(lineStack_reproject, centerStack_h, k_h,l_dis,pixelPitch);
%%% 根据两个空间点和 centerStack_h计算在每个宏像素中的重投影的line feather的方程
%}
%{ 
%Test通过手动调节Point_left和Point_right，观察重投影linefeather的图像，并总结出
%在下一步迭代时，对这两个点的调节不一定全遍历。有可能可以像凸优化似的，一步一步的找下降方向和步长。
%figure;imshow(uint8(CaliImg));hold on;
%color = 'y-';Plot_LineFeather_h(radius, lineStack_h, k_h, centerStack_h,color);hold off;
figure;imshow(uint8(CaliImg));hold on;
[~, ~]= ExportLine3DLocalTerminals(...
    Point_left, Point_right, Selected_corner_org, Selected_corner_left, d_img, true);

[lineStack_h_filtered, k_h_filtered, lineStack_reproject]=...
    reprojection_linefeather(Point_left, Point_right, l_dis, pixelPitch, ...
    centerStack_h, lineStack_h, k_h, radius);
color = 'g-';Plot_LineFeather_h(radius, lineStack_reproject, k_h, centerStack_h,color);
hold off;
disp('Wait');close all;
%}
end
%%====================================================================
%%===== calculating 3D line:initial solution:scheme2,2d point to 3d line
%%====================================================================
%%%%% one point in every micro image
%%%%%方法1，所有光线的交点的点集的svd分解的前两个基向量
%%%%%方法2，子孔径图像中的线作为初始，在raw中每个宏像素中只取一个点，根据子孔径图像的直线方程只拟合z的坐标，然后再
%%%%%综合x的坐标，共同拟合3d的直线

%%====================================================================
%%============= calculating 3D line:iterated refine
%%====================================================================
%%% generate很多的系列“一致性”模板
%%%由初始解，在给定适合的沿z轴和y轴，的步长的情况下，给出一系列“一致性”模板
[Point1_On3DlineOffset, Point2_On3DlineOffset, L_3Dline, MaxNCC]=...
    IteratedRefine3DLine(Point_left, Point_right, l_dis, pixelPitch, ...
    centerStack_h, lineStack_h, k_h, radius, pixelHeight, pixelWidth, CaliImg);
toc %% 计时结束
end
%%%%%%%%%%%%%% Several Sub functions %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%