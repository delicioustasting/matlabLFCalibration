function [Point1_On3Dline, Point2_On3Dline, L_3Dline] = LineFeatherEstablished_h(...
    Selected_corner_org, Selected_corner_left,...
        lineStack_h,centerStack_h,k_h,...
        CaliImg, pixelHeight, pixelWidth, radius, l_dis, pixelPitch)
%%====================================================================
%%===== calculating 3D line:ɸѡ�� filtering
%%====================================================================
tic %%��ʱ��ʼ
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

%%% ʹ��lineStack_h��ѹջ��linefeather�Ĳ���������ռ�3D line
[Point1_On3Dline, Point2_On3Dline, L_3Dline] ...
    =Calculate3Dline(lineStack_h_propotion_filter, centerStack_h_propotion_filter, num_PointsStack,l_dis,pixelPitch);
%%% ���������ռ��� centerStack_h������ÿ���������е���ͶӰ��line feather�ķ���
[Point_left, Point_right]= ExportLine3DLocalTerminals(...
    Point1_On3Dline, Point2_On3Dline, Selected_corner_org, Selected_corner_left, d_img, true);
[lineStack_h_filtered, k_h_filtered, lineStack_reproject]=...
    reprojection_linefeather(Point_left, Point_right, l_dis, pixelPitch, ...
    centerStack_h, lineStack_h, k_h, radius);
%%%%%% ��ͼ
%Plot_3DLineXYProjection(Point1_On3Dline, Point2_On3Dline, Selected_corner_org, Selected_corner_left, d_img);
color = 'y-';Plot_LineFeather_h(radius, lineStack_reproject, k_h, centerStack_h,color);
hold off;
disp('Wait');close all;
%%%%%%
Point_left = [2313;2260;-95];
Point_right = [3332;2187;-115];
%{
   %figure;imshow(uint8(CaliImg));hold on;
%%% ʹ��lineStack_h��ѹջ��linefeather�Ĳ���������ռ�3D line
[Point1_On3Dline, Point2_On3Dline, L_3Dline] ...
    =Calculate3Dline(lineStack_reproject, centerStack_h, k_h,l_dis,pixelPitch);
%%% ���������ռ��� centerStack_h������ÿ���������е���ͶӰ��line feather�ķ���
%}
%{ 
%Testͨ���ֶ�����Point_left��Point_right���۲���ͶӰlinefeather��ͼ�񣬲��ܽ��
%����һ������ʱ������������ĵ��ڲ�һ��ȫ�������п��ܿ�����͹�Ż��Ƶģ�һ��һ�������½�����Ͳ�����
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
%%%%%����1�����й��ߵĽ���ĵ㼯��svd�ֽ��ǰ����������
%%%%%����2���ӿ׾�ͼ���е�����Ϊ��ʼ����raw��ÿ����������ֻȡһ���㣬�����ӿ׾�ͼ���ֱ�߷���ֻ���z�����꣬Ȼ����
%%%%%�ۺ�x�����꣬��ͬ���3d��ֱ��

%%====================================================================
%%============= calculating 3D line:iterated refine
%%====================================================================
%%% generate�ܶ��ϵ�С�һ���ԡ�ģ��
%%%�ɳ�ʼ�⣬�ڸ����ʺϵ���z���y�ᣬ�Ĳ���������£�����һϵ�С�һ���ԡ�ģ��
[Point1_On3DlineOffset, Point2_On3DlineOffset, L_3Dline, MaxNCC]=...
    IteratedRefine3DLine(Point_left, Point_right, l_dis, pixelPitch, ...
    centerStack_h, lineStack_h, k_h, radius, pixelHeight, pixelWidth, CaliImg);
toc %% ��ʱ����
end
%%%%%%%%%%%%%% Several Sub functions %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%