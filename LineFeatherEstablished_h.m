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
%%%����Point1_On3Dline
Z_org1 = Point_left(3,1);
Y_org1 = Point_left(2,1);
[Point1Pixel_Z_stack, Point1PixelY_stack] = GenerateOffsetZX(Z_org1, Y_org1, l_dis, radius, pixelPitch);
%%%����Point2_On3Dline
Z_org2 = Point_right(3,1);
Y_org2 = Point_right(2,1);
[Point2Pixel_Z_stack, Point2PixelY_stack] = GenerateOffsetZX(Z_org2, Y_org2, l_dis, radius, pixelPitch);
clear Z_org1 X_org1 Z_org2 X_org2;
%figure;
MaxNCC = 0;MaxNCC_idx = [0, 0, 0, 0];
NCC = zeros(size(Point1Pixel_Z_stack,2),size(Point1PixelY_stack,2),...
    size(Point2Pixel_Z_stack,2),size(Point2Pixel_Z_stack,2));
for Z1_idx = 1:size(Point1Pixel_Z_stack,2)
    for X1_idx = 1:size(Point1PixelY_stack,2)
        for Z2_idx = 1:size(Point2Pixel_Z_stack,2)
            for X2_idx = 1:size(Point2PixelY_stack,2)
                Point1_On3DlineOffset = Point_left;
                Point2_On3DlineOffset = Point_right;
                %%%%% �滻
                Point1_On3DlineOffset(3,1) = Point1Pixel_Z_stack(1,Z1_idx);
                Point1_On3DlineOffset(2,1) = Point1PixelY_stack(1,X1_idx);
                Point2_On3DlineOffset(3,1) = Point2Pixel_Z_stack(1,Z2_idx);
                Point2_On3DlineOffset(2,1) = Point2PixelY_stack(1,X2_idx);
                %%%%%
                [lineStack_h_filtered, k_h_filtered, lineStack_reproject]=...
                    reprojection_linefeather(Point1_On3DlineOffset, Point2_On3DlineOffset, ...
                    l_dis, pixelPitch, centerStack_h, lineStack_h, k_h, radius);
                %%%%% ��ͼ
                %{
                f_temp = figure;
                imshow(uint8(CaliImg));hold on;
                b=[-radius+0.5,radius-0.5];
                for i = 1:k_h
                    boundary=[b;-(b*lineStack_reproject(1,i)+lineStack_reproject(3,i))/lineStack_reproject(2,i)];
                    plot(boundary(1,:)+centerStack_h(1,i),boundary(2,:)+centerStack_h(2,i),'g-','LineWidth',1);
                end
                hold off;close;
                %}
                %%%%%
                TempNCC = TotalNCC(centerStack_h, lineStack_reproject,...
                    k_h, radius, pixelHeight, pixelWidth, CaliImg);
                NCC(Z1_idx, X1_idx,Z2_idx,X2_idx) = TempNCC;
                if TempNCC >MaxNCC
                    MaxNCC = TempNCC;
                    MaxNCC_idx = [Z1_idx, X1_idx,Z2_idx,X2_idx];
                end
                %{
                RawImgTemplate = GenarateJointTemplate(centerStack_h, lineStack_reproject, k_h, radius,...
                    pixelHeight, pixelWidth);
                    %}
                %imshow(RawImgTemplate);hold on;
            end
        end
    end
end
% ע�����Ļ�����������Ϊ��λ��������⵽���Ⱥ�����ٳ���pixelPitch.
%hold off;
disp('Wait');
%{
maximum_ncc = max(max(max(max(NCC))));
for Z1_idx = 1:size(Point1Pixel_Z_stack,2)
    for X1_idx = 1:size(Point1PixelY_stack,2)
        for Z2_idx = 1:size(Point2Pixel_Z_stack,2)
            for X2_idx = 1:size(Point2PixelY_stack,2)
                if NCC(Z1_idx, X1_idx,Z2_idx,X2_idx) == maximum_ncc
                    maximum_idx = [Z1_idx, X1_idx,Z2_idx,X2_idx];
                end
            end
        end
    end
end
%}
Point1_On3DlineOffset = Point_left;
Point2_On3DlineOffset = Point_right;
%%%%% �滻
Point1_On3DlineOffset(3,1) = Point1Pixel_Z_stack(1,MaxNCC_idx(1,1));
Point1_On3DlineOffset(2,1) = Point1PixelY_stack(1,MaxNCC_idx(1,2));
Point2_On3DlineOffset(3,1) = Point2Pixel_Z_stack(1,MaxNCC_idx(1,3));
Point2_On3DlineOffset(2,1) = Point2PixelY_stack(1,MaxNCC_idx(1,4));
%%%%%
[lineStack_h_filtered, k_h_filtered, lineStack_reproject]=...
    reprojection_linefeather(Point1_On3DlineOffset, Point2_On3DlineOffset, ...
    l_dis, pixelPitch, centerStack_h, lineStack_h, k_h, radius);
[Point1_On3Dline, Point2_On3Dline, L_3Dline] ...
    =Calculate3Dline(lineStack_reproject, centerStack_h, k_h,l_dis,pixelPitch);
toc %% ��ʱ����
%%%%% ��ͼ
f_temp = figure;
imshow(uint8(CaliImg));hold on;
color = 'y-';Plot_LineFeather_h(radius, lineStack_reproject, k_h, centerStack_h,color);
[~, ~]= ExportLine3DLocalTerminals(...
    Point1_On3DlineOffset, Point2_On3DlineOffset, Selected_corner_org, Selected_corner_left, d_img, true);
hold off;close;
end
%%%%%%%%%%%%%% Several Sub functions %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z_stack, X_stack] = GenerateOffsetZX(Z_org, X_org, l_dis, radius, pixelPitch)
Z_org = Z_org*pixelPitch; % ������Ϊ��λת���Գ���Ϊ��λ
X_org = X_org*pixelPitch; % ������Ϊ��λת���Գ���Ϊ��λ

TanTheta = radius*pixelPitch/l_dis;
D_half = Z_org*TanTheta;
maxOffset_org = D_half*((Z_org+l_dis)/Z_org);
%%% ǰ����for��Point1_3D�ı仯��������for��Point2_3D�ı仯
maxOffsetVector = (maxOffset_org - 4*pixelPitch):(0.25*pixelPitch):(maxOffset_org + 4*pixelPitch);
Z_stack = zeros(size(maxOffsetVector));
for i = 1:size(maxOffsetVector,2)
    Z_modified = (l_dis*D_half)/(maxOffsetVector(1,i) - D_half);
    Z_modified_pixel = Z_modified/pixelPitch; % ע�����Ļ�����������Ϊ��λ��������⵽���Ⱥ�����ٳ���pixelPitch.
    
    if Z_modified < 100*l_dis % ��ֹƫ��̫��ʹ�����ɵ�Z_modified�������������޴�
        Z_stack(1,i) = Z_modified_pixel;
    end
end
X_offsetVector = linspace((-4*pixelPitch*abs(Z_modified/l_dis)),(4*pixelPitch*abs(Z_modified/l_dis)),33);
X_stack = zeros(size(X_offsetVector));
for j = 1:size(X_offsetVector,2)
    X_modified = X_org + X_offsetVector(1,j);
    X_modified_pixel = X_modified/pixelPitch; % ע�����Ļ�����������Ϊ��λ��������⵽���Ⱥ�����ٳ���pixelPitch.
    
    X_stack(1,j) = X_modified_pixel;
end
% ��ȷ����ϵ���С�仯�������Ӳ��������Ӳ�ķ�Χ���ֵ�
end