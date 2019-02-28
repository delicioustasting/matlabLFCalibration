function [Point1_On3DlineOffset, Point2_On3DlineOffset, L_3Dline, MaxNCC]=...
    IteratedRefine3DLine(Point_left, Point_right, l_dis, pixelPitch, ...
    centerStack_h, lineStack_h, k_h, radius, pixelHeight, pixelWidth, CaliImg)
%%%对于Point1_On3Dline
Z_org1 = Point_left(3,1);
Y_org1 = Point_left(2,1);
[Point1Pixel_Z_stack, Point1PixelY_stack] = GenerateOffsetZX(Z_org1, Y_org1, l_dis, radius, pixelPitch);
%%%对于Point2_On3Dline
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
                %%%%% 替换
                Point1_On3DlineOffset(3,1) = Point1Pixel_Z_stack(1,Z1_idx);
                Point1_On3DlineOffset(2,1) = Point1PixelY_stack(1,X1_idx);
                Point2_On3DlineOffset(3,1) = Point2Pixel_Z_stack(1,Z2_idx);
                Point2_On3DlineOffset(2,1) = Point2PixelY_stack(1,X2_idx);
                %%%%%
                [lineStack_h_filtered, k_h_filtered, lineStack_reproject]=...
                    reprojection_linefeather(Point1_On3DlineOffset, Point2_On3DlineOffset, ...
                    l_dis, pixelPitch, centerStack_h, lineStack_h, k_h, radius);
                %%%%% 画图
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
% 注意最后的还得是以像素为单位，于是求解到长度后最后再除以pixelPitch.
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
%%%%% 替换
Point1_On3DlineOffset(3,1) = Point1Pixel_Z_stack(1,MaxNCC_idx(1,1));
Point1_On3DlineOffset(2,1) = Point1PixelY_stack(1,MaxNCC_idx(1,2));
Point2_On3DlineOffset(3,1) = Point2Pixel_Z_stack(1,MaxNCC_idx(1,3));
Point2_On3DlineOffset(2,1) = Point2PixelY_stack(1,MaxNCC_idx(1,4));
%%%%%
[~,~,V_space] = svd([[Point1_On3DlineOffset',1]; [Point2_On3DlineOffset',1]]);
L_3Dline = V_space(3:4,:); % 输出最后的L_3Dline modified by jdy 20190227
%%%%%
%{
figure;plot3([Point1_On3DlineOffset(1,1), Point2_On3DlineOffset(1,1)],...
        [Point1_On3DlineOffset(2,1), Point2_On3DlineOffset(2,1)],...
        [Point1_On3DlineOffset(3,1), Point2_On3DlineOffset(3,1)]);grid on;
%}
[lineStack_h_filtered, k_h_filtered, lineStack_reproject]=...
    reprojection_linefeather(Point1_On3DlineOffset, Point2_On3DlineOffset, ...
    l_dis, pixelPitch, centerStack_h, lineStack_h, k_h, radius);
%%%%% 画图
f_temp = figure;
imshow(uint8(CaliImg));hold on;
color = 'y-';Plot_LineFeather_h(radius, lineStack_reproject, k_h, centerStack_h,color);
[~, ~]= ExportLine3DLocalTerminals(...
    Point1_On3DlineOffset, Point2_On3DlineOffset, Selected_corner_org, Selected_corner_left, d_img, true);
hold off;close;
%{
[Point1_On3Dline, Point2_On3Dline, L_3Dline] ...
    =Calculate3Dline(lineStack_reproject, centerStack_h, k_h,l_dis,pixelPitch);
%}
end

function [Z_stack, X_stack] = GenerateOffsetZX(Z_org, X_org, l_dis, radius, pixelPitch)
Z_org = Z_org*pixelPitch; % 以像素为单位转成以长度为单位
X_org = X_org*pixelPitch; % 以像素为单位转成以长度为单位

TanTheta = radius*pixelPitch/l_dis;
D_half = abs(Z_org*TanTheta); % modified by JDY 20190227
maxOffset_org = D_half*((Z_org+l_dis)/Z_org);
%%% 前两个for是Point1_3D的变化，后两个for是Point2_3D的变化
maxOffsetVector = (maxOffset_org - 4*pixelPitch):(0.25*pixelPitch):(maxOffset_org + 4*pixelPitch);
Z_stack = zeros(size(maxOffsetVector));
for i = 1:size(maxOffsetVector,2)
    Z_modified = (l_dis*D_half)/(maxOffsetVector(1,i) - D_half);
    Z_modified_pixel = Z_modified/pixelPitch; % 注意最后的还得是以像素为单位，于是求解到长度后最后再除以pixelPitch.
    
    if abs(Z_modified) < 100*l_dis % 防止偏移太大，使新生成的Z_modified坐标趋向于无限大 % modified by JDY 20190227
        Z_stack(1,i) = Z_modified_pixel;
    end
end
X_offsetVector = linspace((-4*pixelPitch*abs(Z_modified/l_dis)),(4*pixelPitch*abs(Z_modified/l_dis)),33);
X_stack = zeros(size(X_offsetVector));
for j = 1:size(X_offsetVector,2)
    X_modified = X_org + X_offsetVector(1,j);
    X_modified_pixel = X_modified/pixelPitch; % 注意最后的还得是以像素为单位，于是求解到长度后最后再除以pixelPitch.
    
    X_stack(1,j) = X_modified_pixel;
end
% 深度方向上的最小变化步长由视差来定，视差的范围被分到
end