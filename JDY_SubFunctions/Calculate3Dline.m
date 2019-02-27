function [Point1_On3Dline, Point2_On3Dline, L_3Dline] =Calculate3Dline(lineStack_h, centerStack_h, k_h,l_dis,pixelPitch)
%%%%%%%%%%%%%% calculate the 3D line in camera image space
L_plane = zeros(k_h,4);% 平面方程 L = [A;B;C;D];
for i=1:k_h
    %%%每个宏像素的linefeather对应的3D空间中，点XYZ，的平面的方程，其中点XYZ的坐标大小，以PixelPitch为单位长度
    %---------- modified by JDY 20190225 
    %{
    L_plane(i,1) = lineStack_h(1,i);
    L_plane(i,2) = lineStack_h(2,i);
    L_plane(i,3) = lineStack_h(3,i);
    L_plane(i,4) = (-1)*((lineStack_h(:,i))')*[centerStack_h(:,i);(l_dis/pixelPitch + 1)];
    %}
    %---------- modified by JDY 20190225
    L_plane(i,1) = lineStack_h(1,i);
    L_plane(i,2) = lineStack_h(2,i);
    L_plane(i,3) =  - lineStack_h(3,i)*(pixelPitch/l_dis);
    L_plane(i,4) = (-1)*((lineStack_h(1:2,i))')*centerStack_h(1:2,i); % centerStack_h以像素为单位
    %---------- modified by JDY 20190225
end
%%% [A,B,C,D]*[X,Y,Z,1] = 0;点在线上，求解系数矩阵行空间的两个极大线性无关组，作为3D line的表达形式
L_3Dline = zeros(2,4);
Point1_On3Dline = zeros(4,1);
Point2_On3Dline = zeros(4,1);

[U,S,V] = svd(L_plane);
L_3Dline(1,:)=V(:,1)';
L_3Dline(2,:)=V(:,2)';
Point1_On3Dline(:,1)=V(:,3);
Point1_On3Dline=Point1_On3Dline/Point1_On3Dline(end);
Point2_On3Dline(:,1)=V(:,4);
Point2_On3Dline=Point2_On3Dline/Point2_On3Dline(end);
end
