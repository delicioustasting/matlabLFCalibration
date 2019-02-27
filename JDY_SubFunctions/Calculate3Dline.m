function [Point1_On3Dline, Point2_On3Dline, L_3Dline] =Calculate3Dline(lineStack_h, centerStack_h, k_h,l_dis,pixelPitch)
%%%%%%%%%%%%%% calculate the 3D line in camera image space
L_plane = zeros(k_h,4);% ƽ�淽�� L = [A;B;C;D];
for i=1:k_h
    %%%ÿ�������ص�linefeather��Ӧ��3D�ռ��У���XYZ����ƽ��ķ��̣����е�XYZ�������С����PixelPitchΪ��λ����
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
    L_plane(i,4) = (-1)*((lineStack_h(1:2,i))')*centerStack_h(1:2,i); % centerStack_h������Ϊ��λ
    %---------- modified by JDY 20190225
end
%%% [A,B,C,D]*[X,Y,Z,1] = 0;�������ϣ����ϵ�������пռ���������������޹��飬��Ϊ3D line�ı����ʽ
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
