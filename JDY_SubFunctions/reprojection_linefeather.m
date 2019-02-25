function [lineStack_h2, k2_h, lineStack_reproject]=...
    reprojection_linefeather(Point1_On3Dline, Point2_On3Dline, l_dis, pixelPitch, ...
    centerStack_h, lineStack_h, k_h, radius)
%figure;imshow(uint8(CaliImg));hold on;
%b=[-radius+0.5,radius-0.5];
lineStack_h2 = zeros(3, 3000);
lineStack_reproject = zeros(3, size(lineStack_h,k_h));
filtered = false;
k2_h = 0;
for i=1:k_h
    current_center = centerStack_h(:,i);
    %%%计算line上的两个3D点，在某个微透镜图像中的2D投影点的坐标。注意：该坐标为全局坐标，而每个
    %%%宏像素的linefeather的方程，使用的是以每个宏像素的中心为原点的局部坐标。
    %%%并且注意，这里方程中所有量的单位仍然是以pixel为单位，于是后面微透镜同sensor的距离，仍然要
    %%%换算成以pixel为单位，即l_dis/pixelPitch
    Point1_On2Dline = [current_center;0] + ((l_dis/pixelPitch)/Point1_On3Dline(3,1))*...
        ([current_center;0] - Point1_On3Dline(1:3,1));
    Point2_On2Dline = [current_center;0] + ((l_dis/pixelPitch)/Point2_On3Dline(3,1))*...
        ([current_center;0] - Point2_On3Dline(1:3,1));
    %%%利用全局坐标绘图
    %{
    line_feather = cross([Point1_On2Dline(1:2,1);1], [Point2_On2Dline(1:2,1);1]);    
    boundary=[b;-(b*line_feather(1,i)+line_feather(3,i))/line_feather(2,i)];
    plot(boundary(1,:)+centerStack_h(1,i),boundary(2,:)+centerStack_h(2,i),'g-','LineWidth',1);
    %}
    %%%计算局部坐标
    Point1_2Dlocal = Point1_On2Dline(1:2,1) - current_center;
    Point2_2Dlocal = Point2_On2Dline(1:2,1) - current_center;
    line_feather = cross([Point1_2Dlocal;1], [Point2_2Dlocal;1]); 
    lineStack_reproject(:,i) = line_feather;
    %{
    boundary=[b;-(b*line_feather(1,1)+line_feather(3,1))/line_feather(2,1)];
    plot(boundary(1,:)+centerStack_h(1,i),boundary(2,:)+centerStack_h(2,i),'b-','LineWidth',1);
    %}
    %plot(current_center(1,1),current_center(2,1),'y.','MarkerSize',5);
    %%%筛选，原始linefeather的斜率如果比重投影的linefeather的斜率，夹角大于5度，就不要了。
    if filtered
    if abs(lineStack_h(1:2,i)'*line_feather(1:2,1))/(norm(lineStack_h(1:2,i))*norm(line_feather(1:2,1)))...
            >= cos(5*pi/180)
            % 两条线的方向向量的夹角(即内积除以模值，即夹角的余弦)，也即两条线的法向量的夹角
        k2_h = k2_h + 1;
        lineStack_h2(:,k2_h) = lineStack_h(:, i);
    end
    end
end
end