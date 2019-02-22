function [lineStack_h, centerStack_h, lineStack_v, centerStack_v, k_h, k_v]...
    = SelectedCorner(CenterSubImg, CaliImg, CornerIndexList,...
    line_h, center_h, line_v, center_v, world_h, world_v)
%%%the initial location of the corner point
CornerList = [325, 150; 204, 240; 291, 241; 238, 157; 166, 331; 396, 243; 425, 144; 353, 71;];
d_img = 14.01; % Coarse distance between centers of adjacent macro images
if false
figure;imshow(CenterSubImg);hold on;
for i = 1:2%size(CornerList,1)
plot(CornerList(i,1), CornerList(i,2), 'ro','MarkerSize', 5);hold on;
end
hold off;
end
figure;imshow(CaliImg);hold on;
for i = 1:size(CornerList,1)
plot(CornerList(i,1)*d_img, CornerList(i,2)*d_img, 'ro','MarkerSize', 5);hold on;
end
%hold off;
%%%%%%%%%%%%%%
for i = 1:size(CornerIndexList,1)
    CornerIndex = CornerIndexList(i,:);
    %%% 局部要使用的数据结构初始化
    lineStack_h = zeros(3, 3000);
    lineStack_v = zeros(3, 3000);
    centerStack_h = zeros(2, 3000);
    centerStack_v = zeros(2, 3000);
    k_h = 0;k_v = 0; % index of each lineStack
    %%%
    for h_idx = 1:size(world_h, 2)
        if world_h(1, h_idx) == CornerIndex(1, 1) &&...
                world_h(2, h_idx) == CornerIndex(1, 2)
            k_h = k_h+1;
            lineStack_h(:,k_h) = line_h(:,h_idx);
            centerStack_h(:,k_h) = center_h(:,h_idx);            
        end
    end
    for v_idx = 1:size(world_v, 2)
        if world_v(1, v_idx) == CornerIndex(1, 1) &&...
                world_v(2, v_idx) == CornerIndex(1, 2)
             k_v = k_v+1;
            lineStack_v(:,k_v) = line_v(:,v_idx);
            centerStack_v(:,k_v) = center_v(:,v_idx);
           
        end
    end
end
disp('Wait');%close all;
%%%%%%%%% display line feather
radius = 7.0;
%figure;imshow(uint8(CaliImg));hold on;

    b=[-radius+0.5,radius-0.5];
    for i=1:k_h
        boundary=[b;-(b*lineStack_h(1,i)+lineStack_h(3,i))/lineStack_h(2,i)];
        plot(boundary(1,:)+centerStack_h(1,i),boundary(2,:)+centerStack_h(2,i),'g-','LineWidth',1);
    end
    for i=1:k_v
        boundary=[-(b*lineStack_v(2,i)+lineStack_v(3,i))/lineStack_v(1,i);b];
        plot(boundary(1,:)+centerStack_v(1,i),boundary(2,:)+centerStack_v(2,i),'g-','LineWidth',1);
    end
    plot(centerStack_h(1,1:k_h),centerStack_h(2,1:k_h),'r.','MarkerSize',5);
    plot(centerStack_v(1,1:k_v),centerStack_v(2,1:k_v),'r.','MarkerSize',5);
    hold off;
%%%%%%%%% display line feather:: over
disp('Wait');
end