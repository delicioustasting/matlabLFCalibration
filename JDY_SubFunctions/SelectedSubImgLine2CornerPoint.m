function [Selected_corner_org, Selected_corner_left]=...
    SelectedSubImgLine2CornerPoint(CornerIndexList, corner, HV_flag)
if HV_flag == 'h'
%%selected corner point
for i = 1:size(corner,2)
    if (corner(1,i) == CornerIndexList(3,1)) &&...
            (corner(2,i) == CornerIndexList(3,2))
        Selected_corner_org = corner(3:4,i);
    end
    if (corner(1,i) == (CornerIndexList(3,1)-1)) &&...
            (corner(2,i) == CornerIndexList(3,2))
        Selected_corner_left = corner(3:4,i);
    end
end
end
if HV_flag == 'v'
%%selected corner point
for i = 1:size(corner,2)
    if (corner(1,i) == CornerIndexList(3,1)) &&...
            (corner(2,i) == CornerIndexList(3,2))
        Selected_corner_org = corner(3:4,i);
    end
    if (corner(1,i) == CornerIndexList(3,1)) &&...
            (corner(2,i) == (CornerIndexList(3,2) - 1))
        Selected_corner_left = corner(3:4,i);
    end
end
end
%{
figure;imshow(CenterSubImg);hold on;
plot(Selected_corner_org(1,1), Selected_corner_org(2,1), 'ro','MarkerSize', 5);hold on;
plot(Selected_corner_left(1,1), Selected_corner_left(2,1), 'bo','MarkerSize', 5);hold on;
hold off;
figure;imshow(CaliImg);hold on;
d_img = 14.01;
plot(Selected_corner_org_h(1,1)*d_img, Selected_corner_org_h(2,1)*d_img, 'ro','MarkerSize', 5);hold on;
plot(Selected_corner_left_h(1,1)*d_img, Selected_corner_left_h(2,1)*d_img, 'bo','MarkerSize', 5);hold on;

plot(Selected_corner_org_v(1,1)*d_img, Selected_corner_org_v(2,1)*d_img, 'yo','MarkerSize', 5);hold on;
plot(Selected_corner_Up_v(1,1)*d_img, Selected_corner_Up_v(2,1)*d_img, 'go','MarkerSize', 5);hold on;
hold off;
%}
end