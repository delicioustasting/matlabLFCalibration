function [Point_left, Point_right]= ExportLine3DLocalTerminals(...
    Point1_On3Dline, Point2_On3Dline, Selected_corner_org, Selected_corner_left, d_img, Plot_Flag)
%%% Line in raw image
RawSelected_corner_org = Selected_corner_org*d_img;
RawSelected_corner_left = Selected_corner_left*d_img;

%%%%%
Point_left = zeros(3,1);
Point_left(1,1) = RawSelected_corner_left(1,1);
lambda_left = (Point_left(1,1) - Point1_On3Dline(1,1))/(Point2_On3Dline(1,1) - Point1_On3Dline(1,1));
Point_left = Point1_On3Dline + lambda_left*(Point2_On3Dline - Point1_On3Dline);
%%%%%
Point_right = zeros(3,1);
Point_right(1,1) = RawSelected_corner_org(1,1);
lambda_right = (Point_right(1,1) - Point1_On3Dline(1,1))/(Point2_On3Dline(1,1) - Point1_On3Dline(1,1));
Point_right = Point1_On3Dline + lambda_right*(Point2_On3Dline - Point1_On3Dline);
%%%%%
if Plot_Flag == true
    plot([Point_left(1,1), Point_right(1,1)],...
        [Point_left(2,1), Point_right(2,1)],...
        'b-','LineWidth',1);
end
%{
plot([Point1_On3Dline(1,1), Point2_On3Dline(1,1),Point1_On3Dline(1,1)+2*(Point2_On3Dline(1,1)-Point1_On3Dline(1,1))],...
[Point1_On3Dline(2,1), Point2_On3Dline(2,1),Point1_On3Dline(2,1)+2*(Point2_On3Dline(2,1)-Point1_On3Dline(2,1))],...
'b-','LineWidth',1);
%}
end