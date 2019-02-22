function Plot_3DLineXYProjection(Point1_On3Dline, Point2_On3Dline, Selected_corner_org, Selected_corner_left, d_img)
%%% Line in raw image
RawSelected_corner_org = Selected_corner_org*d_img;
RawSelected_corner_left = Selected_corner_left*d_img;

Line_RawImg = zeros(3,1);
Line_RawImg_k = (Point1_On3Dline(2,1) - Point2_On3Dline(2,1)) /...
    (Point1_On3Dline(1,1) - Point2_On3Dline(1,1));
Line_RawImg_b = -Line_RawImg_k*Point2_On3Dline(1,1)+Point2_On3Dline(2,1);
%%%%%
Point_left = zeros(2,1);
Point_right = zeros(2,1);
Point_left(1,1) = RawSelected_corner_left(1,1);
Point_left(2,1) = Line_RawImg_k*Point_left(1,1) + Line_RawImg_b;
Point_right(1,1) = RawSelected_corner_org(1,1);
Point_right(2,1) = Line_RawImg_k*Point_right(1,1) + Line_RawImg_b;
%%%%%
plot([Point_left(1,1), Point_right(1,1)],...
[Point_left(2,1), Point_right(2,1)],...
'b-','LineWidth',1);
%{
plot([Point1_On3Dline(1,1), Point2_On3Dline(1,1),Point1_On3Dline(1,1)+2*(Point2_On3Dline(1,1)-Point1_On3Dline(1,1))],...
[Point1_On3Dline(2,1), Point2_On3Dline(2,1),Point1_On3Dline(2,1)+2*(Point2_On3Dline(2,1)-Point1_On3Dline(2,1))],...
'b-','LineWidth',1);
%}
end