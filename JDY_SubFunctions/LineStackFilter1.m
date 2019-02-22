function [lineStack_h_propotion_filter, centerStack_h_propotion_filter, num_PointsStack] = ...
    LineStackFilter1(Selected_corner_org, Selected_corner_left,...
    lineStack_h,centerStack_h,k_h,...
    l_dis, CaliImg, d_img)
%%% Paras of line which links left point and right point
Line_SubImg = zeros(3,1);
Line_SubImg_k = (Selected_corner_left(2,1) - Selected_corner_org(2,1)) /...
    (Selected_corner_left(1,1) - Selected_corner_org(1,1));
Line_SubImg_b = -Line_SubImg_k*Selected_corner_org(1,1)+Selected_corner_org(2,1);
%{
    figure;imshow(CenterSubImg);hold on;
    x_temp = round(Selected_corner_left(1,1)):1:round(Selected_corner_org(1,1));
    y_temp = Line_SubImg_k*x_temp + Line_SubImg_b;
    plot(x_temp,y_temp,'g-','LineWidth',1);    
%}
%%% Line in raw image
RawSelected_corner_org = Selected_corner_org*d_img;
RawSelected_corner_left = Selected_corner_left*d_img;
Line_RawImg = zeros(3,1);
Line_RawImg_k = (RawSelected_corner_left(2,1) - RawSelected_corner_org(2,1)) /...
    (RawSelected_corner_left(1,1) - RawSelected_corner_org(1,1));
Line_RawImg_b = -Line_RawImg_k*RawSelected_corner_org(1,1)+RawSelected_corner_org(2,1);
%{
    figure;imshow(CaliImg);hold on;
    x_temp = round(RawSelected_corner_left(1,1)):1:round(RawSelected_corner_org(1,1));
    y_temp = Line_RawImg_k*x_temp + Line_RawImg_b;
    plot(x_temp,y_temp,'g-','LineWidth',1);    
%}
%%selected 2D point
Points3D_Stack = zeros(3,k_h);num_PointsStack = 0;
lineStack_h_propotion_filter = zeros(3,k_h);
centerStack_h_propotion_filter = zeros(2,k_h);
Points3D_Proportion = zeros(1,k_h);
figure;imshow(CaliImg);hold on;
for i = 1:k_h % 当前这里面的变量还是以像素为单位的吧
    X_temp = centerStack_h(1,i);
    %%%
    Y_3D = Line_RawImg_k*X_temp + Line_RawImg_b;
    Y_center = centerStack_h(2,i);
    x_local = X_temp - X_temp;
    y_local = (-lineStack_h(1,i)*x_local-lineStack_h(3,i))/lineStack_h(2,i);
    Y_2D = y_local + Y_center;
    %%%
    Proportion = ((Y_2D - Y_3D)/(Y_center - Y_3D));
    disp(Proportion);
    if Proportion > 0 && abs(Proportion -1) > 0.05 &&...
            abs(Y_center - Y_2D) >= 1.6
        
        Z_temp = l_dis /(Proportion - 1);
        num_PointsStack = num_PointsStack + 1;
        Points3D_Stack(:,num_PointsStack) = [X_temp;Y_3D;Z_temp];
        lineStack_h_propotion_filter(:,num_PointsStack) = lineStack_h(:,i);
        centerStack_h_propotion_filter(:,num_PointsStack) = centerStack_h(:,i);
        Points3D_Proportion(1,num_PointsStack) = Proportion;
        %%%
        if 1
            plot(X_temp, Y_3D, 'ro');
            plot(X_temp, Y_center, 'bo');
            plot(X_temp, Y_2D, 'go');
        end
        if num_PointsStack == 23 || num_PointsStack ==24 ...
                || num_PointsStack==36 || num_PointsStack==37 ...
                || num_PointsStack ==47
            plot(X_temp, Y_3D, 'r*');
            plot(X_temp, Y_center, 'b*');
            plot(X_temp, Y_2D, 'g*');
        end
        %%%
    end
end
%{
figure;plot(Points3D_Stack(3,:));
figure;plot(Points3D_Proportion);
%}
disp('Wait');close all;
end