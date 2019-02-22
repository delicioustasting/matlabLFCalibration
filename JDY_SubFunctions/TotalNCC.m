function NCC = TotalNCC(centerStack_h, lineStack_output,...
                    k_h, radius, pixelHeight, pixelWidth, CaliImg)
SumCross = 0;
SumSource = 0;
SumTemplate = 0;
%%%%%%%%%%%%%%%
for i = 1:k_h
    XGrid_coords = centerStack_h(1,i);
    YGrid_coords = centerStack_h(2,i);
    XGrid_integer= round(XGrid_coords);
    YGrid_integer= round(YGrid_coords);
    line_param = zeros(3,1);
    line_param(1:2,1) = lineStack_output(1:2,i);    
    line_param(3,1) = lineStack_output(3,i)...
        - (XGrid_coords - XGrid_integer)*lineStack_output(1,i)...
        - (YGrid_coords - YGrid_integer)*lineStack_output(2,i);% 把line的参数变到以，新的整数格点为原点的坐标系上'
    
    % 模板产生函数貌似只能产生原点在最中间的整像素的模板，因为该函数中的坐标值，貌似只能是整数格点的值
    % 当，原点不设在整数格点的像素时，所有(2*6+1)*(2*6+1)= 13*13个格点的坐标值均成为小数了，该函数不适用。
    % 于是，先把以centerStack为原点的局部坐标中的line_param变成，以整格点为原点的line_param
    template=LinearTemplate(radius,line_param); % 该函数会产生2*（radius+1）+1边长大小的正方形，应该够了吧
    source = zeros(size(template));
    for x_p = -radius-1:1:radius+1
        for y_p = -radius-1:1:radius+1
            X_coords= x_p + XGrid_integer;
            Y_coords= y_p + YGrid_integer;
            
            if (((X_coords -XGrid_coords)*(X_coords -XGrid_coords)...
                    +(Y_coords -YGrid_coords)*(Y_coords -YGrid_coords))...
                    <= (radius*radius))&&...
                    (X_coords>=1)&&(X_coords<=pixelWidth)&&...
                    (Y_coords>=1)&&(Y_coords<=pixelHeight)
                %%%%%%%%%%%% calculate the pixel location coordinates
                source(y_p+radius+1,x_p+radius+1) = CaliImg(Y_coords, X_coords);
            end
        end
    end
    %%%
     SumCross = SumCross + sum(sum(source.*template));
    SumSource = SumSource + sum(sum(source.*source));
    SumTemplate = SumTemplate + sum(sum(template.*template));
end
%%%%%%%%%%%%%%%
NCC = SumCross/sqrt(SumSource * SumTemplate);

end

function template=LinearTemplate(size_half,line_param) % size_half默认等于radius-1，？；

size1=size_half*2+1;
template=zeros(size1,size1);
[I,J]=meshgrid(-size_half-0.5:size_half+0.5,-size_half-0.5:size_half+0.5);
d=line_param(1)*I+line_param(2)*J+line_param(3);
d(d<0)=-1;
d(d>=0)=1;
intersection_x=(-line_param(2)*(-size_half-0.5:size_half+0.5)-line_param(3))/line_param(1);
intersection_y=(-line_param(1)*(-size_half-0.5:size_half+0.5)-line_param(3))/line_param(2);
for j=1:size1
    for i=1:size1
        sum=d(j,i)+d(j,i+1)+d(j+1,i)+d(j+1,i+1);
        if sum==4
            template(j,i)=1;
        elseif sum==-4
            template(j,i)=0;
        elseif abs(sum)==2
            if d(j,i)*sum<0
                temp=abs(intersection_x(j)+size_half+1.5-i)*abs(intersection_y(i)+size_half+1.5-j)*0.5;
            elseif d(j,i+1)*sum<0
                temp=(1-abs(intersection_x(j)+size_half+1.5-i))*abs(intersection_y(i+1)+size_half+1.5-j)*0.5;
            elseif d(j+1,i)*sum<0
                temp=abs(intersection_x(j+1)+size_half+1.5-i)*(1-abs(intersection_y(i)+size_half+1.5-j))*0.5;
            else
                temp=(1-abs(intersection_x(j+1)+size_half+1.5-i))*(1-abs(intersection_y(i+1)+size_half+1.5-j))*0.5;
            end
            if sum>0
                template(j,i)=1-temp;
            else
                template(j,i)=temp;
            end
        else
            if d(j,i)+d(j,i+1)==0
                temp=(abs(intersection_x(j)+size_half+1.5-i)+abs(intersection_x(j+1)+size_half+1.5-i))*0.5;
            else
                temp=(abs(intersection_y(i)+size_half+1.5-j)+abs(intersection_y(i+1)+size_half+1.5-j))*0.5;
            end
            if d(j,i)>0
                template(j,i)=temp;
            else
                template(j,i)=1-temp;
            end
        end
    end
end

end