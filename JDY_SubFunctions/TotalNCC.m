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
        - (YGrid_coords - YGrid_integer)*lineStack_output(2,i);% ��line�Ĳ����䵽�ԣ��µ��������Ϊԭ�������ϵ��'
    
    % ģ���������ò��ֻ�ܲ���ԭ�������м�������ص�ģ�壬��Ϊ�ú����е�����ֵ��ò��ֻ������������ֵ
    % ����ԭ�㲻����������������ʱ������(2*6+1)*(2*6+1)= 13*13����������ֵ����ΪС���ˣ��ú��������á�
    % ���ǣ��Ȱ���centerStackΪԭ��ľֲ������е�line_param��ɣ��������Ϊԭ���line_param
    template=LinearTemplate(radius,line_param); % �ú��������2*��radius+1��+1�߳���С�������Σ�Ӧ�ù��˰�
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

function template=LinearTemplate(size_half,line_param) % size_halfĬ�ϵ���radius-1������

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