function Plot_LineFeather_h(radius, lineStack_reproject, k_h, centerStack_h,color)

b=[-radius+0.5,radius-0.5];
for i = 1:k_h
boundary=[b;-(b*lineStack_reproject(1,i)+lineStack_reproject(3,i))/lineStack_reproject(2,i)];
plot(boundary(1,:)+centerStack_h(1,i),boundary(2,:)+centerStack_h(2,i),color,'LineWidth',1);
end
plot(centerStack_h(1,1:k_h),centerStack_h(2,1:k_h),'r.','MarkerSize',5);

end