function SetCameraParas
if exist ('CameraParas.mat')
    %load('CameraParas.mat');
    return
end
%%%%%========================================================input paras
%%% device :: sensor
pixelPitch = 1.4;%1.3999999999999999366e-06;
pixelWidth = 7728;
bitsPerPixel = 10;
pixelHeight = 5368;
%%% device::mla
mla_rotation  = -0.0025808198843151330948;
lensPitch = 20;%2.0000000000000001636e-05;
sensorOffset=[0,0,40];%[ -4.4266057014465332508e-06, 4.2581391334533687899e-06, 3.6999999999999998114e-05];
l_dis=sensorOffset(1,3); % l_dis, distance between MLA and sensor;L_dis, distance between MLA and mainlens
%%% device::Main lens
focalLength = 0.011832640061537134935;
opticalCenterOffset=[-5.8207762776874005795e-05,-1.9448425518930889666e-05];

%%%calibration image generation
%%%人为设置标定板黑白格子的边长
l_CaliUnit = round(pixelWidth*pixelPitch/16);
l_half_CaliUnit = l_CaliUnit*0.5;

%%% 
L_dis = 160*l_dis;
x_offset = pixelWidth*pixelPitch*0.5 + sensorOffset(1,1);
y_offset = pixelHeight*pixelPitch*0.5 + sensorOffset(1,2);
z_offset = -l_dis; % equals to sensorOffset(1,3)
radius = 7;
%%%%%%=======================================================input paras::over
save('CameraParas.mat');
end