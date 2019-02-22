function LoadNecessaryData
if exist ('PreComputedData.mat')    
    return
end
%%%%%%%%%%%%%%%%%%%%%%%% input necessary data {{
%%% centers of the macro images calculated by white images
load('G:\ZSP_white and cali plane\JinBo20190118\microlens_center_list.mat','center_list');
%%% calibration images, center view sub-image
CaliImg = imread('G:\ZSP_white and cali plane\JinBo20190118\IMG_1333-1.png');
CenterSubImg = imread('G:\ZSP_white and cali plane\JinBo20190118\CI_IMG_1333-1.bmp');
%figure;imshow(CenterSubImg);
%%% line feather info
load('G:\ZSP_white and cali plane\JinBo20190118\L_IMG_1333-1.mat');
%%% Corner info of calibration image
load('G:\ZSP_white and cali plane\JinBo20190118\CI_IMG_1333-1.mat');
%%%%%%%%%%%%%%%%%%%%%%%% input necessary data }}
save('PreComputedData.mat');
end