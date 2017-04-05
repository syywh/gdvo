path = '/home/dxq/data/datasets/kitti/2011_09_26/2011_09_26_drive_0001_extract/';

leftImageFile = [path 'image_00/data/0000000000.png'];
rightImageFile = [path 'image_01/data/0000000000.png'];

leftImage = imread(leftImageFile);
rightImage = imread(rightImageFile);

PL = [7.215377e+02 0.000000e+00 6.095593e+02 0.000000e+00;
    0.000000e+00 7.215377e+02 1.728540e+02 0.000000e+00;
    0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00];
RL = [9.999239e-01 9.837760e-03 -7.445048e-03;
    -9.869795e-03 9.999421e-01 -4.278459e-03;
    7.402527e-03 4.351614e-03 9.999631e-01];
TL = [2.573699e-16 -1.059758e-16 1.614870e-16];
% KL = PL(:, 1:3) * RL';
KL = [9.842439e+02 0.000000e+00 6.900000e+02 0.000000e+00 9.808141e+02 2.331966e+02 0.000000e+00 0.000000e+00 1.000000e+00];
KL = reshape(KL, 3,3)

PR = [7.215377e+02 0.000000e+00 6.095593e+02 -3.875744e+02;
   0.000000e+00 7.215377e+02 1.728540e+02 0.000000e+00;
   0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00];
RR = [9.996878e-01 -8.976826e-03 2.331651e-02;
    8.876121e-03 9.999508e-01 4.418952e-03;
   -2.335503e-02 -4.210612e-03 9.997184e-01];
TR = [-5.370000e-01 4.822061e-03 -1.252488e-02];
% KR = PR(:, 1:3) * RR';
KR=[9.895267e+02 0.000000e+00 7.020000e+02 0.000000e+00 9.878386e+02 2.455590e+02 0.000000e+00 0.000000e+00 1.000000e+00];
KR = reshape(KR, 3,3);

cameraParamsL = cameraParameters('IntrinsicMatrix', KL);
cameraParamsR = cameraParameters('IntrinsicMatrix', KR);

rotateR2toR1 = RL'*RR;

stereoParams = stereoParameters(cameraParamsL, cameraParamsR, rotateR2toR1, (TR-TL)');
% load('webcamsSceneReconstruction.mat');
[J1, J2] = rectifyStereoImages(leftImage, rightImage, stereoParams);
disparityMap = disparity(J1, J2 ,'Method', 'BlockMatching');

% figure; imshow(cat(3, J1(:,:,1), J2(:,:)), 'InitialMagnification', 50);
figure; imshow(cat(2, J1(:,:,1), J2(:,:)), 'InitialMagnification', 50);
% figure; imshow(disparityMap, [0, 64], 'InitialMagnification', 50);
% reconstructScene(disparityMap, stereoParams);
% [J1, J2] = rectifyStereoImages(leftImage, rightImage, stereoParams);
% disparityMap = disparity(rgb2gray(J1), rgb2gray(J2));
% figure; imshow(cat(3, J1(:,:,1), J2(:,:,2:3)), 'InitialMagnification', 20);

imshow(disparityMap, [0, 64]);
title('Disparity Map');
colormap jet
colorbar
pointCloud = reconstructScene(disparityMap, stereoParams);
Z = pointCloud(:, :, 3);
mask = repmat(Z > 10 & Z < 50, [1, 1, 3]);
J1(~mask) = 0;
figure
imshow(J1, 'InitialMagnification', 50);
