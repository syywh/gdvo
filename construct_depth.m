function [ depth ] = construct_depth( left, right, stereoParams, baseline_focal )

% [J1, J2] = rectifyStereoImages(left, right, stereoParams);

% imshow(J1);
% figure(1);
% imshow(J2);
% imshowpair(left, right);
% figure(2);
% imshowpair(J1,J2);
% disparityMap = disparity(J1, J2 ,'Method', 'BlockMatching');
% disparityMap = disparity(J1, J2);
% disparityMap = disparity(left, right,'Method', 'BlockMatching');
disparityRange = [0 96];
% disparityMap = disparity(left, right,'BlockSize',5,'DisparityRange',disparityRange, 'DistanceThreshold',100,'UniquenessThreshold',30);
disparityMap = disparity(left, right,'BlockSize',5,'DisparityRange',disparityRange, 'UniquenessThreshold',70,'ContrastThreshold',0.1);
% figure(3);
% imshow(disparityMap, [0, 64]);

depth = baseline_focal./disparityMap;

% pointCloud = double(reconstructScene(disparityMap, stereoParams));
% depth = pointCloud(:, :, 3);
% mask = repmat(Z > 10 & Z < 50, [1, 1, 3]);
% mask = (depth>10 & depth<50);
% J1(~mask) = 0;
% figure(4);
% imshowpair(left, depth);
% figure(2)
% imshow(J1, 'InitialMagnification', 50);
end

