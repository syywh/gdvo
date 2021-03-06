clear;clc;cla;close all
% path = 'C:/MATLAB/libra/DynamicSLAM/rgbd_dataset_freiburg1_xyz/';
% path = '';
path = '/home/dxq/data/datasets/kitti_odo/data_odometry_gray/sequences/05';
path_left = [path '/image_0/'];
path_right = [path '/image_1/'];
% fn = importdata([path 'rgbd.txt']);
fn = dir(path_left);

% f_true = fopen('../00.txt', 'rt');
% 
% TT = [];
% for i = 1:3
%     T = zeros(3,4);
%     for j = 1:12
%          fprintf(f_true, '%g',T(j));
%     end
%     TT = [TT T];
% end

cur_left = [];      cur_right = [];
ref_left = [];      ref_right = [];
depth = [];
pose = [];
traj = [];
trajT = [];
trajframeT = [];
ts = [];
kf_list_x = [];
kf_list_y = [];
kf_list_depth = [];
KF_LIST_LENGTH = 7;

KL = [9.799200e+02 0.000000e+00 6.900000e+02 0.000000e+00 9.741183e+02 2.486443e+02 0.000000e+00 0.000000e+00 1.000000e+00];
KL = reshape(KL, 3,3);
KR=[9.903522e+02 0.000000e+00 7.020000e+02 0.000000e+00 9.855674e+02 2.607319e+02 0.000000e+00 0.000000e+00 1.000000e+00];
KR = reshape(KR, 3,3);
TL = [-9.251859e-17 8.326673e-17 -7.401487e-17];
TR = [-5.370000e-01 5.964270e-03 -1.274584e-02];
delT =[-3.861448e+02/9.799200e+02  0 0] ;
%  delT =[0  0 0] ;

cameraParamsL = cameraParameters('IntrinsicMatrix', KL);
cameraParamsR = cameraParameters('IntrinsicMatrix', KR);
RL = [9.999454e-01 7.259129e-03 -7.519551e-03 -7.292213e-03 9.999638e-01 -4.381729e-03 7.487471e-03 4.436324e-03 9.999621e-01];
RR = [9.996568e-01 -1.110284e-02 2.372712e-02 1.099810e-02 9.999292e-01 4.539964e-03 -2.377585e-02 -4.277453e-03 9.997082e-01];
% RL = [1.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00];
% RR = [9.993513e-01 1.860866e-02 -3.083487e-02 -1.887662e-02 9.997863e-01 -8.421873e-03 3.067156e-02 8.998467e-03 9.994890e-01];
RL = reshape(RL,3,3);
RR = reshape(RR,3,3);
rotateR2toR1 = RL'*RR
% stereoParams = stereoParameters(cameraParamsL, cameraParamsR, rotateR2toR1, (TR-TL)');
 stereoParams = stereoParameters(cameraParamsL, cameraParamsR, eye(3), delT');%已经是rectify的了

% for i = 3:length(fn)  
for i = 3:11
    disp(['Process frame ' num2str(i-3)]);
    cur_left =  im2double(imread([path_left fn(i).name]));
    cur_right =  im2double(imread([path_right fn(i).name]));
    [curdepth] = construct_depth(cur_left, cur_right, stereoParams, 3.861448000000e+02);
    curdepth = medfilt2(curdepth);
%     imshow(curdepth,[0 100])
    width = size(cur_left,1);%512
    height = size(cur_left,2);%1392
%     curdepth =  imread([path fndepth]);
%     cur_left = rgb2gray(cur_left);
    [curgradientx, curgradienty] = imgradientxy(cur_left);
%     figure; imshow(cur);
%     figure;  imshow(curgradientx);
%     figure;  imshow(curgradienty);
    if ~isempty(ref_left)
        tic

            kf_ref_x = kf_list_x(:, 1:height);
            kf_ref_y = kf_list_y(:, 1:height);
            kf_ref_depth = kf_list_depth(:, 1:height);
            Tkf_r = inv(trajT(:,1:4)) * trajT(:,size(trajT,2)-3:size(trajT,2));
            if size(trajT,2)/4<2
                initV = [0.999997800000000,0.000527262800000000,-0.00206693500000000,-0.0469029400000000,-0.000529650600000000,0.999999200000000,-0.00115486500000000,-0.0283992800000000,0.00206632400000000,0.00115595800000000,0.999997100000000,0.858694100000000];
                initV = reshape(initV,4,3)';
                initV(4,:) = zeros(1,4);
                initV(4,4) = 1;
                initV=eye(4);
%                 imshow(refdepth);
                T = estimate_ds_dense_gradient_Joint(cur_left,refgradientx,refgradienty,kf_ref_x, kf_ref_y,curgradientx,curgradienty,refdepth,kf_ref_depth,curdepth,initV,initV, true);
            else
                estimateT = inv(trajT(:,size(trajT,2)-7:size(trajT,2)-4)) * trajT(:,size(trajT,2)-3:size(trajT,2));
                Tc = trajT(:,size(trajT,2)-3:size(trajT,2))*estimateT;
                T = estimate_ds_dense_gradient_Joint(cur_left,refgradientx,refgradienty,kf_ref_x, kf_ref_y,curgradientx,curgradienty,refdepth,kf_ref_depth,curdepth,estimateT,Tkf_r,false);
            end
        toc
        if norm(T(1:3,4))>0.15 || norm(T(1:3,1:3)-eye(3),'fro')>0.15  % returns the Frobenius norm of X.
            kf = true;
        end
        
        cpose = cpose*T;
        ref_left = cur_left;
        refgradientx=curgradientx;
        refgradienty = curgradienty;
    	refdepth = curdepth;
        kf_list_x = [kf_list_x curgradientx];
        kf_list_y = [kf_list_y curgradienty];
        kf_list_depth = [kf_list_depth curdepth];
        trajT = [trajT cpose];
        if size(kf_list_x,2)/height >KF_LIST_LENGTH
            kf_list_x(:,1:height) = [];
            kf_list_y(:,1:height) = [];
            kf_list_depth(:,1:height) = [];
            trajT(:, 1:4) = [];
        end
    else   % first pose
        T = eye(4);
        pose = eye(4);
        cpose = pose;
        kf = true;
        ref_left = cur_left;
        refgradientx=curgradientx;
        refgradienty = curgradienty;
    	refdepth = curdepth;
        kf_list_x = [kf_list_x curgradientx];
        kf_list_y = [kf_list_y curgradienty];
        kf_list_depth = [kf_list_depth curdepth];
        trajT = [trajT cpose*T];
        if size(kf_list_x,2)/height >KF_LIST_LENGTH
            kf_list_x(:,1:height) = [];
            kf_list_y(:,1:height) = [];
            kf_list_depth(:,1:height) = [];
            trajT(:, 1:4) = [];
        end


    end
    trajframeT = [trajframeT cpose];
    if kf   %kf 之后更新pose
        pose = pose*T;
        pose(1:3,1:3) = quat2rotm(rotm2quat(pose(1:3,1:3)));
        T = eye(4);
        kf = false;
        q = rotm2quat(pose(1:3,1:3));
        traj = [traj;pose(1:3,4).' q(2:4) q(1)];
%         plot3(traj(end,1),traj(end,2),traj(end,3),'ro');
        traj = [traj;pose(1:3,4).' q(2:4) q(1)];


    end
%     plot3(cpose(1,4),cpose(2,4),cpose(3,4),'b*');hold on;
%     axis equal
%     drawnow
end
temp = [ts repmat(' ',size(ts,1),1) num2str(traj)];
dlmwrite('res.txt',temp,'delimiter','');
frame_num = size(trajframeT,2)/4;
trajframeT = reshape(trajframeT, 4,4 , frame_num);
ft = fopen('result.txt','wt');
T_real = load('../05.txt');
error_t = [];
for i = 1:frame_num
     tempT = trajframeT(:,:,i)';
     for j = 1:12
           fprintf(ft, '%g ', tempT(j));
     end
    fprintf(ft,'\n');
    real_t = [T_real(i,4) T_real(i,8) T_real(i,12)];
    error_t =  [error_t norm(real_t -  tempT(4,1:3))];
    
end
figure;
plot(error_t);




