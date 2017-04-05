clear;clc;cla;close all
% path = 'C:/MATLAB/libra/DynamicSLAM/rgbd_dataset_freiburg1_xyz/';
path = '';
fn = importdata([path 'rgbd.txt']);
cur = [];
ref = [];
depth = [];
pose = [];
traj = [];
trajT = [];
ts = [];
kf_list_x = [];
kf_list_y = [];
kf_list_depth = [];
KF_LIST_LENGTH = 5;
for i = 1:length(fn)  
    disp(['Process frame ' num2str(i)]);
    fnrgb = fn{i}(19:43);
    fndepth = fn{i}(63:89);
    
    cur =  im2double(imread([path fnrgb]));
    curdepth =  imread([path fndepth]);
    cur = rgb2gray(cur);
    [curgradientx, curgradienty] = imgradientxy(cur);
%     figure; imshow(cur);
%     figure;  imshow(curgradientx);
%     figure;  imshow(curgradienty);
    if ~isempty(ref)
        tic
%         T = estimate_ds_dense(ref,cur,refdepth,curdepth,T);%这里的T可以处理一下
    if size(kf_list_depth,2)/640 ==1
            T = estimate_ds_dense_gradient(cur,refgradientx,refgradienty,curgradientx,curgradienty,refdepth,curdepth,T);
            kf_list_x = [kf_list_x curgradientx];
            kf_list_y = [kf_list_y curgradienty];
            kf_list_depth = [kf_list_depth curdepth];
            trajT = [trajT pose];
            if size(kf_list_x,2)/640 >5
                kf_list_x(:,1:640) = [];
                kf_list_y(:,1:640) = [];
                kf_list_depth(:,1:640) = [];
                trajT(:, 1:4) = [];
            end
            size(trajT,2)/4
    else
            kf_ref_x = kf_list_x(1:480, 1:640);
            kf_ref_y = kf_list_y(1:480, 1:640);
            kf_ref_depth = kf_list_depth(1:480, 1:640);
            Tkf_r = trajT(:,1:4)' * trajT(:,size(trajT,2)-3:size(trajT,2));
            T = estimate_ds_dense_gradient_Joint(cur,refgradientx,refgradienty,kf_ref_x, kf_ref_y,curgradientx,curgradienty,refdepth,kf_ref_depth,curdepth,T,Tkf_r);
    end
%         T = estimate_ds_ic_dense(ref,cur,refdepth,T);
        toc
        if norm(T(1:3,4))>0.15 || norm(T(1:3,1:3)-eye(3),'fro')>0.15  % returns the Frobenius norm of X.
            kf = true;
            tic
%             T = estimate_icp(ref,cur,refdepth,curdepth,T);
%             T = estimate_ds_dense(ref,cur,refdepth,curdepth,T);
%             T = estimate_ds_dense_gradient(cur,refgradientx,refgradienty,curgradientx,curgradienty,refdepth,curdepth,T);
%             T = estimate_ds_ic_dense_kf(ref,cur,refdepth,curdepth,T);
%             T = estimate_ds_ic_kf(ref,cur,refdepth,curdepth,T);
            toc
        end
        cpose = pose*T;
    else   % first pose
        T = eye(4);
        pose = eye(4);
        cpose = pose;
        kf = true;
    end
    if kf   %kf 之后更新pose
        ts = [ts; fn{i}(1:17)];
        pose = pose*T;
        pose(1:3,1:3) = quat2rotm(rotm2quat(pose(1:3,1:3)));
        T = eye(4);
        kf = false;
        ref = cur;
        refgradientx=curgradientx;
        refgradienty = curgradienty;
    	refdepth = curdepth;
        q = rotm2quat(pose(1:3,1:3));
        plot3(traj(end,1),traj(end,2),traj(end,3),'ro');
        traj = [traj;pose(1:3,4).' q(2:4) q(1)];

    end
    plot3(cpose(1,4),cpose(2,4),cpose(3,4),'b*');hold on;
    axis equal
    drawnow
end
temp = [ts repmat(' ',size(ts,1),1) num2str(traj)];
dlmwrite('res.txt',temp,'delimiter','');
