
path_image = '/home/dxq/data/datasets/kitti_odo/data_odometry_gray/sequences/00';
path_velo = '/home/dxq/data/datasets/kitti_odo/dataset_odo_velodyne/sequences/00/velodyne/';
path_left = [path_image '/image_0/'];
path_right = [path_image '/image_1/'];

fn_image = dir(path_left);
fn_velo = dir(path_velo);

% T = [4.276802385584e-04 -9.999672484946e-01 -8.084491683471e-03 -1.198459927713e-02 -7.210626507497e-03 8.081198471645e-03 -9.999413164504e-01 -5.403984729748e-02 9.999738645903e-01 4.859485810390e-04 -7.206933692422e-03 -2.921968648686e-01];
% T = reshape(T,4,3);
% Tr = eye(4);
% Tr(1:3, :) = T';
Tr = [
    0.0080   -1.0000   -0.0008   -0.0138
   -0.0028    0.0008   -1.0000   -0.0554
    1.0000    0.0080   -0.0028   -0.2919
         0         0         0    1.0000];

P0P = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02 0.000000000000e+00 0.000000000000e+00 7.188560000000e+02 1.852157000000e+02 0.000000000000e+00 0.000000000000e+00 0.000000000000e+00 1.000000000000e+00 0.000000000000e+00];
P0P = reshape(P0P,4,3);
P0 = eye(4);
P0(1:3,:) = P0P';

for i = 3:3
    cur_left =  im2double(imread([path_left fn_image(i).name]));
    cur_right =  im2double(imread([path_right fn_image(i).name]));
    [curdepth] = construct_depth(cur_left, cur_right, 1, 3.861448000000e+02);
    fid = fopen([path_velo fn_velo(i).name],'rb');
    velo = fread(fid,[4 inf],'single')';
    velo = velo(1:10:end,:); % remove every 5th point for display speed
    pad = ones(size(velo,1),1);
    velo(:,4) = pad;
    fclose(fid);
    
    width = size(cur_left,1);
    height = size(cur_left,2);
    image = P0P'*Tr*velo';
    image_depth = Tr*velo';
%     figure;
%     for i = 1: size(velo,1)
% %         plot3(image_depth(1,i), image_depth(2,i), image_depth(3,i), '.b');
%           plot3(velo(i,1), velo(i,2), velo(i,3), '.b');
%         hold on;
%     end
    uNew = round(image(1,:)./image(3,:))+1;
    vNew = round(image(2,:)./image(3,:))+1;

    far = 80;
    valid =   (uNew <= height) & uNew > 0 & vNew <= width & vNew > 0 & image_depth(3,:)>0 &image_depth(3,:)<far ;
    sum(valid)
    
    newimage = zeros(width,height);
    
    num = 0;
    for m = 1:size(valid,2)
            if ~valid(m)
                continue;
            end   

            if newimage(vNew(m), uNew(m))> image(3,m) || newimage(vNew(m), uNew(m)) == 0
                newimage(vNew(m), uNew(m)) = image_depth(3,m);   
                
            end
    end
%     a = (sum(newimage,2))
figure(1);
    imshowpair(curdepth, newimage);
    
    curdd = curdepth;
    for i = 1:size(curdepth,1)
        for j  =1:size(curdepth,2)
            if(newimage(i,j) == 0 || curdepth(i,j) == 0)
                curdd(i,j) = 0;
                newimage(i,j) = 0;
            end

            if(curdepth(i,j) < 0 || curdepth(i,j) >far)
                curdd(i,j) =0;
                newimage(i,j) = 0;
            end
%             if(abs(curdepth(i,j) - newimage(i,j)) > 10 )
%                 curdd(i,j) =0;
%                 newimage(i,j) = 0;
%             end
        end
    end

    figure(2);
    imshow(curdd);
    delta = abs(curdd-newimage);
%     max(max(delta))
[M,I] = max(delta);
[mm, ii] = max(M);
[mmm,iii] = max(delta(:,ii))
    max(max(delta))
        figure(3);
    imshow(delta);
    figure(4);
    imshow(cur_left);
    hold on;
    plot(ii,iii,'or');
    hold on;
    
    b = curdd>0;
    sum(sum(b))
    a = delta<.2& delta>0 ;
    sum(sum(a))

    
end