function T = estimate_ds_dense_gradient_Joint(cur,ref_gra_x,ref_gra_y, kf_gra_x, kf_gra_y,cur_gra_x,cur_gra_y,refdepth, kfdepth,curdepth,T, Tkf_r)
dscale = 0.0002./1.032;
pcx = 318.6; pcy = 255.3; pfx = 517.3; pfy = 516.5;                
% pcx = 335.963; pcy = 190.625; pfx = 351.498; pfy = 351.498; 
s0 = 8;
near = 0.5;     %保留的点
far = 4;
iteration_s = [30 30 ];
refdepth = dscale*double(refdepth);
curdepth = dscale*double(curdepth);
kfdepth = dscale*double(kfdepth);
% ref_gra_x = im2double(ref_gra_x);
% ref_gra_y = im2double(ref_gra_y);
% cur_gra_x = im2double(cur_gra_x);
% cur_gra_y = im2double(cur_gra_y);
weight = 1;

%% Process
s = s0;
for n = 1:length(iteration_s)
    % pyramid
    cx = pcx/s; cy = pcy/s; fx = pfx/s; fy = pfy/s;
    K = [fx 0 cx;0 fy cy;0 0 1];
    refdepthp = imresize(refdepth,1/s,'nearest');
    curdepthp = imresize(curdepth,1/s,'nearest');
    kfdepthp = imresize(kfdepth, 1/s, 'nearest');
    curp = imresize(cur,1/s);
    refp_gra_x = imresize(ref_gra_x,1/s);
    refp_gra_y = imresize(ref_gra_y, 1/s);
    curp_gra_x = imresize(cur_gra_x,1/s);
    curp_gra_y = imresize(cur_gra_y, 1/s);
    kfp_gra_x = imresize(kf_gra_x, 1/s);
    kfp_gra_y = imresize(kf_gra_y, 1/s);
    [uref,vref] = meshgrid(1:size(refdepthp,2),1:size(refdepthp,1));    % 填充1到col，1到row
    [ukf, vkf] = meshgrid(1:size(kfdepthp,2),1:size(kfdepthp,1));
%     [refgxx,refgxy] = imgradientxy(refp_gra_x);
%     [refgyx,refgyy] = imgradientxy(refp_gra_y);
    gref = sqrt(refp_gra_x.^2+refp_gra_y.^2);
    gref = gref(:).';
    
    gkf = sqrt(kfp_gra_x.^2 + kfp_gra_y.^2);
    gkf = gkf(:).';
    
    % feature selection
    uref = uref(:).'-1;%按着列，再行横着来
    vref = vref(:).'-1;
    ukf = ukf(:).' -1;
    vkf =  vkf(:).'-1;
    
    zref = refdepthp(:).';
    zkf = kfdepthp(:).';
    
    Pref = K\[zref.*uref;zref.*vref;zref];
    refvalid = Pref(3,:)>near & Pref(3,:)<far & gref>0.5/s;
    
    Pkf = K\[zkf.*ukf;zkf.*vkf;zkf];
   kfvalid = Pkf(3,:)>near & Pkf(3,:)<far & gkf>0.5/s;
    
    % jacobian    
    maxiter = iteration_s(n);%每次迭代的次数
    err = zeros(1,maxiter);
    RK = 4.6851;
% RK = 50;
    for m = 1 : maxiter
        %% for reference image
        render = zeros(size(curp_gra_x));
        Rotate = T(1:3,1:3);
        Translate = T(1:3,4);
        ENewT = Pref - repmat(Translate,1,size(Pref,2));
        ENewRT = Rotate.'*ENewT;
        uvNew = K*ENewRT;
        uNew = round(uvNew(1,:)./uvNew(3,:))+1;
        vNew = round(uvNew(2,:)./uvNew(3,:))+1;
        uu = uNew;  vv = vNew;  uv = uvNew;
        idx = (uNew(:)-1)*size(refdepthp,1)+vNew(:);
        valid = uvNew(3,:)>0.5 & uvNew(3,:)<4 & uNew <= size(curdepthp,2) & uNew > 0 ...
            & vNew <= size(curdepthp,1) & vNew > 0 & refvalid;
        validvalid = valid;
%         render(idx(valid)) = refp(valid);
        Vsum = 0;
        V = [];    J = [];
        Vdsum = 0;
        H = zeros(6,6);
        b = zeros(6,1);
        [gdxx,gdxy] = imgradientxy(curp_gra_x);  
        [gdyx,gdyy] = imgradientxy(curp_gra_y);  
%         [gdxd,gdyd] = imgradientxy(curdepthp);
         

        for i =1 : size(Pref,2)
            if ~valid(i)
                continue;
            end             
            refPatch_x = refp_gra_x(vref(i)+1,uref(i)+1);
            refPatch_y = refp_gra_y(vref(i)+1,uref(i)+1);
            curPatch_x = curp_gra_x(vNew(i),uNew(i));
            curPatch_y = curp_gra_y(vNew(i),uNew(i));
            gdxxPatch = gdxx(vNew(i),uNew(i));
            gdxyPatch = gdxy(vNew(i),uNew(i));
            gdyxPatch = gdyx(vNew(i),uNew(i));
            gdyyPatch = gdyy(vNew(i),uNew(i));
           
            Vi_x = refPatch_x - curPatch_x;     
            Vi_y = refPatch_y - curPatch_y;
%             Vx = [Vx; Vi_x];
%             Vy = [Vy; Vi_y];
            V= [V;[Vi_x; Vi_y]];
            Vsum = Vsum+(Vi_x).^2+(Vi_y).^2;

            EiNewT = ENewT(:,i);        
            EiNewRT = ENewRT(:,i);
            firstDiffx = [gdxxPatch gdxyPatch];
            firstDiffy = [gdyxPatch gdyyPatch];
%             firstDiffd = [gdxPatchd gdyPatchd];
            secondDiff = [ fx / EiNewRT(3), 0, -fx * EiNewRT(1) * EiNewRT(3).^(-2);
                                  0, fy / EiNewRT(3), -fy * EiNewRT(2) * EiNewRT(3).^(-2) ];
            EiNewTx = [ 0, -EiNewT(3), EiNewT(2);
                               EiNewT(3), 0, -EiNewT(1);
                               -EiNewT(2), EiNewT(1), 0];
            thirdDiff = Rotate' * [ -eye(3), EiNewTx];
            Jix = firstDiffx * secondDiff * thirdDiff;
            Jiy = firstDiffy * secondDiff * thirdDiff;

            Ji = [Jix; Jiy];
            J = [J ;Ji];

        end

     %% for keyframe image
       Tkf_c = Tkf_r*T;
        Rotate = Tkf_c(1:3,1:3);
        Translate = Tkf_c(1:3,4);
        ENewT = Pkf - repmat(Translate,1,size(Pkf,2));
        ENewRT = Rotate.'*ENewT;
        uvNew = K*ENewRT;
        uNew = round(uvNew(1,:)./uvNew(3,:))+1;
        vNew = round(uvNew(2,:)./uvNew(3,:))+1;
        uu = [uu uNew];     vv = [vv vNew];     uv = [uv uvNew];
        valid = uvNew(3,:)>0.5 & uvNew(3,:)<4 & uNew <= size(curdepthp,2) & uNew > 0 ...
            & vNew <= size(curdepthp,1) & vNew > 0 & kfvalid;
        validvalid = [validvalid valid];
%         [gdxx,gdxy] = imgradientxy(curp_gra_x);  
%         [gdyx,gdyy] = imgradientxy(curp_gra_y);  
%         [gdxd,gdyd] = imgradientxy(curdepthp);
         

        for i =1 : size(Pkf,2)
            if ~valid(i)
                continue;
            end             
            kfPatch_x = kfp_gra_x(vkf(i)+1,ukf(i)+1);
            kfPatch_y = kfp_gra_y(vkf(i)+1,ukf(i)+1);
            curPatch_x = curp_gra_x(vNew(i),uNew(i));
            curPatch_y = curp_gra_y(vNew(i),uNew(i));
            gdxxPatch = gdxx(vNew(i),uNew(i));
            gdxyPatch = gdxy(vNew(i),uNew(i));
            gdyxPatch = gdyx(vNew(i),uNew(i));
            gdyyPatch = gdyy(vNew(i),uNew(i));
           
            Vi_x = kfPatch_x - curPatch_x;     
            Vi_y = kfPatch_y - curPatch_y;
%             Vx = [Vx; Vi_x];
%             Vy = [Vy; Vi_y];
            V= [V;[Vi_x; Vi_y]];
            Vsum = Vsum+(Vi_x).^2+(Vi_y).^2;

            EiNewT = ENewT(:,i);        
            EiNewRT = ENewRT(:,i);
            firstDiffx = [gdxxPatch gdxyPatch];
            firstDiffy = [gdyxPatch gdyyPatch];
%             firstDiffd = [gdxPatchd gdyPatchd];
            secondDiff = [ fx / EiNewRT(3), 0, -fx * EiNewRT(1) * EiNewRT(3).^(-2);
                                  0, fy / EiNewRT(3), -fy * EiNewRT(2) * EiNewRT(3).^(-2) ];
            EiNewTx = [ 0, -EiNewT(3), EiNewT(2);
                               EiNewT(3), 0, -EiNewT(1);
                               -EiNewT(2), EiNewT(1), 0];
            thirdDiff = Rotate' * [ -eye(3), EiNewTx];
            Jix = firstDiffx * secondDiff * thirdDiff;
            Jiy = firstDiffy * secondDiff * thirdDiff;

            Ji = [Jix; Jiy];
            J = [J ;Ji];

        end

        %???是否一起算？
        median_v = median(V);
        V_n = V - median_v;
        V_n = V_n./std(V_n);%标准差

        W = diag( g_robust_kernal(V_n, RK).^2 );
        S = (J' * W *J )\( J'*W*V);
%         S = (H\b).';
        t = S(1 : 3); w = S(4 : 6); 
        wx = [ 0, -w(3), w(2); 
                  w(3), 0, -w(1);
                  -w(2), w(1), 0 ];
        deltaRT = expm([wx t;0 0 0 0]);
        q = rotm2quat(deltaRT(1:3,1:3));
        q = q./norm(q);
        deltaRT(1:3,1:3) = quat2rotm(q);
        T = T * deltaRT;
        
        err(m) = Vsum;%+Vdsum;        
        figure(2)
        plot(err);
        figure(3)
        cm = colormap(hsv);
        imshow(curp);hold on;
%         scatter(uNew(:,valid),vNew(:,valid),5,cm(round(uvNew(3,valid)./4*size(cm,1)),:),'fill');hold off;
        scatter(uu(:,validvalid),vv(:,validvalid),5,cm(round(uv(3,validvalid)./4*size(cm,1)),:),'fill');hold off;
        drawnow;
    end
    s = s /2;
end
% plot(err);
end