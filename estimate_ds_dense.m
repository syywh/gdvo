function T = estimate_ds_dense(ref,cur,refdepth,curdepth,T)
dscale = 0.0002./1.032;
pcx = 318.6; pcy = 255.3; pfx = 517.3; pfy = 516.5;                
% pcx = 335.963; pcy = 190.625; pfx = 351.498; pfy = 351.498; 
s0 = 8;
near = 0.5;
far = 4;
ocluthres = 0.5;
iteration_s = [30 30];
refdepth = dscale*double(refdepth);
curdepth = dscale*double(curdepth);
ref = im2double(ref);
cur = im2double(cur);
weight = 1;

%% Process
s = s0;
for n = 1:length(iteration_s)
    % pyramid
    cx = pcx/s; cy = pcy/s; fx = pfx/s; fy = pfy/s;
    K = [fx 0 cx;0 fy cy;0 0 1];
    refdepthp = imresize(refdepth,1/s,'nearest');
    curdepthp = imresize(curdepth,1/s,'nearest');
    refp = imresize(ref,1/s);
    curp = imresize(cur,1/s);
    [u,v] = meshgrid(1:size(refdepthp,2),1:size(refdepthp,1));
    [refgx,refgy] = imgradientxy(refp);
    g = sqrt(refgx.^2+refgy.^2);
    g = g(:).';
    
    % feature selection
    u = u(:).'-1;
    v = v(:).'-1;
    z = refdepthp(:).';
    P = K\[z.*u;z.*v;z];
    refvalid = P(3,:)>near & P(3,:)<far & g>0.2/s;
    
    % jacobian    
    maxiter = iteration_s(n);
    err = zeros(1,maxiter);
    for m = 1 : maxiter
        render = zeros(size(curp));
        Rotate = T(1:3,1:3);
        Translate = T(1:3,4);
        ENewT = P - repmat(Translate,1,size(P,2));
        ENewRT = Rotate.'*ENewT;
        uvNew = K*ENewRT;
        uNew = round(uvNew(1,:)./uvNew(3,:))+1;
        vNew = round(uvNew(2,:)./uvNew(3,:))+1;
        idx = (uNew(:)-1)*size(refdepthp,1)+vNew(:);
        valid = uvNew(3,:)>0.5 & uvNew(3,:)<4 & uNew <= size(curdepthp,2) & uNew > 0 ...
            & vNew <= size(curdepthp,1) & vNew > 0 & refvalid;
%         render(idx(valid)) = refp(valid);
        Vsum = 0;
        Vdsum = 0;
        H = zeros(6,6);
        b = zeros(6,1);
        [gdx,gdy] = imgradientxy(curp);  
%         [gdxd,gdyd] = imgradientxy(curdepthp);
         
        for i =1 : size(P,2)
            if ~valid(i)
                continue;
            end             
            refPatch = refp(v(i)+1,u(i)+1);
            curPatch = curp(vNew(i),uNew(i));
            gdxPatch = gdx(vNew(i),uNew(i));
            gdyPatch = gdy(vNew(i),uNew(i));
            
            refPatchd = ENewRT(3,i);
            curPatchd = curdepthp(vNew(i),uNew(i));
%             gdxPatchd = gdxd(vNew(i),uNew(i));
%             gdyPatchd = gdyd(vNew(i),uNew(i));
           
            Vi = refPatch - curPatch;
%             Vdi = sqrt(weight)*(refPatchd - curPatchd);
            Vsum = Vsum+(Vi).^2;
%             Vdsum = Vdsum+(Vdi).^2;
%             if abs(Vdi)>ocluthres
%                 continue;
%             end
            EiNewT = ENewT(:,i);        
            EiNewRT = ENewRT(:,i);
            firstDiff = [gdxPatch gdyPatch];
%             firstDiffd = [gdxPatchd gdyPatchd];
            secondDiff = [ fx / EiNewRT(3), 0, -fx * EiNewRT(1) * EiNewRT(3).^(-2);
                                  0, fy / EiNewRT(3), -fy * EiNewRT(2) * EiNewRT(3).^(-2) ];
            EiNewTx = [ 0, -EiNewT(3), EiNewT(2);
                               EiNewT(3), 0, -EiNewT(1);
                               -EiNewT(2), EiNewT(1), 0];
            thirdDiff = Rotate' * [ -eye(3), EiNewTx];
            Ji = firstDiff * secondDiff * thirdDiff;
%             Jdi = sqrt(weight)*firstDiffd * secondDiff * thirdDiff;
            H = H + Ji.'*exp(-(maxiter-m)/maxiter*norm(Vi))*Ji;% + Jdi.'*exp(-(maxiter-m)/maxiter*norm(Vdi))*Jdi;
            b = b + Ji.'*exp(-(maxiter-m)/maxiter*norm(Vi))*Vi;% + Jdi.'*exp(-(maxiter-m)/maxiter*norm(Vdi))*Vdi;
        end
        S = (H\b).';
        t = S(1 : 3)'; w = S(4 : 6)'; 
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
        scatter(uNew(:,valid),vNew(:,valid),5,cm(round(uvNew(3,valid)./4*size(cm,1)),:),'fill');hold off;
        drawnow;
    end
    s = s /2;
end
% plot(err);
end