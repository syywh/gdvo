% frame_num = size(trajframeT,2)/4;
% trajframeT = reshape(trajframeT, 4,4 , frame_num);
% T_real = load('../00.txt');

figure;
for i = 1:frame_num
    plot3(trajframeT(1,4,i), trajframeT(2,4,i), trajframeT(3,4,i),'b*');    hold on;
    plot3(T_real(i,4), T_real(i,8), T_real(i,12),'r*');    hold on;
end
 axis equal
 
 figure;
 for i = 1:frame_num
    plot3(trajframeT(1,4,i)-T_real(i,4), trajframeT(2,4,i)-T_real(i,8), trajframeT(3,4,i)-T_real(i,12),'b*');    hold on;
   
end
 axis equal
 grid on