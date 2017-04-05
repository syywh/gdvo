function [ ans ] = g_robust_kernal( t, K )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

idx = abs(t)>K;
ans = (1- (t./K).^2).^2;
ans(idx) = 0;
end

