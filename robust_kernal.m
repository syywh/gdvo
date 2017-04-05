function ans = robust_kernal(t, K)
    
if(abs(t) < K)
    ans = K*K*(1-(1- (t/K)^2)^3 )/6;
else
    ans = K*K/6;
end

end