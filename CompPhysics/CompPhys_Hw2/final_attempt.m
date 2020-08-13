function out = final_attempt(alpha, phistar,mstar);
load('N_mgal.mat');
load('Mgal.mat');
out = 0;
for mgal = logspace(9.5,11.5,1000)
    %calculate n for a bunch of mgal
    n=@(mgal) log(10)*phistar*((mgal./mstar).^(alpha+1)).*exp(-mgal./mstar);%theory
    n_t = n(mgal);%theory
    n_o = spline(10.^Mgal,N_mgal,mgal);%obs
    out = out + ((n_o-n_t)^2)/abs((n_t)); 
end
%essentially, we are outputting the chi-squared value
end