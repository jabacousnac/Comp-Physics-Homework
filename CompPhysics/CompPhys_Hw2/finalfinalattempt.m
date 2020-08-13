function finalfinal_attempt
%look at finalattempt.m, first.
%%
% load observed data:
load('N_mgal.mat');
load('Mgal.mat');
load('N_mgalerr.mat');

[alpha, phistar, mstar] = deal(-0.9, 10^-2.5, 10^11.5); %initial values
dx = 0.001; %increment

%define three matrices we will update to keep track of indices:
alpha_matrix = [alpha; alpha+dx];
phistar_matrix = [phistar; phistar+dx];
mstar_matrix = [mstar; mstar+dx];
gamma = 0.001; %learning rate
tol = 1; %tolerance
i = 1; %index
while abs(tol) > 0.0001
%update alpha
    i = i+1;
    alpha = alpha - (gamma*(final_attempt(alpha_matrix(i),phistar_matrix(i-1),mstar_matrix(i-1))-...
        final_attempt(alpha_matrix(i-1),phistar_matrix(i-1),mstar_matrix(i-1)))/...
        (alpha_matrix(i)-alpha_matrix(i-1)));
%phistar
    phistar = phistar - (gamma*(final_attempt(alpha_matrix(i-1),phistar_matrix(i),mstar_matrix(i-1))-...
        final_attempt(alpha_matrix(i-1),phistar_matrix(i-1),mstar_matrix(i-1)))/...
        (phistar_matrix(i)-phistar_matrix(i-1)));
%mstar
    mstar = mstar - (gamma*(final_attempt(alpha_matrix(i-1),phistar_matrix(i-1),mstar_matrix(i))-...
        final_attempt(alpha_matrix(i-1),phistar_matrix(i-1),mstar_matrix(i-1)))/...
        (mstar_matrix(i)-mstar_matrix(i-1)));
%add variables to matrix:
    alpha_matrix = [alpha_matrix;alpha];
    phistar_matrix = [phistar_matrix;phistar];
    mstar_matrix = [mstar_matrix;mstar];
%display chi-squared value and the other parameters
    chisq = final_attempt(alpha, phistar, mstar);
    disp(chisq);%DOES NOT MINIMIZE :/
%update tolerance
    tol = alpha_matrix(i)-alpha_matrix(i-1);
end
end