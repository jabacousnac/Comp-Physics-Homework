function hw2_gdescent
%% PART ONE: IMPLEMENT ON SIMPLER FUNCTION
%work on function f(x,y) = (x-2)^2 + (y-2)^2
f=@(x,y) (x-2)^2+(y-2)^2;%my function
gamma = 0.01;
[x, y, dx, dy, tol1, tol2] = deal(6,10, 0.01, 0.01, 1, 1);
xmatrix = [x; x+dx];
ymatrix = [y; y+dy];
fmatrix = [f(x,y); f(x+dx,y+dy)];
i = 1; %some counter
while (abs(tol1))&&(abs(tol2)) > 0.00001
    i = i+1;
    x = x - (gamma*(f(xmatrix(i),ymatrix(i-1))-f(xmatrix(i-1),ymatrix(i-1)))/...
        (xmatrix(i)-xmatrix(i-1)));
    y = y - (gamma*(f(xmatrix(i-1),ymatrix(i))-f(xmatrix(i-1),ymatrix(i-1)))/...
        (ymatrix(i)-ymatrix(i-1)));
    xmatrix = [xmatrix;x];
    ymatrix = [ymatrix;y];
    tol1 = xmatrix(i)-xmatrix(i-1);
    tol2 = ymatrix(i)-ymatrix(i-1);
    fmatrix = [fmatrix; f(x,y)];
end

figure(1)
plot(ymatrix(1:14:end),'r*');
hold on;
plot(xmatrix(1:15:end), 'b*');
xlim([2, inf]);
xlabel('# iterations');
ylabel('parameter');
figure(2)
plot((1:length(fmatrix)), fmatrix);
xlabel('# iterations');
ylabel('f(x,y)');

%{ 
%% PART TWO: IMPLEMENT ON SCHECHTER FUNCTION (INITIAL IMPLEMENTATION --
SCRAPPED!)
[alpha, phistar, mstar] = deal(-0.9, 10^-2.5, 10^10.5);
n=@(mgal,alpha,phistar,mstar) log(10)*phistar*((mgal./mstar).^(alpha+1)).*exp(-mgal./mstar);%Schechter function
%Mgal is a saved variable (as Mgal.m), and likewise for N_mgal.m:
load('N_mgal.mat');
load('Mgal.mat');
load('N_mgalerr.mat');

%first, let us create a cubic spline interpolation from observed data
n_obs = spline(10.^Mgal,N_mgal,logspace(9.5,11.8,20));%change 20 to something
    %to match the value of h that will be hereafter chosen
%figure(2);
%plot(linspace(9.5,11.8,20), log10(n_obs),'r*') %check if the spline function works!
%hold on;

%when calculating chi-sq., we will subtract from our spline. Okay?
%Now, let us evaluate the Schechter function with small step sizes
max_it = 1000;
%n_theory = n(logspace(9.5,11.8,max_it));
%plot(linspace(9.5,11.8,max_it), log(n_theory), 'b-');

dx = 0.001; %increment all quantities by that (delta x)
%some matrices that we will continuously update upon each evaluation, later
%on.
alpha_matrix = [alpha; alpha+dx];
phistar_matrix = [phistar; phistar+dx];
mstar_matrix = [mstar; mstar+dx];
chisq=@(mgal, alpha, phistar, mstar)... %Begin function:
    (((log(10)*phistar*((mgal./mstar).^(alpha+1)).*exp(-mgal./mstar))...%model
    - spline(10.^Mgal,N_mgal,mgal))^2) ...%observed
    /(spline(10.^Mgal,N_mgalerr,mgal))^2;
%now, let us minimize the chi-squared fit, using gradient descent
gamma = 0.001;
tol1 = 1;
mgal = [10^9];
i = 1;
sum = 0;
sum_matrix = [];
while abs(tol1) > 0.00001
    i = i+1;
    chisq=@(mgal, alpha, phistar, mstar)... %Begin function:
    (((log(10)*phistar*((mgal./mstar).^(alpha+1)).*exp(-mgal./mstar))...%model
    - spline(10.^Mgal,N_mgal,mgal))^2) ...%observed
    /(spline(10.^Mgal,N_mgalerr,mgal))^2;
%alpha:
    alpha = alpha - (gamma*(chisq(mgal, alpha_matrix(i),phistar_matrix(i-1),mstar_matrix(i-1))-...
        chisq(mgal, alpha_matrix(i-1),phistar_matrix(i-1),mstar_matrix(i-1)))/...
        (alpha_matrix(i)-alpha_matrix(i-1)));
%phistar
    phistar = phistar - (gamma*(chisq(mgal, alpha_matrix(i-1),phistar_matrix(i),mstar_matrix(i-1))-...
        chisq(mgal, alpha_matrix(i-1),phistar_matrix(i-1),mstar_matrix(i-1)))/...
        (phistar_matrix(i)-phistar_matrix(i-1)));
%mstar
    mstar = mstar - (gamma*(chisq(mgal, alpha_matrix(i-1),phistar_matrix(i-1),mstar_matrix(i))-...
        chisq(mgal, alpha_matrix(i-1),phistar_matrix(i-1),mstar_matrix(i-1)))/...
        (mstar_matrix(i)-mstar_matrix(i-1)));
%add variables to matrix:
    alpha_matrix = [alpha_matrix;alpha];
    phistar_matrix = [phistar_matrix;phistar];
    mstar_matrix = [mstar_matrix;mstar];
    
    %Now calculate the chi-squared:
    sum = 0;
    for mgal = logspace(9.5,11.5,1000)
        sum = sum+(((log(10)*phistar*((mgal./mstar).^(alpha+1)).*exp(-mgal./mstar))...%model
    - spline(10.^Mgal,N_mgal,mgal))^2) ...%observed
    /(spline(10.^Mgal,N_mgalerr,mgal))^2;
    end
    %disp(sum); %this is the chi-square value that should be decreasing
    
    tol1 = alpha_matrix(i)-alpha_matrix(i-1);
end
%}

end