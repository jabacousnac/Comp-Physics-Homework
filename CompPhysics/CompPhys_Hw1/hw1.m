function hw1

%% QUESTION 1
% Initialize.
numdtpt = 20; %total number of points that will be plotted
x = single(0.1);
err_matrix = zeros(numdtpt,4); %this is an error matrix that we will populate
err_matrix = single(err_matrix);
%h = single(linspace(0,1,1000000)');
h = single(logspace(-6,0,numdtpt));
err_matrix(:,1) = log10(h);

%do it for the first method (forward)
y1_cos = (cos(x+h)-cos(x))./h;
%second method (central)
y2_cos = (cos(x+h)-cos(x-h))./(2*h);
%third method (extrapolated differences)
y3_cos = (-cos(x+2*h)+8*cos(x+h)-8*cos(x-h)+cos(x-2*h))./(12*h);

%{
%populate the error matrices
err_matrix(:,2) = log10(abs((y1_cos-(-sin(x)))./(-sin(x))));
err_matrix(:,3) = log10(abs((y2_cos-(-sin(x)))./(-sin(x))));
err_matrix(:,4) = log10(abs((y3_cos-(-sin(x)))./(-sin(x))));
%}

%same for exponential expressions
y1_exp = (exp(x+h)-exp(x))./h;
y2_exp = (exp(x+h)-exp(x-h))./(2*h);
y3_exp = (-exp(x+2*h)+8*exp(x+h)-8*exp(x-h)+exp(x-2*h))./(12*h);
%populate the error matrices
err_matrix(:,2) = log10(abs((y1_exp-(exp(x)))./(exp(x))));
err_matrix(:,3) = log10(abs((y2_exp-(exp(x)))./(exp(x))));
err_matrix(:,4) = log10(abs((y3_exp-(exp(x)))./(exp(x))));
err_matrix(:,1) = log10(h);

% Let the plotting begin...
figure(1)
plot(err_matrix(:,1), err_matrix(:,2), 'r-');
hold on;
plot(err_matrix(:,1), err_matrix(:,3), 'g-');
plot(err_matrix(:,1), err_matrix(:,4), 'b-');
xlabel('log_{10}|h|');
ylabel('log_{10}|Error in d/dx[e(x)]|'); %change label for cos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% QUESTION 2
% In this question, we integrate exp(-t) from 0 to 1
a = single(0);
b = single(1);

% First method -- Midpoint Method
numdtpt = 20;
N = logspace(0,7,numdtpt);
I_mid_matrix = single(zeros(length(N),3));
I_mid_matrix(:,1) = N; %first column will be number of bins
for i = 1:length(N)
    numbins = N(i);
    h = single((b-a)/numbins);
    t = a + h/2;%initial value of t
    sum  = single(0);
    while t<1
        sum = (sum + exp(-t));%evaluate function at each t
        t = t+h; %increment t by h
    end
    I_mid_matrix(i,2) = sum*(b-a)/(numbins);
end
I_true = single((exp(1)-1)/exp(1));
I_mid_matrix(:,3) = abs((I_mid_matrix(:,2)-I_true)/I_true);
    
% Second method -- Trapezoid Rule
N = logspace(0,7,numdtpt);
I_trap_matrix = single(zeros(length(N),3));
I_trap_matrix(:,1) = N;
for i = 1:length(N)
    numbins = N(i);
    h = single((b-a)/numbins);
    t = h; %initial value of t has to be at the right edge of the first bin
    sum = single(0);
    while t<1
        t_prev = t-h;
        sum = (sum + exp(-t_prev) + exp(-t));
        t = t+h; %increment t by h
    end
    I_trap_matrix(i,2) = 0.5*h*sum; 
end
I_trap_matrix(:,3) = abs((I_trap_matrix(:,2)-I_true)/I_true);

% Third method -- Simpson's Rule
N = logspace(0,7,numdtpt);
I_simp_matrix = single(zeros(length(N),3));
I_simp_matrix(:,1) = N;
for i = 1:length(N)
    numbins = N(i); %we will keep addding with this number of bins
    h = single((b-a)/numbins);
    sum1 = single(0);
    sum2 = single(0);
    for n1 = 1:2:numbins %note how we go over every 2nd term
        sum1 = sum1 + exp(-(a+(n1*h)));
    end
    for n2 = 2:2:numbins-2
        sum2 = sum2 + exp(-(a+(n2*h)));
    end
    %Now, we can add the terms
    I_simp_matrix(i,2) = (1/3)*h*(exp(-1*a) + exp(-1*b) + 4*sum1 + 2*sum2);
end
I_simp_matrix(:,3) = abs((I_simp_matrix(:,2)-I_true)/I_true);

figure(2)
plot(log10(I_mid_matrix(:,1)), log10(I_mid_matrix(:,3)), 'r-');
hold on;
plot(log10(I_trap_matrix(:,1)), log10(I_trap_matrix(:,3)), 'g-');
plot(log10(I_simp_matrix(:,1)), log10(I_simp_matrix(:,3)), 'b-');
xlabel('log_{10}|Number of bins|');
ylabel('log_{10}|\epsilon|');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% QUESTION 3
% We have a bunch of values k values and a bunch of P(k) values. 
% It is given as data. We need to interpolate so we have 'continuum' of
% values of P(k) for any k. 

% STEP 1: save and load k and P_k
k_query = logspace(-4,1,100);
load('k.mat','k');
load('P_k.mat','P_k');

% STEP 2: Interpolate for P(k)
s = spline(k,P_k,k_query); %MATLAB's built-in cubic spline interpolation.
    %Thank you, MATLAB.
figure(3)
plot(log10(k_query),log10(s),'k-');
ylabel('log_{10}|P_{k}|');
xlabel('log_{10}|k|');

% STEP 3: Calculate the integrand
N = 100; %this will be the number of r values between [5,120]
r_matrix = linspace(0,120,N);
xi_matrix = zeros(1,length(r_matrix));
numbins = 200000;
a = 0.0001;
b = 1000;
h = (b-a)/numbins;
my_spline = spline(k,P_k,(a+(h/2):h:b));


for index = 1:length(r_matrix)
    %for some r:
    r = r_matrix(index);
    x = a+(h/2);
    I = 0;
    for i = 1:b %from 1 to some large value
        p = my_spline(i);
        x = x+h; %increment my k by h
        I = I + (x^2*p*sin(x*r)/(x*r)); %evaluate the integrand
    end
    xi_matrix(index) = (I*(b-a)/numbins)/(2*pi^2);
end

figure(4)
plot(r_matrix, (r_matrix).^2.*xi_matrix,'k-')
xlabel('r');
ylabel('r^{2}\xi(r)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end