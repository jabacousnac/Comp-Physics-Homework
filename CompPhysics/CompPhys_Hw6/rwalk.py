import numpy as np
import math
from math import exp, pi, cos, sin, log
from numpy import empty,zeros, log10
from pylab import plot, imshow, show, xlabel, ylabel, colorbar
import matplotlib.pyplot as plt
from cmath import exp,pi
import random

def ANSI():
    X_n = random.randint(1,10000); #random seed
    a = 1103515245;
    c = 12345;
    m = 2**31;
    return ((a*X_n+c)%m)/(m+1)

def gaussian(): #output n gaussian numbers
    n = 10000;
    counter = 0;
    sigma = 1;
    random_list = list();
    while counter<n:
        counter+=1;
        r = (-2*sigma*sigma*log(1-ANSI()))**0.5;
        theta = 2*pi*(ANSI());
        x = r*cos(theta);
        #y = r*sin(theta);
        random_list.append(x);
    num_bins = 20;
    n, bins, patches = plt.hist(random_list, num_bins, facecolor='blue', alpha=0.5)
    xlabel("$X_{n}$");
    ylabel("Counts");
    plt.show()
    return random_list

def dft(): #code from homework #3
    #datalist = gaussian();
    datalist = rwalk();
    coeffsq = [];
    N = len(datalist);
    c = zeros(N//2+1,complex);#normally, we only need to calculate it for half of our sample if sample is real
    for k in range(N//2+1):
        for n in range(N):
            c_k = datalist[n]*exp(-2j*pi*k*n/N);
            c[k]+=c_k; #our fourier coefficients
    for coeff in c:
        coeffsq.append(log10(coeff.real**2+coeff.imag**2))
#plot stuff:
    fig = plt.figure();
    ax1 = fig.add_subplot(111);#in case I want to divide my plot
    plt.xlabel("$log |k|$");
    plt.ylabel("$log (|C_{k}|^2)$");
    line = ax1.plot(list(log10(range(1,1+N//2+1))),coeffsq, lw=2);
    plt.show()
    return 

def rwalk():
    n = 10000;
    datalist = gaussian();
    x = np.zeros(n) 
    for i in range(1,n):
        x[i] = x[i-1]+datalist[i-1];
    plot(x);
    plt.title("Random Walk ($n = " + str(n) + "$ steps)") 
    show() 
    return x

        
    
        
        
    
