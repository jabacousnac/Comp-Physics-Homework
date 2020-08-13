import numpy as np
import matplotlib.pyplot as plt
from numpy import zeros
from cmath import exp,pi

def get_data():
    datalist = [];
    with open("sunspots.txt", 'r') as text: #open sunspots file
        for line in text:
            data = line.strip().split(',');
            datalist.append(float(data[1]));#append second element because that is the sunspot number
            #So, now we have datalist, a list of that contains all the sunspot numbers without commas and all the unnecessary garbage
        return datalist

def dft(datalist):
    coeffsq = [];
    N = len(datalist);
    c = zeros(N//2+1,complex);#normally, we only need to calculate it for half of our sample if sample is real
    for k in range(N//2+1):
        for n in range(N):
            c_k = datalist[n]*exp(-2j*pi*k*n/N);
            c[k]+=c_k; #our fourier coefficients
    for coeff in c:
        coeffsq.append(coeff.real**2+coeff.imag**2)
#plot stuff:
    fig = plt.figure();
    ax1 = fig.add_subplot(111);#in case I want to divide my plot
    plt.xlabel("$k$");
    plt.ylabel("$|C_{k}|^2$");
    line = ax1.plot(list(range(1,1+N//2+1)),coeffsq, lw=2);
    plt.show()
    return 


