import numpy as np
import math
from math import exp
from numpy import empty,zeros
from pylab import plot, imshow, show, xlabel, ylabel, colorbar
import matplotlib.pyplot as plt
import random as ran
from random import random, randrange
#from random import random, shuffle, randrange
import os

M = 50; #size of array
Tmax = 100;
Tmin = 10**-4;
k = 1;
tau = 10**6; #experiment with different tau's
#start by choosing a random site
def pick_site(L): #L=[i,j]
    i = L[0];
    j = L[1];
    myL = [[i+1,j], [i-1,j], [i,j+1], [i,j-1]];
    return ran.choice(myL) #pick a neighboring site

def sim_annealing():
    totals_list = list(); #contains the total energy
    lattice = empty((M,M),np.float64);
    T = Tmax;
    t = 0; #time
    q = 0; #a counter
    while T>Tmin:
        t += 1;
        T = Tmax*exp(-t/tau)#cool time
        #now, pick 2 adjacent sites
        i1 = randrange(1,M-1);
        j1 = randrange(1,M-1);
        L = pick_site([i1,j1]);
        i2 = L[0];
        j2 = L[1];
        totals_list.append(sum(sum(lattice)));
        if lattice[i1][j1] == 0 and lattice[i2][j2] == 0: #then, add a dimer
            lattice[i1][j1]-=1; #dimers carry 'negative' energy.
            lattice[i2][j2]-=1;
        elif lattice[i1][j1] == -1 and lattice[i2][j2] == -1: #dimer already present. Do probability thing. 
            lattice_next_step = lattice.copy();
            lattice_next_step[i1][j1]*=0;
            lattice_next_step[i2][j2]*=0;
            prob = exp(-1/(k*T));
            if random()<prob:
                lattice = lattice_next_step; #then we remove the dimers
            else:
                lattice = lattice;
        else:
            pass
        print(T);
        if q==1 or q==10 or q==100 or q==1000 or q==10000 or q==100000 or q==1000000:
            plt.imshow(lattice, cmap='gray', interpolation='nearest');
            plt.xlabel("x");
            plt.ylabel("y");
            plt.colorbar();
            title = "tau=" + str(tau) + " , " + "T = " + str(T);
            plt.title(title);
            filename = "dimers_tau=" + str(tau) + "_T=" + str(T) + ".png";
            plt.savefig(filename);
            plt.close();
        q+=1;
    plt.imshow(lattice, cmap='gray', interpolation='nearest');
    title = "tau=" + str(tau) + " , " + "T = " + str(T);
    plt.title(title);
    filename = "dimers_tau=" + str(tau) + "_T=" + str(T) + ".png";
    plt.xlabel("x");
    plt.ylabel("y");
    plt.colorbar();
    plt.savefig(filename);
    plt.close();
    #another plot
    plot(totals_list);
    ylabel("Energy");
    xlabel("time");
    show();

    return [tau, totals_list[-1]]
                
if __name__ == "__main__":
    sim_annealing()
        
    
        
        
    
