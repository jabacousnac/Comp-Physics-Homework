import numpy as np
from math import exp
from numpy import empty,zeros
from pylab import plot, imshow, show, xlabel, ylabel, colorbar
import matplotlib.pyplot as plt
import random
from random import random, shuffle, randrange
import os

M = 20; #size of array
J = -1; #can remove this and include a term later if J is not simply a constant
kT = 1;

def get_energy(lattice): #function to calculate the total energy
    energy = 0;
    for row in range(0,M):
        for col in range(0,M):
            if (row == M-1 and col != M-1): #this is the last row
                energy += J*lattice[row][col]*lattice[row][col+1]; #right
            elif (col == M-1 and row != M-1): #this is the last column
                energy += J*lattice[row][col]*lattice[col][row+1]; #bottom
            else:
                if (row == M-1 and col == M-1):
                    energy += 0;
                else:
                    energy += J*lattice[row][col]*(lattice[row+1][col] + lattice[row][col+1]); #bottom and right
    return energy

def ising(N): #N is the number of steps for which we make the system evolve
    '''
    #If you want exactly half and half spins:
    energy_list = list();
    lattice = list();
    spins = int(M**2/2)*([1] + [-1]);
    shuffled = spins.copy();
    shuffle(shuffled);
    for i in range(len(shuffled)):
        if  i%M == 0:
            lattice.append(shuffled[i:i+M])
    #print(lattice)
    '''

    #1. initialize the array with roughly half and half spins
    energy_list = list(); #total energy list
    mag_list = list(); #magnetization list
    lattice = empty((M,M),int);
    for i in range(M):
        for j in range(M):
            if random()<0.5:
                lattice[i,j]=1;
            else:
                lattice[i,j]=-1;
    
    #2. choose random spin to flip and accept or reject
    for i in range(N):
        #save figure for comparison
        if i%100000 == 0 or i == 10 or i == 100 or i == 1000 or i == 10000:
            print(i)


            #create folder and save figure
            myString = "Ising_kT=" + str(kT);
            directory = './'+myString+'/';
            filename = "ising_step_kT=" + str(kT) + "_step_number_"  + str(i);
            file_path = os.path.join(directory, filename)
            if not os.path.isdir(directory):
                os.mkdir(directory)
            #create plot
            plt.imshow(lattice, cmap='gray', interpolation='nearest');
            plt.xlabel("x");
            plt.ylabel("y");
            save_results_to = file_path;
            plt.savefig(save_results_to, dpi = 300);
            
        #flip
        x = randrange(M);
        y = randrange(M);
        E_i = get_energy(lattice);
        energy_list.append(E_i);
        mag_list.append(sum(sum(lattice)));
        next_step_lattice = lattice.copy();
        next_step_lattice[x][y]*=-1; #this is a spin flip
        E_j = get_energy(next_step_lattice);
        if E_j < E_i:
            lattice = next_step_lattice;
        else:
            prob = exp(-(E_j - E_i)/kT);
            if random()<prob:
                lattice = next_step_lattice;
    plt.close();
    plot(range(N),mag_list);
    xlabel("# Time Steps");
    ylabel("Magnetization");
    show();
    return 
