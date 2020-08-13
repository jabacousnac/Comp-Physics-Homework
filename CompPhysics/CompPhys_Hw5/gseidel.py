import numpy as np
import sys
from numpy import zeros, empty, array, arange,max, polyfit
import matplotlib.pyplot as plt
from pylab import plot, imshow, show, xlabel, ylabel, colorbar, arrow
np.set_printoptions(threshold=sys.maxsize) #this line causes python to nottruncate our array when printing

'''
PART 2(a)
'''
def cmap(): #generate the charge density field colormap. 
    dx = 1; #dx is the length of the unit cell
    #get (x,y) values
    x_list = list();
    y_list = list();
    with open("particles.txt", "r") as my_file:
        for line in my_file:
            string = line.split();
            x_list.append(float(string[0]));
            y_list.append(float(string[1]));
    #now we have two lists, x_list and y_list.
    #create a 2D array
    M = 100;
    grid = zeros([M+1,M+1],float);
    #now, we assign a charge density value to each point in the 100x100 grid
    #loop through the list of x- and y-values
    for index in range(len(x_list)):
        i = x_list[index];
        j = y_list[index];
        #print([i,j]);
        
        #METHOD 2: LEGIT? #A-D are the four vertices of the moving grid
        A = [i-0.5,j-0.5];
        B = [i+0.5,j-0.5];
        C = [i+0.5,j+0.5];
        D = [i-0.5,j+0.5];
        
        grid[int(A[0]),int(A[1])] -= abs(0.25*(int(A[0]+1)-A[0])*(int(A[1])+1-A[1])); #for top left
        grid[int(A[0]),int(A[1])+1] -= abs(0.25*(int(A[0]+1)-A[0])*(int(A[1])+1-A[1]));
        grid[int(A[0])+1,int(A[1])] -= abs(0.25*(int(A[0]+1)-A[0])*(int(A[1])+1-A[1]));
        grid[int(A[0])+1,int(A[1])+1] -= abs(0.25*(int(A[0]+1)-A[0])*(int(A[1])+1-A[1]));

        grid[int(B[0]),int(B[1])] -= abs(0.25*(int(B[0])-B[0])*(int(B[1])+1-B[1])); #for top right
        grid[int(B[0]),int(B[1])+1] -= abs(0.25*(int(B[0])-B[0])*(int(B[1])+1-B[1]));
        grid[int(B[0])+1,int(B[1])] -= abs(0.25*(int(B[0])-B[0])*(int(B[1])+1-B[1]));
        grid[int(B[0])+1,int(B[1])+1] -= abs(0.25*(int(B[0])-B[0])*(int(B[1])+1-B[1]));

        #for bottom right:
        grid[int(C[0]),int(C[1])] -= abs(0.25*(int(C[0])-C[0])*(int(C[1])-C[1])); 
        grid[int(C[0]),int(C[1])+1] -= abs(0.25*(int(C[0])-C[0])*(int(C[1])-C[1]));
        grid[int(C[0])+1,int(C[1])] -= abs(0.25*(int(C[0])-C[0])*(int(C[1])-C[1]));
        grid[int(C[0])+1,int(C[1])+1] -= abs(0.25*(int(C[0])-C[0])*(int(C[1])-C[1]));

        grid[int(D[0]),int(D[1])] -= abs(0.25*(int(D[0])+1-D[0])*(int(D[1])-D[1])); #for bottom left
        grid[int(D[0]),int(D[1])+1] -= abs(0.25*(int(D[0])+1-D[0])*(int(D[1])-D[1]));
        grid[int(D[0])+1,int(D[1])] -= abs(0.25*(int(D[0])+1-D[0])*(int(D[1])-D[1]));
        grid[int(D[0])+1,int(D[1])+1] -= abs(0.25*(int(D[0])+1-D[0])*(int(D[1])-D[1]));
        
        '''
        #METHOD 1: PIECE OF TRASH: 
        #assign charge to all four corners of same cell
        grid[int(i),int(j)] -= 0.25; #top left
        grid[int(i),int(j)+1] -= 0.25;
        grid[int(i)+1,int(j)] -= 0.25;
        grid[int(i)+1,int(j)+1] -= 0.25; #bottom right
        '''
    #plt.imshow(grid, cmap='viridis', interpolation='nearest');
    #plt.colorbar();
    #plt.xlabel("x");
    #plt.ylabel("y");
    #plt.show();
    return grid

'''
PART 2(b)
'''
def relaxation(): #I stole my own code for problem 1 and modified it for this one
    eps = 1; #epsilon_0. set to 1, for now
    rho = 1; #wtf is rho, here?
    M = 100; 
    h = int(100/M); #h is the separation between adjacent grids, aka dt.
    tol = 1e-10; #goddamn.
    phi_array = zeros([M+1,M+1],float); #create a grid to hold the phi values
    phi_prime_array = empty([M+1,M+1],float);
    #BCs: V = 0 on all 4 sides.
    phi_array[0,:] = 0;
    phi_array[:,0] = 0;
    phi_array[M,:] = 0;
    phi_array[:,M] = 0;
    delta = 1.0;
    rho = cmap();
    while delta > tol:
        #Get the new values of phi
        for x in range(M+1):
            for y in range(M+1):
                if x == 0 or y == 0 or x == M or y == M:
                    phi_prime_array[x,y] = phi_array[x,y];
                else:
                    phi_prime_array[x,y] = 0.25*( \
(rho[x,y]*(h**2)/eps) \
+ phi_array[x-h,y] + phi_array[x+h,y] + phi_array[x,y-h] \
+ phi_array[x,y+h]);
            #get delta so we know if we can stop
        dif = abs(phi_array - phi_prime_array);
        delta = max(dif);
        print(delta);
        if delta == 0:
            delta = to1*10;#small but not small enough
        #switch the two arrays
        phi_array, phi_prime_array = phi_prime_array, phi_array;
    #plot stuff
    plt.imshow(phi_array, cmap='viridis', interpolation='nearest')
    plt.xlabel("x");
    plt.ylabel("y");
    plt.colorbar()
    plt.show();
    return

'''
part2(c)
'''
def gseidel(omega):
    counter = 0; #this will keep count of the number of iterations it takes to converge
    eps = 8.85e-12; #permittivity in free space
    M = 100; 
    #omega = 0.96; #0 would be no overrelaxation
    h = int(100/M); 
    tol = 1e-10; 
    phi_array = zeros([M+1,M+1],float); #only one array, now.
    delta = 1.0;
    rho = cmap();
    rho*= 1.6e-19; #multiply by unit charge
    while delta > tol:
        delta = 0.0;
        counter += 1; #increment this guy
        for x in range(M):
            for y in range(M):
                phi_array_old = phi_array[x,y];
                phi_array_new = (1+omega)*(0.25)*(phi_array[x+h,y] + phi_array[x-h,y] + phi_array[x,y+h] + phi_array[x,y-h] + rho[x,y]*(h**2)/eps)\
- omega*phi_array_old; #eq 9.17 in the book.
                phi_array[x,y] = phi_array_new;
                dif = abs(phi_array_new - phi_array_old); 
                if dif > delta:
                    delta = dif;
    #plot stuff
    #plt.imshow(phi_array, cmap='viridis', interpolation='nearest')
    #plt.xlabel("x");
    #plt.ylabel("y");
    #plt.colorbar()
    #plt.show();
    return counter #return the number of iterations, so we can use it for the golden ratio search.

if "__name__" == "__main__":
    relaxation()

def grs(): 
    Z = ((5**(1/2))+1)*0.5;
    x_1 = 0.8;
    x_4 = 0.99;
    delta = 1.0;
    omega_list = np.array;
    iteration_list = list();
    omega_list = list();
    while delta > 0.001:
        x_3 = x_1 + (x_4-x_1)/Z;
        x_2 = x_4 - (x_4-x_1)/Z;
        if gseidel(x_2)<gseidel(x_3):
            #define a new x4
            x_4 = x_3;
        else:
            #new x1
            x_1 = x_2;
        omega = (x_4+x_1)/2;
        num_it = gseidel(omega)
        iteration_list.append(num_it);
        omega_list.append(omega);
        delta = x_4-x_1;
        print([omega,num_it]);
    plot(omega_list, iteration_list, 'k.')
    xlabel("$\omega$");
    ylabel("# iterations");
    show() 
    return

if __name__ == "__main__":
    grs()
