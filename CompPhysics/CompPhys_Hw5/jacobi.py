import numpy as np
from numpy import zeros, empty, array, arange
from pylab import plot, imshow, show, xlabel, ylabel

def jacobi():
    g = 9.81; #gravity on earth
    M = 100; 
    h = 10.0/M; #h is the separation between adjacent grids, aka dt.
    tol = 10**-6;
    z_array = zeros([M+1], float); #create a grid to hold the z values
    t_array = zeros([M+1], float);
    z_prime_array = empty([M+1], float);
    #BCs:
    z_array[0] = 0;
    z_array[100] = 0;
    delta = 1.0;
    while delta > tol:
        #Get the new values of z
        for t in range(M+1):
            t_array[t] = t/10;
            if t==0 or t==M:
                z_prime_array[t] = z_array[t];
            else:
                z_prime_array[t] = 0.5*((g*(h)**2)+z_array[t-1]+z_array[t+1]);
        #get delta so we know if we can stop
        delta=max(abs(z_array-z_prime_array));
        print(delta);
        #switch the two arrays
        z_array, z_prime_array = z_prime_array, z_array;
    #plot stuff
    plot(t_array,z_array);
    xlabel("t (s)");
    ylabel("z (m)");
    show();
    return

if __name__ == "__main__":
    jacobi()


    

    
