import numpy as np
from numpy import empty, array, arange
from pylab import plot,show,xlabel,ylabel

#write a function that takes an array of x,y values and outputs their derivatives
def evaluator(X):
    a = 1;
    b = 3;
    x = X[0];
    y = X[1];
    xprime = 1-(b+1)*x + a*x**2*y;
    yprime = b*x - a*x**2*y;
    Y = array([xprime, yprime],np.float64);
    return Y

#a bunch of global variables
N = 1000;
H = 20/N;
tpoints = arange(0,20,H);
tol = 10**-10;
xpoints, ypoints = ([] for i in range(2));

def bs():
    X = array([0,0],np.float64);#change initial conditions here
    for t in tpoints:
        print(t);
        xpoints.append(X[0]);
        ypoints.append(X[1]);
        #first and second evaluation using midpoint step
        n = 1;
        X1 = X + 0.5*H*evaluator(X);
        X2 = X + H*evaluator(X1);
        #define a bunch of arrays that contain the Richardson extrapolation
        X1_arr = empty([n,2],np.float64);
        X1_arr[0] = 0.5*(X1 + X2 + 0.5*H*evaluator(X2));
        #evaluate error
        err = 2*H*tol;
        while err>H*tol:
            n+=1;
            h = H/n;
            X1 = X + 0.5*h*evaluator(X);
            X2 = X + h*evaluator(X1);
            for i in range(n-1):
                X1 += h*evaluator(X2);
                X2 += h*evaluator(X1);
            #get extrapolation estimates
            X2_arr = X1_arr;
            X1_arr = empty([n,2],np.float64);
            X1_arr[0] = 0.5*(X1 + X2 + 0.5*h*evaluator(X2));
            for m in range(1,n):
                epsilon = (X1_arr[m-1]-X2_arr[m-1])/((n/(n-1))**(2*m)-1);
                X1_arr[m] = X1_arr[m-1] + epsilon;
            err = abs(epsilon[0]);
        X = X1_arr[n-1];
    #plot results
    plot(tpoints,xpoints,tpoints,ypoints);
    xlabel('t');
    ylabel('Concentrations');
    show();
    return
        
        

    

    
