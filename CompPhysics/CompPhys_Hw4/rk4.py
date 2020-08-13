from math import pi
import numpy as np
from numpy import array, empty, arange
from matplotlib.pyplot import plot, show

def f(r): #r is an array (x,x',y,y')
    A = 1;
    B = 1;
    x = r[0];
    v_x = r[1];
    y = r[2];
    v_y = r[3];
    R = array([x,y],np.float64);
    V = array([v_x,v_y],np.float64);
    #want to return the velocities and accelerations (x', x'', y', y'')
    '''without dynamical friction
    a = (-0.25*R)/((x**2 + y**2)**(3/2)); #a=r''=-1/4*r/(|r|^3)'''
    #with dynamical friction:
    a = (-0.25*R)/((x**2 + y**2)**(3/2))-((A/(B+(v_x**2 + v_y**2)**(3/2)))*V);
    return array([v_x, a[0], v_y, a[1]],float)

def get_rho(r_big,r_small,h,delta): #delta is target accuracy. 10^-4?
    err = ((r_big[0]-r_small[0])**2+(r_big[2]-r_small[2])**2)**(1/2);
    rho = 30*h*delta/err;
    return rho

def rk4(): 
    r = array([1.0,0.0,0.0,0.4],np.float64); #initial conditions
    t0 = 0.0;
    t = t0;
    tf = 5.0; #80 gives approximately 10 orbits
    N =1000;
    h = (tf-t0)/N;
    tpoints = [t0];
    rpoints = [1.0];
    xpoints = [1.0];
    ypoints = [0.0];

    r_val = 1; #|r|
    r_min = 1;
    while r_val>10**-7: #condition with dynamical friction: r_val>10**-7:
        #do adaptive step size
        #One big step of size 2h (to check for adaptive):
        #print("{%.10f}" %t);
        
        #output some stuff so we know it's running

        if r_val<r_min:
            print(np.log10(r_min)); #need this to get to 10^-7.
            r_min = r_val;

        H = 2*h;
        k1 = H*f(r);
        k2 = H*f(r+0.5*k1);
        k3 = H*f(r+0.5*k2);
        k4 = H*f(r+k3);
        r_big = r+1/6*(k1 + 2*k2 + 2*k3 + k4);
        
        #Two small steps of sizes h (to check for adaptive):
        #first small step
        k1 = h*f(r);
        k2 = h*f(r+0.5*k1);
        k3 = h*f(r+0.5*k2);
        k4 = h*f(r+k3);
        r_small = r+1/6*(k1 + 2*k2 + 2*k3 + k4);
        #second small step
        k1 = h*f(r_small);
        k2 = h*f(r_small+0.5*k1);
        k3 = h*f(r_small+0.5*k2);
        k4 = h*f(r_small+k3);
        r_small += 1/6*(k1 + 2*k2 + 2*k3 + k4);

        #now, do adaptive stuff
        rho = get_rho(r_big,r_small,h,10**-6); #check for the quadratures
        if rho>1: #error is small
            if rho**(1/4)>2:
                rho = 16;
            r = r_small;
            #get magnitude |r|:
            r_val = (r[0]**2+r[2]**2)**(1/2);
            rpoints.append(np.log10(r_val));
            #update time.
            t += 2*h;
            tpoints.append(t);
            xpoints.append(r[0]);
            ypoints.append(r[2]);
            h*=rho**(1/4);
        else:
            h*=rho**(1/4);
        #print(rho);

        
    #plot(xpoints, ypoints); 
    plot(tpoints, rpoints)
    show();
    return

if __name__ == "__main__": 
    rk4()
