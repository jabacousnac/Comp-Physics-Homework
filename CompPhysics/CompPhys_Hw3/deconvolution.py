import numpy as np
import matplotlib.pyplot as plt
from scipy import loadtxt, exp, empty
from numpy.fft import rfft2, irfft2

sigma = 25;
def showblurred(): #this is my function to display the blurred image
    arr = np.loadtxt("blur.txt",delimiter=",");
    #plt.imshow(arr, cmap='gray', interpolation='nearest')
    #plt.show();
    return arr

def gaussian(y,x): #calculate gaussian
    sigma = 25;
    return exp(-(x**2+y**2)/(2*sigma**2))
    
def psf(): #this is the function to generate the point spread function 
    arr = np.loadtxt("blur.txt",delimiter=",");
    dim = arr.shape;
    gauss_array = empty([dim[1], dim[0]], float)  
    for y in range(0,dim[0]):
        for x in range(0,dim[1]):
            gauss_array[y,x] = gaussian((y+dim[0]/2)%dim[0]-dim[0]/2,(x+dim[1]/2)%dim[1]-dim[1]/2);
    plt.imshow(gauss_array, cmap='gray',interpolation='nearest');
    plt.show();
    return gauss_array

def decon():
    #first, perform FTs on the two arrays
    conv = rfft2(showblurred());
    gauss = rfft2(psf());
    arr = np.loadtxt("blur.txt",delimiter=",");
    dim = arr.shape;
    counter =0; #this is to see how many zeros we encounter in psf
    #now, we can perform deconvolution
    sharp_array = empty([dim[0], dim[1]//2 + 1], complex);
    epsilon = 10 ** -3;
    for x in range(dim[1]//2 + 1): #nearest integer rounded down
        for y in range(dim[0]):
            #take care of case when gaussian is zero
            if abs(gauss[y,x]) < epsilon: 
                gauss[y,x]=1
                counter+=1;
            #deconvolve:
            sharp_array[y,x] = conv[y,x]/gauss[y,x];
    final_arr = irfft2(sharp_array); #IFT
    plt.imshow(final_arr, cmap = 'gray', interpolation = 'nearest');
    plt.show();
    print(counter);
    return

