#Learning Numpy

import numpy as np
from scipy.spatial import distance
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

A=np.array([[0, 1, 2, 3, 4],[10,20,30,40,50]]) #Make array with list of elements
print A
print A[1][3] #Calls row-1,column-3
print A.shape #Gives you dimensions of array

B=np.arange(2,20,3) #Creates evenly spaced array (start,end,step size), interval does not include last value
print B

C=np.zeros((4,4)) #Emtpy 4x4 matrix of zeros - (4,4) is a tuple
print C

D=np.eye(3) #Diagonal array of 1s
print D

E=np.diag(np.array([1,2,3,4])) #Calls for diagonal array of elements in np.array
print E

F=np.linspace(2,4,8)#Simlar to np.arange but (start, end, number of elements in array
print F

G=np.random.rand(3,5) #To obtain number of random values in the range [0,1]
print G

H=np.arange(0,16,1).reshape(4,4) #First creates one-dimensional array, then makes it into 2-dimensional array with reshape
print H

np.fill_diagonal(H,0) #Function updates array with values in diagonal (array, value to be added to diagonal)
print H

I=np.array([[0,1,2,3],  #To use distance.squareform we need to have a symmetric matrix
           [1,0,6,7],
           [2,6,0,8],
           [3,7,8,0]])
print I
J=distance.squareform(I) #Converts the redundant squareform matrix (RSM) to the condensed distance matrix (CDM)
print J

hier=hierarchy.average(J)
hierarchy.dendrogram(hier)
plt.show()