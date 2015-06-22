#Learning Hierarchical Clustering and Heatmap
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
from scipy.cluster import hierarchy

os.chdir("NETWORK/DIRECTED_NETWORK")

#Open sample file
FILE_IN=open("hier.sample").read().splitlines()

#Save the headers into array
colheaders=FILE_IN[0].strip().split()[1:]
print colheaders #Conditions
rowheaders=[]

datamatrix=[]

for line in FILE_IN[1:]:
    data=line.strip().split()
    rowheaders.append(data[0])
    datamatrix.append([float(x) for x in data[1:]])

print datamatrix
print rowheaders #genes

#Convert native data array (python) into numpy array
datamatrix=np.array(datamatrix)
print datamatrix

#Now that we have the numpy array we can get the pairwise distances of the data - THIS IS THE CONDENSED DISTANCE MATRIX
distancematrix=distance.pdist(datamatrix, "euclidean") #IMPORTANT - The pairwise distances of two vectors are the 
                                                        #the distances based on the attributes (i.e. coordinates)  
                                                        #contained in the vector. This can be calculated using euclidean
                                                        #distances of vectors or other formulas (look up pdist for more
                                                        #info)
print "distance matrix:",distancematrix #Since we had vectors in this example (5 genes), there are exactly 10 possible
                                        #pairwise comparissons (N*(N-1)/2). THIS IS THE "CONDENSED DISTANCE MATRIX"
                                        
#IMPORTANT - Apparently it is wrong to convert this to a redundant squareform matrix using distance.squareform(CDM). 
#The subsequent function's documentation hierarchy.linkage() until recent has wrongly stated that it can take a redundant
#squareform matrix (RSM), apparently converting the CDM to RSM will make the subsequent linkage function think that the 
#CDM is the original set of observations, which is wrong. We will directly input the CDM into linkage function.
#http://blog.nextgenetics.net/?e=44

#Now perform hierarchical/agglomerative clustering in the CDM
linkagematrix=hierarchy.linkage(distancematrix, "single") #the method chosen varies, default=single, look in documentation
print "linkage matrix:\n", linkagematrix
#Output looks like this (explanation of index next to "->" not in output):    
#[[ 0.          1.          1.94190068  2.        ] ->Index 0 + 1 make cluster Index 5
#[ 3.          4.          2.0903242   2.        ]  ->Index 3 + 4 make cluster Index 6
#[ 2.          5.          3.30441329  3.        ]  ->Index 2 + 5 make cluster Index 7
#[ 6.          7.          3.81052537  5.        ]] ->Index 6 + 7 makes ROOT (which contains all clusters, all vectors)
#First two columns specify the indexes, the third column the linkage distance and the last one the number of vectors(genes)
#in the original dataset in that leaf(cluster) of the tree. The indexes are labeled based on their position of the original 
#dataset, we had five vectors so they are numbered [0,1,2,3,4].The third as a linkage distance between vector 2 and cluster
#5, which contains 2 vectors, making the total number of vectors in that cluster equal to 3. 

#PLOT - Now we can plot directly the linkage matrix, this is what we plot
hierarchy.dendrogram(linkagematrix, labels=rowheaders) #rowheaders match to it because it is done by matching the indeces
                                                        #of the data presented in the linkage matrix
plt.show()

#THE FOLLOWING SECTION IS TO ORDER THE DATA TO MATCH IT TO HYPOTHETICAL HEATMAP PLOTTED NEXT TO IT
#IF, we are required to generate a heatmap then we would need to reorder the rows of the original data according to the 
#clustering(leaves in dendogram, look at plot):
#gene    D   E   C   A   B
#index   3   4   2   0   1

heatmaporder=hierarchy.leaves_list(linkagematrix) #Takes linkage matrix (what we use for the dendogram) and produces an array
                                                    #of the order of leaves in it (leaves_list)
print heatmaporder

#Now order originally data according to the leaves list. This works because of numpy slicing annotation 
ordereddatamatrix=datamatrix[heatmaporder,:] #Rows ordered as vector, all columns
print ordereddatamatrix

#Do the same for row headers
rowheaders=np.array(rowheaders)
orderedrowheaders=rowheaders[heatmaporder,:]
print orderedrowheaders #genes

#CREATE HEATMAP
#Start arguments
fig,ax=plt.subplots() #subplots() returns 2 tuples (that's why we need 2), "fig" is the figure object and "ax" would be
                        #a single axis object or an array of axis objects, NOW WE CAN CHANGE THEM
                        #For documentation and usage look up matplotlib.figure and matplotlib.axes

#Set labels for each tick in box (column and row in heatmaps) #Look at your matrix to decide where things go
ax.set_xticklabels(colheaders) #Condition label
ax.set_yticklabels(orderedrowheaders) #Gene labels

#Set ticks at the middle of the box (so that the label is in the center)
ax.xaxis.set_ticks_position("bottom") #with respect to frame
ax.yaxis.set_ticks_position("right") #with respect to frame
ax.set_xticks(np.arange(len(colheaders))+0.5) # np.arange creates an array of values from zero till before the length 
                                                #specified, each block is treated with coordinates as part of the axis,
                                                #having a unit of 1, so the function is saying:
                                                #"place the ticks at 0(0+0.5), 1.5,2.5 and 3.5, the middle of each block 
ax.set_yticks(np.arange(len(orderedrowheaders))+0.5)
ax.xaxis.set_ticks_position("none") #To delete the little ticks
ax.yaxis.set_ticks_position("none")

#Set titles
ax.set_title("heatmap of gene expression")
ax.set_xlabel("conditions")
ax.set_ylabel("genes")

heatmap=plt.pcolormesh(ordereddatamatrix,cmap=plt.cm.Blues, norm=None)
plt.show()

