#LEARNING TO MAKE HEAT MAPS

#To makea  heatmap you need to have a matrix of values, such as gene expression levels against different types of cancers

import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

#These are the labels, would need to know how to use them later
column_labels=list("ABCD")
row_lables=list("XWYZ")

#Create random array for example
data=np.random.rand(10,10) #random.rand from numpy creates an array of an unifrom distribution over [0,1] in the given matrix shape

print data
print data[2][3] #[row][column] format

#Make heatmap of array
heatmap=plt.pcolor(np.flipud(data), cmap=cm.Blues) #first arguement is your array and second argument is the color map argument
                                                    #IMPORTANT - Use np.flipud() for visualization purposes of ordered rows. The heatmap is
                                                    #normally plotted from the [0] row to the [n] row, which would look upside down when 
                                                    #comparing the array row order to the graph. np.flipud() flips the row order in the up/down
                                                    #direction for easier visualization
plt.show()