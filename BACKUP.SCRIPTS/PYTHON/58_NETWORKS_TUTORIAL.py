import networkx as NX
import matplotlib.pyplot as plt

#Creating a random network
G1=NX.erdos_renyi_graph(50,0.3) # n number of nodes, probability of p of edge between nodes independent of every other edge

#Draw network
NX.draw(G1, NX.shell_layout(G1))
plt.show()

#Get node degree distribution
DD=[]
for node in G1.nodes():
    DD.append(G1.degree(node))
plt.hist(DD)
plt.show()
    
#Get average shortest path length(Characteristic path length (L))
L=NX.average_shortest_path_length(G1)
print L

#Get the average clustering coefficient of the graph (CC)
CC=NX.average_clustering(G1)
print CC


#Creating a small-world network
SW=NX.watts_strogatz_graph(50,6,0.3) #model for small world network with parameters (nodes, number of nearest neighbors each node is connected
                                    #to, probability of rewiring each edge)

#Draw network
NX.draw(SW, NX.shell_layout(SW))
plt.show()

#Get L
L2=NX.average_shortest_path_length(SW)
print L2

#Get CC
CC2=NX.average_clustering(SW)
print CC2

#Get node degree distribution
DD=[]
for node in SW.nodes():
    DD.append(SW.degree(node))
plt.hist(DD)
plt.show()


#Creating a scale-free graph
SF=NX.scale_free_graph(50) #number of nodes

#Draw network
NX.draw(SF, NX.shell_layout(SF))
plt.show()

#Get L
L2=NX.average_shortest_path_length(SF)
print L2    

#Get node degree distribution
DD=[]
for node in SF.nodes():
    DD.append(SF.degree(node))
plt.hist(DD)
plt.show()