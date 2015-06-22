#Learning NETWORKX

import networkx as NX
import matplotlib.pyplot as plt

#Create graph where to put nodes
G=NX.Graph()

#Try adding nodes - UNIQUE IDs
G.add_node("John")
G.add_node("Jose")
G.add_node("Mara")

#Now add edges to connect nodes - State what is it adding between, takes two elements
G.add_edge("John", "Mara")
G.add_edge("Jose", "Mara")

#Print some network info
print G.number_of_nodes()
print G.number_of_edges()

print G.nodes()
print G.edges()

#Calculate degree of specific nodes, or list degrees of all nodes - DICTIONARY STYLE or int() if it is just one that is being asked for
print G.degree("Mara") #TO GET EDGE DENSITY 
print type(G.degree("Mara"))
print G.degree()

######CREATING A DIRECTED GRAPH
D=NX.DiGraph()

#Add edges in batches as list of tuples - ADDING NODES ON THE FLY
D.add_edges_from([("A", "B"), ("C", "A")]) #FORMAT IN TUPLE IS FROM FIRST ELEMENT TO SECOND ELEMENT (out-degree, in-degree)

print D.in_degree()
print D.out_degree()

#Get neighbors for specific nodes as List, SINCE DIRECTED GRAPH, will print out-degree neighbors
print D.neighbors("A")
print D.neighbors("B")

######TRANSFORM TO UNDIRECTED GRAPH
UD=D.to_undirected()

print UD.neighbors("B") #Now since it is an undirected graph it gets all in- and out- degree neighbots, if it has them
print UD.neighbors("A")

######CREATING NETWORKS DIRECTLY FROM TEXT FILES

#FIRST - Edge list format (simple 2 columns files, one edge per line, first in line is source node)
FIRST=NX.read_edgelist("DATABASES/test.edge")
print FIRST.edges() #NOT SURE WHAT "u" CHARACTERS ARE, BUT APPARENTLY DON'T AFFECT ITERATING OF STRUCTURE, SO FOR NOW IGNORE THEM

#SECOND - Adjecent list format (first in line is source node, rest are considered target nodes)
SECOND=NX.read_adjlist("DATABASES/test.adj1")
print SECOND.edges() 

######CREATING WEIGHTED GRAPHS
W=NX.Graph()

W.add_edge("a", "b", weight=5) #weights are addedd as a third argument to "on-the-fly" node addition
W.add_edge("a", "c", weight=10)
W.add_edge("a", "d", weight=4)
W.add_edge("b", "d", weight=1)

print W.edges(data=True) #'data=True' to include "weight" as a third component

#One pythonic way of filter edges that have weights greater than a threshold
print filter(lambda x: int(x.split()[2])>3, NX.generate_edgelist(W, data=["weight"])) #generate_edgelist, makes list out of network graph elemnts

#NX way to filter edges that pass a threshold
W_T=[(u,v) for (u,v,d) in W.edges(data=True) if d["weight"]>3]
print W_T

#get degree as usual, regardless of weight - TO GET EDGE DENSITY
print W.degree("a") #When no weight stated then all edge's weights are equal to 1

#get weighted degree which is the sum of the weights of all edges
print W.degree("a",weight="weight") #specify with second argument that you need weight

#DRAWING NETWORK
NX.draw(W, NX.shell_layout(W))
plt.show()