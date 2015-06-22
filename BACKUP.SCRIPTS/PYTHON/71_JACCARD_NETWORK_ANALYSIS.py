#CONTINUATION OF NETWORK ANALYSIS - JACCARD NETWORK

#Create Network based on Jaccard
import networkx as NX
import matplotlib.pyplot as plt
import community

#Get files and process them
NODE_EDGE_WEIGHTS=open("NETWORK/DIRECTED_NETWORK/NODE_EDGE_WEIGHTS").read().splitlines()
ZERO_DEGREE_NODES=open("NETWORK/DIRECTED_NETWORK/ZERO_DEGREE_NODES").read().splitlines()

EDGES=[]
for line in NODE_EDGE_WEIGHTS:
    EDGES.append([line.split()[0],line.split()[1],float(line.split()[2])])

#Make network
METABO_JC=NX.Graph()
#    METABO_JC.add_nodes_from(ZERO_DEGREE_NODES) #UNLOCK FOR ZERO DEGREE NODES
METABO_JC.add_weighted_edges_from(EDGES)
print METABO_JC.degree("Adenosine_monophosphate", weight="weight") #Make sure to add weight="weight" to get weighted degree
print len(METABO_JC.nodes())
"""
#Get Degree distribution weighted and unweighted
FILE_OUT1=open("NETWORK/DIRECTED_NETWORK/DEGREE_DISTRIBUTION_WEIGHTED","w")
FILE_OUT2=open("NETWORK/DIRECTED_NETWORK/DEGREE_DISTRIBUTION_UNWEIGHTED","w")
for node in METABO_JC.nodes():
    FILE_OUT1.write(str(METABO_JC.degree(node, weight="weight"))+"\n")
    FILE_OUT2.write(str(METABO_JC.degree(node))+"\n")
FILE_OUT1.close()
FILE_OUT2.close()    

#Look for connected components
print NX.is_connected(METABO_JC)
print NX.number_connected_components(METABO_JC)
Components=NX.connected_component_subgraphs(METABO_JC) #Breaking network into subgraph of found independent connect components, each
                                                        #component is a subgraph
print len(Components[0]), len(Components[1])
print Components[1].nodes()

#We could look at the Network clustering coefficients - Tells us how concentrated the neighboorhood of a node is
#First we could look at the global clustering coefficient
FILE_OUT3=open("NETWORK/DIRECTED_NETWORK/CC_DISTRIBUTION_WEIGHTED", "w")
CCL=NX.clustering(METABO_JC, weight="weight").values() #Clustering coefficient list of all nodes, function returns dictionary of each metabolite as 
                                        #key and its CC as value, so we only care about values
                                        #We can use unweighted and weighted as weight="weight"
for value in CCL: FILE_OUT3.write(str(value)+"\n")
FILE_OUT3.close()
print sum(CCL)/len(CCL)

#Get Network diameter
DIAMETER_LC=NX.diameter(Components[0])
DIAMETER_SC=NX.diameter(Components[1])
print "Diameters:", DIAMETER_LC, DIAMETER_SC

#Get Betweness Centrality and Eigenvector Centrality
FILE_OUT4=open("NETWORK/DIRECTED_NETWORK/BETWEEN_CENTRALITY_WEIGHTED","w")
BC=NX.betweenness_centrality(METABO_JC, weight="weight", normalized=True) #Returns dictionary of nodes as keys and centrality
                                                                            # as values, normalized is optional, weighted
for met,value in BC.iteritems():
    FILE_OUT4.write(met+" "+str(value)+"\n")
FILE_OUT4.close()

FILE_OUT5=open("NETWORK/DIRECTED_NETWORK/EIGEN_CENTRALITY_WEIGHTED","w")
EC=NX.eigenvector_centrality(METABO_JC) #Returns dictionary of nodes as keys and centrality as values
for met,value in EC.iteritems():
    FILE_OUT5.write(met+" "+str(value)+"\n")
FILE_OUT5.close()

#Get Maximal cliques in Network
MC=NX.find_cliques(METABO_JC) #List of all cliques as lists of all nodes in a click

print NX.graph_number_of_cliques(METABO_JC) #Number of maximal cliques in graph
print NX.graph_clique_number(METABO_JC) #Size of largest clique

MNC=NX.node_clique_number(METABO_JC) #Returns the size of the largest clique for each given node, DICTIONARY
FILE_OUT6=open("NETWORK/DIRECTED_NETWORK/MAXIMUM_CLIQUE_PER_NODE", "w") #To get the distribution of the sizes of the 
                                                                        #largest maximal clique per node
for node, clique_value in MNC.iteritems():
    FILE_OUT6.write(node+" "+str(clique_value)+ "\n")
FILE_OUT6.close()

TNC=NX.number_of_cliques(METABO_JC) #Returns the total number of maximal cliques per nodes, DICTIONARY
FILE_OUT7=open("NETWORK/DIRECTED_NETWORK/NUMBER_MAXIMAL_CLIQUE_PER_NODE", "w") #To get distribution of the number of 
                                                                                #maximal cliques per node 
for node, number_of_mc in TNC.iteritems():
    FILE_OUT7.write(node+" "+str(number_of_mc)+"\n")
FILE_OUT7.close()
"""
#Do Hierarchical Clustering #NO NEED TO AT THE MOMENT

#Do Modularity Implementation for Community Detection (Partitioning the graph)
PARTITION=community.best_partition(METABO_JC) #Returns dictionary
print PARTITION #Dictionary of metabolite as key and community number they belong to as value

for i in set(PARTITION.values()): #Iterate over values (communities
    print "Community", i
    members = [nodes for nodes in PARTITION.keys() if PARTITION[nodes] == i] #Get the members of each community into list
    print len(members), members

#DRAW NETWORK

plt.figure(figsize=(25,12))
pos=NX.spring_layout(METABO_JC)

NX.draw_networkx_nodes(METABO_JC, pos, node_color="r", node_size=100)
NX.draw_networkx_edges(METABO_JC, pos, width=0.2)
plt.show()


