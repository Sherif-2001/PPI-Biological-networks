import numpy as np
import pandas as pd
import functions as func
import networkx as nx

# the new function

newedges = pd.read_csv("finished.csv",header=0,usecols=["head","tail","edge_weight"])

head = newedges["head"].to_list()
tail = newedges["tail"].to_list()
weight = newedges["edge_weight"].to_list()

tempData = []
for i in np.arange(100):
    tempData.append((head[i],tail[i],weight[i]))

# Initialize a graph with edges, name, and graph attributes
network = nx.DiGraph(name="PPI Network")

# Add given edges data to the network
network.add_weighted_edges_from(tempData)

# Draw the given network 
func.drawNetworkGraph(network)

# Get the shortest path(s) between the two given proteins
func.getShortestPaths(network,"Q16787","P24043")