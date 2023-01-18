import functions as func
import networkx as nx

edgesData = func.getNetworkData("Alpha6Beta4Integrin-edges.txt")

network = nx.Graph(name="PPI Network")
network.add_edges_from(edgesData,weight=1)

func.drawNetworkGraph(network)
func.getShortestPaths(network,"Q16787","P24043")

















# DON'T DELETE

# networkData = open("PathLinker_2018_human-ppi-weighted-cap0_75.txt", "r")
# newNetwork = nx.Graph(name="NEw")
# for data in networkData:
#     if data.split("\t")[0] == "#tail":
#         continue
#     weight = float(data.split("\t")[2])
#     newNetwork.add_edge(data.split("\t")[0],data.split("\t")[1],weight= weight)

# nx.draw(newNetwork, with_labels=True)
# plt.figure(figsize=(3.841, 7.195), dpi=100)
# plt.savefig('NEW Network.png')
