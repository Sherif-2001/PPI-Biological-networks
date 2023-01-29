from matplotlib import pyplot as plt
import networkx as nx
import numpy as np

def getNetworkData(filePath):
    """
    Get the network data edges from given file

    Parameters
    ----------
    filePath : str
        the path of the data file

    Returns
    -------
    list
        list of edges of proteins from the text file
    """
    data = np.loadtxt(filePath, delimiter=',', skiprows=1, dtype=str)
    tempData = []
    for row in data:
        tempData.append((row.split("\t")[0],row.split("\t")[1]))
    return tempData

def drawNetworkGraph(network):
    """
    Draw the given network and save it in image file

    Parameters
    ----------
    network : NetworkX graph or list of nodes
        the network of proteins to be drawn
    """
    pos = nx.spring_layout(network,k=1)
    nx.draw(network, with_labels=True,pos=pos)
    plt.show()

def makeSubNetwork(paths):
    """
    Make a subnetwork from the shortest paths of two proteins

    Parameters
    ----------
    paths : list
        list of protein paths
    """
    subNetwork = nx.DiGraph(name="SubNetwork")
    tempEdges = []
    for path in paths:
        for i in range(len(path)-1):
            tempEdges.append(tuple(path[i:i+2]))
    subNetwork.add_edges_from(tempEdges)
    pos = nx.spring_layout(subNetwork,k=1)
    nx.draw(subNetwork, with_labels=True,pos=pos)
    plt.show()

def writePathsToFile(data, weight):
    """
    Write the shortest paths, their total weights & each edge weight to a text file

    Parameters
    ----------
    data : list of nodes
        list of shortest paths data

    weight : double
        the total weight of the paths
    """
    f = open("ShortestPaths.txt", "w")
    f.write("#Path\t\t\tTotal_Weight\tEdges_Weights\n")
    for path in data:
        pathJoined = ",".join(path)
        edgesWeights = ",".join(["1" for i in range(len(path)-1)])
        f.write(pathJoined+"\t")
        f.write(str(weight)+"\t\t\t")
        f.write(edgesWeights + "\n")
    f.close()

def getShortestPaths(network,head,tail):
    """
    Get the shortest path(s) between the two proteins

    Parameters
    ----------
    network : NetworkX graph or list of nodes
        the network of proteins to be drawn

    head : str
        the first node of the path
    
    tail : str
        the last node of the path
    """
    try:
        shortestPaths = nx.shortest_simple_paths(network,head,tail)   
    except:
        print(f"No path between {head} and {tail}")
        return 
    if head == tail:
        print("Please use two different proteins")
        return
    tempPaths = []
    tempWeight = 100
    for path in shortestPaths:
        if nx.path_weight(network, path, 'weight') <= tempWeight:
            tempWeight = nx.path_weight(network, path, 'weight')
            tempPaths.append(path)
        elif nx.path_weight(network, path, 'weight') > tempWeight:
            break
    writePathsToFile(tempPaths, tempWeight)
    makeSubNetwork(tempPaths)
