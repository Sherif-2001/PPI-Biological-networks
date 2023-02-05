from matplotlib import pyplot as plt
import networkx as nx
import numpy as np

from path_class import PathData

def getNetworkData(filePath):
    """
    Get the network data edges from given [filePath]

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
        tempData.append((row.split("\t")[0],row.split("\t")[1],int(row.split("\t")[2])))
    return tempData

def drawNetworkGraph(network):
    """
    Draw the given [network] and save it in image file

    Parameters
    ----------
    network : NetworkX graph or list of nodes
        the network of proteins to be drawn
    """
    pos = nx.spring_layout(network,k=1)
    nx.draw(network, with_labels=True,pos=pos)
    plt.show()

def getShortestPaths(network,head,tail):
    """
    Get the shortest path(s) between [head] and [tail] in the given [network]

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
    tempPathsData = []
    tempTotalWeight = 100
    for path in shortestPaths:
        if nx.path_weight(network, path, 'weight') > tempTotalWeight:
            break
        tempEdges = []
        for i in range(len(path) - 1):
            edge = path[i:i+2]
            tempEdges.append(tuple(edge + [nx.path_weight(network, edge, 'weight')]))
        tempTotalWeight = nx.path_weight(network, path, 'weight')
        tempPathsData.append(PathData(path, tempTotalWeight, tempEdges))
    makeSubNetwork(tempPathsData)
    writePathsToFile(tempPathsData)

def makeSubNetwork(pathsData):
    """
    Make a sub-network from the given [pathsData]

    Parameters
    ----------
    paths : list
        list of protein paths
    """
    subNetwork = nx.DiGraph(name="SubNetwork")
    tempEdges = []
    for data in pathsData:
        tempEdges += data.edges
    subNetwork.add_weighted_edges_from(tempEdges)
    pos = nx.spring_layout(subNetwork,k=1)
    nx.draw(subNetwork, with_labels=True,pos=pos)
    plt.show()

def writePathsToFile(pathsData):
    """
    Write the shortest paths, their total weights & each edge weight from [pathsData] to a text file

    Parameters
    ----------
    data : list of nodes
        list of shortest paths data

    weight : double
        the total weight of the paths
    """
    f = open("ShortestPaths.txt", "w")
    f.write("#Path\t\t\tTotal_Weight\tEdges_Weights\n")
    for data in pathsData:
        f.write(",".join(data.path)+"\t")
        f.write(str(data.total_weight)+"\t\t\t")
        edgesWeights = []
        for edge in data.edges:
            edgesWeights.append(str(edge[-1]))
        f.write(",".join(edgesWeights) + "\n")
    f.close()

def getList(data,protein):
    """
    Obtain a list of the protein's connected counterparts

    Parameters
    ----------
    data : list
        List of protein edges
    protein : str
        The protein's name
    """
    i = 0
    newFile = open("ListOfProteins.txt","w")
    newFile.write("Tail\t\tHead\tEdge_Weight\n")
    for m in np.arange(len(data)):
        if protein == data[m][0] == data[m][1]:
            continue
        if (protein == data[m][0]) or (protein == data[m][1]):
            newFile.write(data[m][0]+"\t\t"+data[m][1]+"\t"+str(data[m][2])+"\n")
            i+=1
    newFile.write("\nDegree = "+str(i))
    newFile.close

def getDegreeAndHistogram(data,list):
    """
    Draw a list's histogram after finding its proteins' degree

    Parameters
    ----------
    data : list
        List of protein edges
    list : list
        Collection of proteins

    Return
    ------
    degree:
        Degree of the proteins in this list
    """
    degree = []
    array = []
    for n in np.arange(len(list)):
        i=0
        for m in np.arange(len(data)):
            if (list[n] == data[m][0]) or (list[n] == data[m][1]):
                i+=1
                array.append(list[n])
        degree.append([list[n],i])
    plt.hist(array)
    return degree

def getOrderedDegree(array):
    """
    Obtain a text file with proteins that are ranked by degree

    Parameters
    ----------
    array : list
        Proteins listed along with their degrees
    """
    array.sort(key=lambda x:x[1],reverse=True)
    newFile = open("OrderedDegreeList.txt","w")
    newFile.write("Protein\tDegree\n")
    for n in np.arange(len(array)):
        newFile.write(array[n][0]+"\t\t"+str(array[n][1])+"\n")
    newFile.close
