import functions as func
from archive.network_sh import edgesData

# The list of all proteins that are linked to a particular protein
func.getList(edgesData, "Q16787")

# Create a histogram for a protein list
degrees = func.getDegreeAndHistogram(
    edgesData, ["P23229", "P00519", "P31749", "P42574"])

# Sort the list by the degree of its proteins
func.getOrderedDegree(degrees)
