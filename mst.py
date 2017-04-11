import fiona
from shapely.geometry import shape
from sys import argv
from collections import deque

def extract_shapes(edgeName, nodeName):
    edges = {line['id']: shape(line['geometry']) for line in list(fiona.open(edgeName) ) }
    nodes = [shape(point['geometry']) for point in list(fiona.open(nodeName) ) ]
    return edges,nodes


def make_network(edges):
    """returns a datastructure that contains edge ids along with information about connections"""
    edgeIds = [key for key in edges.keys()]
    network = {edge: [] for edge in edgeIds}
    
    for e1 in range(len(edgeIds)):
        e1Id = edgeIds[e1]
        for e2 in range(e1+1, len(edgeIds)):
            e2Id = edgeIds[e2]
            if edges[e1Id].touches(edges[e2Id]):
                network[e1Id].append(e2Id)
                network[e2Id].append(e1Id)

    return network
           

def get_important_edges(edges, nodes):
    """returns the edges that are closest to the points that need to be hit in the MST"""
    importantEdges = set()
    
    for node in nodes:
        closest = None
        minDist = None
        for edge in edges:
            distance = node.distance(edges[edge])
            if (minDist == None or distance < minDist):
                minDist = distance
                closest = edge
                #TODO: if there are multiple segments in the edge find which segment we want / maybe not depends on what fiona considers a segment
        importantEdges.add(closest)
    return importantEdges


def find_paths(network, importantEdges, edges):
    """finds the path and cost of that path between all of the edges we need in our MST"""
    paths = {}#(e1, e2): {"cost": cost, "path": [edges]}
    for e in importantEdges:
        seenEdges = {e2: 0 for e2 in network[e]}
        seenEdges[e] = -1
        edgeQueue = deque([(e2,0,[]) for e2  in network[e]])#(the node to visit, what it cost to get there, what we've visited so far)
        while(len(edgeQueue)):
            current = edgeQueue.popleft()
            if (current[0] in importantEdges):
                if tuple(sorted([e, current[0]])) not in paths:
                    paths[tuple(sorted([e, current[0]]))] = {'cost': current[1], 'path': current[2]}
            else:
                accessable = network[current[0]]
                for nextEdge in accessable:
                    costToNext = current[1] + edges[current[0]].length
                    if nextEdge not in seenEdges or costToNext < seenEdges[nextEdge]:
                        seenEdges[nextEdge] = costToNext
                        edgeQueue.append((nextEdge, costToNext, current[2]+[current[0]]))
    return paths


def get_edges_in_route(network, edges, start, mst):
    """used to connect a 'node' to a pre-existing mst, returns the shortest path to get to any point that is already in the mst"""
    seenEdges = {e2: 0 for e2 in network[start]}
    seenEdges[start] = -1
    edgeQueue = deque([(e2,0,[]) for e2  in network[start]])#(the node to visit, what it cost to get there, what we've visited so far)
    while(len(edgeQueue)):
        current = edgeQueue.popleft()
        if (current[0] in mst):
            return current[2]
        else:
            accessable = network[current[0]]
            for nextEdge in accessable:
                costToNext = current[1] + edges[current[0]].length
                if nextEdge not in seenEdges or costToNext < seenEdges[nextEdge]:
                    seenEdges[nextEdge] = costToNext
                    edgeQueue.append((nextEdge, costToNext, current[2]+[current[0]]))
    return []

def find_mst(paths, importantEdges, network, edges):#Kruskal's algorithm
    """returns the edges we need to have for a minimum spanning tree"""
    sortedEdges = deque(sorted([{"connects":e,"cost":paths[e]["cost"],"path":paths[e]["path"]} for e in paths], key= lambda x: x["cost"]))
    forest = []
    for e in importantEdges:
        mst = set(e)
        forest.append({"nodes":[e],"mst":mst})
            
    while len(sortedEdges):
        potentialEdge = sortedEdges.popleft()
        n1 = potentialEdge["connects"][0]
        n2 = potentialEdge["connects"][1]

        #forest indices
        f1 = None
        f2 = None

        other = None
        for i in range(len(forest)):
            if n1 in forest[i]:
                other = n2
                f1 = i
            if n2 in forest[i]:
                other = n1
                f1 = i
            if other is not None:
                for j in range(i, len(forest)):
                    if other in forest[j]:
                        f2 = j
                        break
                break
        
        if f1 is None and f2 is not None:
            forest[f2]["nodes"].append(other)
            forest[f2]["mst"].update(get_edges_in_route(network, edges, other, forest[f2]["mst"]))
        if f2 is None and f1 is not None:
            forest[f1]["nodes"].append(other)
            forest[f1]["mst"].update(get_edges_in_route(network, edges, other, forest[f1]["mst"]))
        if f1 is not None and f2 is not None and f1 != f2:
            forest[f1]["nodes"].extend(forest[f2]["nodes"])
            forest[f1]["mst"].update(forest[f2]["mst"])
            del forest[f2]
            forest[f1]["mst"].update(potentialEdge["path"])

        if len(forest[0]["nodes"]) == len(importantEdges):
            break

    return forest[0]["mst"]


def write_file(outFileName, edgeIds, edges, copyParamsFromFile):
    with fiona.open(copyParamsFromFile) as source:
        sourceDriver = source.driver
        sourceCrs    = source.crs
        sourceSchema = source.schema

    with fiona.open(outFileName, 'w',
                    driver=sourceDriver,
                    crs=sourceCrs,
                    schema=sourceSchema) as outFile:
        for edge in edgeIds:
            outFile.write(edges[edge])

if (len(argv) != 4):
    print("usage: python mst.py edges nodes outfile")
    exit(2)
edgesFileName = argv[1]
nodesFileName = argv[2]
outFileName   = argv[3]

edgeLines, nodes = extract_shapes(edgesFileName, nodesFileName)
edges = {line['id']: line for line in list(fiona.open(edgesFileName) ) }

print("Calculating the network, please be patient")
network = make_network(edgeLines)#this is the slow part

print("Determining which edges need to be in the network")
importantEdges = get_important_edges(edgeLines, nodes)

print("Finding the paths that connect the edges")#essentially turning the network into a graph
paths = find_paths(network, importantEdges, edgeLines)

print("Finding the minimum spanning tree")
mst = find_mst(paths, importantEdges, network, edgeLines)

print("writing to "+outFileName)
write_file(outFileName, mst, edges, edgesFileName)
