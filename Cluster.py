#!/usr/bin/env python
import math
import sys

# Euclidean distance between 2 points in n-dimension space
def distance(point1, point2) :
	dim = len(point1)
	dist = 0
	for i in range(dim) :
		dist += (point1[i]-point2[i])**2
	return math.sqrt(dist)

# Minimum among distances of a point to different centers 
def minDistance(point, Centers) :
	distances = [distance(point, c) for c in Centers]
	return min(distances)

# An algorithm to find N centers to partition given points into N clusters.
# Pick the first center arbitrarily, then keep adding one center at a time until we get N centers.
# The point farthest to the current set of centers is selected as a center. 
def FarthestFirstTraversal(Points, numClusters) :
	Centers = [Points[0]]
	Points = Points[1:]

	for i in range(numClusters-1) :
		farDist = 0.0
		for point in Points :
			minDist = minDistance(point, Centers)
			if minDist > farDist :
				farDist = minDist
				newCenter = point
		Centers.append(newCenter)
		Points.remove(newCenter)

	return Centers

# Distortion is a measure of quality of centers chosen
def Distortion(Points, Centers): 
	distortion = 0
	for point in Points :
		distortion += (minDistance(point, Centers))**2
	return distortion/len(Points)

# This clustering assigns every point to a center closest to it.
def Cluster(Points, Centers) :
	numClusters = len(Centers)

	Clusters = [ [] for i in range(numClusters) ]
	for point in Points :
		minDist = sys.float_info.max
		for i in range(numClusters) :
			dist = distance(point, Centers[i])
			if dist < minDist :
				minDist = dist
				k = i
		Clusters[k].append(point)

	return Clusters

# Center of Gravity of a cluster is the average of coordinates of its members
def CenterOfGravity(Cluster) :
	size = len(Cluster)
	dim = len(Cluster[0])

	CoG = [0 for i in range(dim)]

	for i in range(dim) :
		for point in Cluster :
			CoG[i] += point[i]
		CoG[i] = CoG[i]/size

	return CoG	
		
	
# Lloyd algorithm to partition a given set of Points to N clusters.
# The algorithm iterates 2 steps:
#	1. Given Points and Centers, determine the Clusters.
#	2. For the Clusters, determine the new Centers
def Lloyd(Points, numClusters) :
	Centers = Points[:numClusters]
	dim = len(Points[0])

	converged = 0
	while not converged :
		converged = 1
		Clusters = Cluster(Points, Centers)
		newCenters = [CenterOfGravity(cluster) for cluster in Clusters]
		for i in range(numClusters) :
			for j in range(dim) :
				if abs(Centers[i][j] - newCenters[i][j]) > 0.0 :
					converged = 0
					Centers = newCenters
					break
			if not converged : break

	return Centers

# (Hard) Clustering assigns every point to one cluster.
# Soft Clustering assigns a point to multiple centers, each with a probablity.
# The probablity is called 'responsibility' and is computed as an exponential function of the distances.
def SoftClustering(Points, Centers, softClusters, Beta) :
	numClusters = len(Centers)
	numPoints = len(Points)

	for i in range(numPoints) :
		for j in range(numClusters) :
			dist = distance(Points[i], Centers[j])
			softClusters[i][j] = math.exp(-Beta * dist)
			# print "Point, Center, Distance and Force:", Points[i], Centers[j], dist, softClusters[i][j]

		total = sum(softClusters[i])
		for j in range(numClusters) :
			softClusters[i][j] = softClusters[i][j] / total
		# print softClusters[i]

	return softClusters

# Soft Centers are similar to (hard) centers, but we use the 'responsibility' weights to compute them.
def SoftCenters(Points, SoftClusters) :
	numPoints = len(Points)
	numCenters = len(SoftClusters[0])
	dim = len(Points[0])

	softCenters = []

	for j in range(numCenters) :

		totalWeight = 0.0
		for i in range(numPoints) :
			totalWeight += SoftClusters[i][j]

		center = [0 for k in range(dim)]
		for i in range(numPoints) :
			weight = SoftClusters[i][j] / totalWeight
			for k in range(dim) :
				center[k] += Points[i][k] * weight 
		softCenters.append(center)

	return softCenters

# Soft Lloyd is Lloyd for soft clustering	
def SoftLloyd(Points, numClusters, Beta) :
	numPoints = len(Points)
	Centers = Points[:numClusters]
	dim = len(Points[0])

	softClusters = [ [0 for i in range(numClusters)] for j in range(numPoints) ]

	converged = 0
	numIter = 0
	while not converged :
		numIter += 1
		if numIter > 100 : break
		converged = 1
		softClusters = SoftClustering(Points, Centers, softClusters, Beta)
		newCenters = SoftCenters(Points, softClusters) 
		# print "Iter: ", numIter, "Centers:", newCenters
		# print
		for i in range(numClusters) :
			for j in range(dim) :
				if abs(Centers[i][j] - newCenters[i][j]) > 1e-4 :
					converged = 0
					Centers = newCenters
					break
			if not converged : break

	return Centers

# Hierarchical clustering clusters 2 sub-clusters in each step.
# Distance between 2 subclusters is computed as the average of distances between 
# every point of one subcluster to every point in another subcluster.
def HierarchicalClustering(DistanceMatrix) :
	numPoints = len(DistanceMatrix)
	Clusters = [ [i] for i in range(1, numPoints+1) ]
	numClusters = numPoints

	while numClusters > 1 :
		# print DistanceMatrix
		# print Clusters

		smallest = sys.float_info.max
		for i in range(numClusters) :
			for j in range(i+1, numClusters) :
				if DistanceMatrix[i][j] < smallest :
					smallest = DistanceMatrix[i][j]
					I, J = i, j

		for k in Clusters[I] : print k,
		for k in Clusters[J] : print k,
		print

		for k in range(numClusters) :
			if k == I : continue
			if k == J : continue
			size1 = len(Clusters[I]) * len(Clusters[k])
			size2 = len(Clusters[J]) * len(Clusters[k])
			dist = DistanceMatrix[k][I] * size1
			dist += DistanceMatrix[k][J] * size2
			dist = dist/(size1 + size2)

			DistanceMatrix[k][I] = dist
			DistanceMatrix[I][k] = dist
			DistanceMatrix[k] = DistanceMatrix[k][:J] + DistanceMatrix[k][J+1:]
		DistanceMatrix[I] = DistanceMatrix[I][:J] + DistanceMatrix[I][J+1:]
		DistanceMatrix.remove(DistanceMatrix[J])

		Clusters[I] += Clusters[J]
		Clusters.remove(Clusters[J])
		numClusters -= 1

#-----------------------------------------------------------------
# Find centers based on the Farthest First Traversal heuristic 
def p6_1(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	numClusters, dimension = fp.readline().rstrip().split()
	numClusters = int(numClusters)

	Points = []
	while True :
		point = fp.readline()
		if not point : break

		point = point.rstrip().split()
		point = [float(coord) for coord in point]
		Points.append(point)

	Centers = FarthestFirstTraversal(Points, numClusters)

	sys.stdout = fpOut

	for center in Centers :
		for coord in center :
			print coord, 
		print


#-----------------------------------------------------------------
# Calculate distortion of a set of points w.r.to a given set of centers
def p6_2(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	numClusters, dimension = fp.readline().rstrip().split()
	numClusters = int(numClusters)

	Centers = []
	for i in range(numClusters) :
		center = fp.readline().rstrip().split()
		center = [float(coord) for coord in center]
		Centers.append(center)
	
	fp.readline()
	
	Points = []
	while True :
		point = fp.readline()
		if not point : break

		point = point.rstrip().split()
		point = [float(coord) for coord in point]
		Points.append(point)


	sys.stdout = fpOut

	print Distortion(Points, Centers)

#-----------------------------------------------------------------
# Lloyd algorithm for clustering 
def p6_3(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	numClusters, dim = fp.readline().rstrip().split()
	numClusters = int(numClusters)
	dim = int(dim)

	Points = []
	while True :
		point = fp.readline()
		if not point : break

		point = point.rstrip().split()
		point = [float(coord) for coord in point]
		Points.append(point)


	Centers = Lloyd(Points, numClusters)
	sys.stdout = fpOut
	for center in Centers :
		for i in range(dim) :
			print center[i],
		print	

#-----------------------------------------------------------------
# Lloyd algorithm for clustering 
def p6_4(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	numClusters, dim = fp.readline().rstrip().split()
	numClusters = int(numClusters)
	dim = int(dim)

	beta = float(fp.readline().rstrip())

	Points = []
	while True :
		point = fp.readline()
		if not point : break

		point = point.rstrip().split()
		point = [float(coord) for coord in point]
		Points.append(point)


	Centers = SoftLloyd(Points, numClusters, beta)
	sys.stdout = fpOut
	for center in Centers :
		for i in range(dim) :
			print center[i],
		print	

#-----------------------------------------------------------------
# Hierarchical Clustering
def p6_5(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	numPoints = int(fp.readline().rstrip())

	distMatrix = []
	while True :
		row = fp.readline()
		if not row : break

		row = row.rstrip().split()
		row = [float(cell) for cell in row]
		distMatrix.append(row)

	sys.stdout = fpOut

	HierarchicalClustering(distMatrix)

"""
p6_1('./Regression/6.1.1', './Solutions/6.1.1')
p6_1('./Regression/6.1.2', './Solutions/6.1.2')
p6_1('./Regression/6.1.3', './Solutions/6.1.3')
p6_2('./Regression/6.2.1', './Solutions/6.2.1')
p6_2('./Regression/6.2.2', './Solutions/6.2.2')
p6_2('./Regression/6.2.3', './Solutions/6.2.3')
p6_3('./Regression/6.3.1', './Solutions/6.3.1')
p6_3('./Regression/6.3.2', './Solutions/6.3.2')
p6_3('./Regression/6.3.3', './Solutions/6.3.3')
p6_4('./Regression/6.4.1', './Solutions/6.4.1')
p6_4('./Regression/6.4.2', './Solutions/6.4.2')
p6_4('./Regression/6.4.3', './Solutions/6.4.3')
p6_5('./Regression/6.5.1', './Solutions/6.5.1')
p6_5('./Regression/6.5.2', './Solutions/6.5.2')
"""
p6_5('./Regression/6.5.3', './Solutions/6.5.3')
