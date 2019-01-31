#!flask/bin/python
from flask import Flask
from flask import Flask, jsonify, request
from flask_cors import CORS
import os
import numpy as np
import Levenshtein
import time
import sys
from decimal import Decimal



dir_path = os.path.dirname(os.path.realpath(__file__))


app = Flask(__name__)
CORS(app)

asd = {
		"x": ["a", "z"],
		"a": ["x", "z", "y"],
		"z": ["a", "x"],
		"y": ["a"]
}


def generate_edges(graph):
	edges = []
	for node in graph:
			for neighbour in graph[node]:
					edges.append((node, neighbour))
	return edges

def generate_edges_from_nodes(graph, nodes):
	edges = []
	for node in nodes:
			for neighbour in graph[node]:
					edges.append((node, neighbour))
	return edges

def filter_edges(edges):
	filteredEdges = []
	for (x,y) in edges:
		temp= ''
		temp = x + y
		backtemp = temp[::-1]
		if ((temp not in filteredEdges) and (backtemp not in filteredEdges)):
			filteredEdges.append(str(temp))
	return filteredEdges



def find_isolated_nodes(graph):
	""" returns a list of isolated nodes. """
	isolated = []
	for node in graph:
			if not graph[node]:
					isolated += node
	return isolated

def get_vertices(graph):
	""" returns the vertices of a graph """
	return list(graph.keys())

def get_degree_sequence(graph):
	""" calculates the degree sequence """
	seq = []
	for vertex in graph:
			seq.append(vertex_degree(graph,vertex))
	seq.sort(reverse=True)
	return tuple(seq)

def min_delta(graph):
	""" the minimum degree of the vertices """
	min = 100000000
	for vertex in graph:
			min_vertex_degree = vertex_degree(graph, vertex)
			if min_vertex_degree < min:
					min = min_vertex_degree
	return min
				
def max_Delta(graph):
	""" the maximum degree of the vertices """
	max = 0
	for vertex in graph:
			 max_vertex_degree = vertex_degree(graph, vertex)
			 if max_vertex_degree > max:
				 max = max_vertex_degree
	return max

def vertex_degree(graph, vertex):
	""" The degree of a vertex is the number of edges connecting
			it, i.e. the number of adjacent vertices. Loops are counted 
			double, i.e. every occurence of vertex in the list 
			of adjacent vertices. """ 
	adj_vertices =  graph[vertex]
	degree = len(adj_vertices) + adj_vertices.count(vertex)
	return degree

def get_density(graph):
	""" method to calculate the density of a graph """
	g = graph
	V = len(g.keys())
	E = len(generate_edges(g))
	return 2.0 * E / (V *(V - 1))

def find_all_paths(graph, start_vertex, end_vertex, path=[]):
	""" find all paths from start_vertex to 
			end_vertex in graph """
	graph = graph 
	path = path + [start_vertex]
	if start_vertex == end_vertex:
			return [path]
	if start_vertex not in graph:
			return []
	paths = []
	for vertex in graph[start_vertex]:
			if vertex not in path:
					extended_paths = find_all_paths(graph,vertex,end_vertex,path)
					for p in extended_paths: 
							paths.append(p)
	return paths

def get_diameter(graph):
	""" calculates the diameter of the graph """
	graph = graph
	v = get_vertices(graph)
	pairs = [ (v[i],v[j]) for i in range(len(v)) for j in range(i+1, len(v)-1)]
	smallest_paths = []
	for (s,e) in pairs:
			paths = find_all_paths(graph,s,e)
			smallest = sorted(paths, key=len)[0]
			smallest_paths.append(smallest)

	smallest_paths.sort(key=len)

	# longest path is at the end of list, 
	# i.e. diameter corresponds to the length of this path
	diameter = len(smallest_paths[-1]) - 1
	return diameter

def levenshteinDist(seq1, seq2):
	size_x = len(seq1) + 1
	size_y = len(seq2) + 1
	matrix = np.zeros((size_x, size_y))
	for x in xrange(size_x):
			matrix[x, 0] = x
	for y in xrange(size_y):
			matrix[0, y] = y

	for x in xrange(1, size_x):
			for y in xrange(1, size_y):
					if seq1[x-1] == seq2[y-1]:
							matrix[x, y] = min(
									matrix[x-1, y] + 1,
									matrix[x-1, y-1],
									matrix[x, y-1] + 1
							)
					else:
							matrix[x, y] = min(
									matrix[x-1, y] + 1,
									matrix[x-1, y-1] + 1,
									matrix[x, y-1] + 1
							)
	print (matrix)
	return (matrix[size_x - 1, size_y - 1])

@app.route('/')
def index():
	return "Hello, World!"


@app.route('/graph-distance/api/v1.0/tasks', methods=['GET'])
def get_tasks():
	return jsonify({'graph': generate_edges(asd)})


@app.route('/graph-distance/api/v1.0/graph', methods=['POST'])
def create_task():
	print '###############'
	graph = request.json['sourceobj']
	nodes = request.json['nodes']
	print "NODES: ", nodes
	edges = generate_edges(graph)
	isolated = find_isolated_nodes(graph)
	vertices = get_vertices(graph)
	min_degree = min_delta(graph)
	max_degree = max_Delta(graph)
	degree_sequence = get_degree_sequence(graph)
	density = get_density(graph)
	diameter = get_diameter(graph)
	print({'graph': graph, 'edges': edges, 'isolated': isolated, 'vertices': vertices, 'degree_sequence': degree_sequence, 'min_degree': min_degree, 'max_degree': max_degree, 'density': density, 'diameter':diameter,'nodes':nodes})
	print '###############'
	return jsonify({'graph': graph, 'edges': edges, 'isolated': isolated, 'vertices': vertices, 'degree_sequence': degree_sequence, 'min_degree': min_degree, 'max_degree': max_degree,'density': density, 'diameter':diameter,'nodes':nodes}), 201


@app.route('/graph-distance/api/v1.0/calculate_graphs_score', methods=['POST'])
def similarity_score():
	print "****************************"
	sourcegraph = request.json['sourcegraph']
	sourcenode = request.json['sourcenode']
	matchgraph = request.json['matchgraph']
	matchnode = request.json['matchnode']
	sourceEdge = generate_edges_from_nodes(sourcegraph, sourcenode)
	matchEdge = generate_edges_from_nodes(matchgraph, matchnode)
	sourceedge_filtered = filter_edges(sourceEdge)
	matchedge_filtered = filter_edges(matchEdge)
	sourceedge_filtered_new = sourceedge_filtered[:]
	SEDStartTimer = time.time()
	roundFigure = 4
	for se in sourceedge_filtered_new:
		es = se[::-1]
		if (se in matchedge_filtered):
			matchedge_filtered.remove(se)
			sourceedge_filtered.remove(se)
		elif (es in matchedge_filtered):
			matchedge_filtered.remove(es)
			sourceedge_filtered.remove(se)
		else:
			continue
	print '*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*'
	if (len(sourceedge_filtered) == 0 and len(matchedge_filtered) == 0):
		print "SAME: 1"
		SEDStopTimer = time.time()
		SEDtotalTimer = SEDStartTimer - SEDStopTimer
		print "TIME1: ", SEDtotalTimer
		print "TIME1: ", '{0:.{1}f}'.format(SEDtotalTimer, roundFigure)
		return jsonify({'result': 1, 'sedtime': '{0:.{1}f}'.format(SEDtotalTimer, roundFigure)}), 201
	elif ((len(matchedge_filtered) != 0) and (len(matchedge_filtered) >= len(sourceedge_filtered))):
		asd = []
		scu = ''
		mct = ''
		for medg in matchedge_filtered:
			mct += medg
			# print "EDG1: ", medg, "MCT: ", mct
		for sedg in sourceedge_filtered:
			scu += sedg
			# print "EDG2: ", sedg, "MCT: ", scu
		# for idx,li in enumerate(sourceedge_filtered):
			# score = Levenshtein.ratio(li, matchedge_filtered[idx])
			# asd.append({'sourceedge_filtered': li, 'matchedge_filtered': matchedge_filtered[idx], 'score': score })
		# print "sourceedge_filtered 2: ", asd
		print "SOURCE: ", scu, "MATCH: ", 
		score = Levenshtein.ratio(scu, mct)
		print "score 2: ", score
		SEDStopTimer = time.time()
		SEDtotalTimer = SEDStartTimer - SEDStopTimer
		print "TIME2: ", SEDtotalTimer
		return jsonify({'result': asd, 'sedtime': SEDtotalTimer}), 201
	elif ((len(sourceedge_filtered) != 0) and (len(sourceedge_filtered) >= len(matchedge_filtered))):
		asd = []
		for idx,li in enumerate(matchedge_filtered):
			score = Levenshtein.ratio(li, sourceedge_filtered[idx])
			asd.append({'matchedge_filtered': li, 'sourceedge_filtered': sourceedge_filtered[idx], 'score': score })
		print "matchedge_filtered 3: ", asd
		SEDStopTimer = time.time()
		SEDtotalTimer = SEDStartTimer - SEDStopTimer
		print "TIME3: ", SEDtotalTimer
		return jsonify({'result': asd, 'sedtime': SEDtotalTimer}), 201
	else:
		for idx,li in enumerate(sourceedge_filtered):
			asd = []
			score = Levenshtein.ratio(li, sourceedge_filtered[idx])
			asd.append({'matchedge_filtered': li, 'sourceedge_filtered': sourceedge_filtered[idx], 'score': score })
		print "sourceedge_filtered else: ", asd
		SEDStopTimer = time.time()
		SEDtotalTimer = SEDStartTimer - SEDStopTimer
		print "TIME4: ", SEDtotalTimer
		return jsonify({'result': asd, 'sedtime': SEDtotalTimer}), 201
	print '*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*'
	print "****************************"


if __name__ == '__main__':
		app.run(debug=True)