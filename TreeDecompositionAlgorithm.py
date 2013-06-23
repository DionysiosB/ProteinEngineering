import networkx as nx
import random
from operator import itemgetter
from numpy import *


#############################################################################################

def sort_by_degree(G):
    return sorted(G.degree(with_labels=True).items(),key = itemgetter(1))

##############################################################################################

def my_very_simple_dict_reverse_lookup(input_dictionary, input_value):
    for dict_index in input_dictionary:
	if input_dictionary[dict_index]==input_value:
	    return dict_index
	else:
	    print 'Did not find the requested value in the dictionary'
	    return -1


##############################################################################################

def my_very_simple_tuple_intersection(tuple1,tuple2):
    return tuple(set(tuple1)& set(tuple2))

##############################################################################################




def make_pairs(input_list):
    output_list=[]
    for k in range(len(input_list)):
	for m in range(k+1,len(input_list)):
	    output_list.append((input_list[k],input_list[m]))

    return output_list


######################################################################################

def dimacs2nx(filename):
    G = nx.Graph()
    for line in open(filename).readlines():
	l = line.split()
	if l[0]=='p':
	    N = int(l[2])
	    for n in range(N):
		G.add_node(n)
	if l[0]=='e':
	    G.add_edge(int(l[1]),int(l[2]))
	if l[0]=='c': continue
    return G


#######################################################################################
#######################################################################################

def tree_decomposition(input_graph):

    current_graph=input_graph.copy()
    decomposition_tree_vertices=list()
    counter=0;
    decomposition_tree=nx.Graph()
    tree_connectivity_dictionary=dict()
    for graph_vertex in current_graph.nodes():
	tree_connectivity_dictionary[graph_vertex]=[]


    while current_graph.order()>0:
	print current_graph.order()
	nodes_sorted_by_degree=sort_by_degree(current_graph)
	print 'nodes_sorted_by_degree', nodes_sorted_by_degree
	minimum_degree_vertex=nodes_sorted_by_degree[0][0]
	print 'Minimum Degree_vertex' , minimum_degree_vertex
	cliques_of_minimum_degree_vertex=nx.cliques_containing_node(current_graph,minimum_degree_vertex)
	print 'cliques_of_minimum_degree_vertex',cliques_of_minimum_degree_vertex
	number_of_cliques_containing_vertex=len(cliques_of_minimum_degree_vertex)
	print 'number_of_cliques_containing_vertex', number_of_cliques_containing_vertex
	minimum_degree_vertex_neighbors=current_graph.neighbors(minimum_degree_vertex)
	print 'minimum_degree_vertex_neighbors', minimum_degree_vertex_neighbors
	new_tree_vertex=[minimum_degree_vertex]
	print 'new_tree_vertex First element: ',new_tree_vertex
	new_tree_vertex.extend(minimum_degree_vertex_neighbors)
	new_tree_vertex=tuple(new_tree_vertex)
	decomposition_tree.add_node(new_tree_vertex)
	print 'decomposition_tree_vertices',decomposition_tree.nodes()
	if number_of_cliques_containing_vertex>1:
	    print 'Not Clique, will remove only one vertex'
	    pairs_of_neighbors=make_pairs(minimum_degree_vertex_neighbors)
	    print 'pairs_of_neighbors',pairs_of_neighbors
	    for additional_edge in pairs_of_neighbors:current_graph.add_edge(additional_edge[0],additional_edge[1])
	    toberemoved=[minimum_degree_vertex]
	    print 'toberemoved ', toberemoved
	else:
	    toberemoved=[minimum_degree_vertex]
	    print 'Clique detected, will try to remove more than one vertex'
	    number_of_clique_edges_per_vertex=len(minimum_degree_vertex_neighbors)
	    print 'number_of_clique_edges_per_vertex',number_of_clique_edges_per_vertex
	    print 'Checking all the vertex`s neighbors...'
	    print 'minimum_degree_vertex_neighbors', minimum_degree_vertex_neighbors
	    for temp_vertex in minimum_degree_vertex_neighbors:
		if current_graph.degree(temp_vertex)==number_of_clique_edges_per_vertex:
		    toberemoved.append(temp_vertex)
		    print 'Will ALSO remove vertex  ', temp_vertex
	for graph_vertex in new_tree_vertex:
	    if graph_vertex in toberemoved:
		current_graph.remove_node(graph_vertex)
		print 'Removed original graph vertex', graph_vertex
		tree_vertices_waiting=tree_connectivity_dictionary[graph_vertex]
		print 'For the removed node, tree_vertices_waiting: ' , tree_vertices_waiting
		for tree_vertex_waiting in tree_vertices_waiting:
		    print 'New Tree vertex:  ' , new_tree_vertex
		    print 'Tree Vertex waiting:', tree_vertex_waiting
		    decomposition_tree.add_edge(new_tree_vertex,tree_vertex_waiting)
		    print 'Connected tree vertices', new_tree_vertex, 'and   ' ,   tree_vertex_waiting
		    print 'The tree edges are now:    ', decomposition_tree.edges()
		    print 'THE NUMBER OF TREE EDGES ARE NOW:   ', len(decomposition_tree.edges())
		for tree_vertex_waiting in tree_vertices_waiting:
		    common_graph_nodes_between_tree_vertices=list(my_very_simple_tuple_intersection(new_tree_vertex,tree_vertex_waiting))
		    for graph_vertex in common_graph_nodes_between_tree_vertices:
			tree_connectivity_dictionary[graph_vertex].remove(tree_vertex_waiting)
			print 'Removed from dictionary entry', graph_vertex , 'tree node ', tree_vertex_waiting
			print 'Now the new dictionary is:  ' , tree_connectivity_dictionary


	    else:
		tree_connectivity_dictionary[graph_vertex].append(new_tree_vertex)
		print 'New tree_connectivity_dictionary node appended. New tree_connectivity_dictionary ', tree_connectivity_dictionary
	print 'tree_connectivity_dictionary:  '	 , tree_connectivity_dictionary
	print 'decomposition_tree.nodes:     ', decomposition_tree.nodes()
	print 'decomposition_tree.edges:     ', decomposition_tree.edges()



    return decomposition_tree








##########
##############



G=nx.Graph()
G.add_node('a');G.add_node('b');G.add_node('c');G.add_node('d');G.add_node('e');G.add_node('f');
G.add_node('g');G.add_node('h');G.add_node('i');G.add_node('j');G.add_node('k');G.add_node('l');G.add_node('m');

G.add_edge('a','b');G.add_edge('a','c');G.add_edge('b','d');G.add_edge('c','d');G.add_edge('c','e');G.add_edge('c','k');
G.add_edge('c','l');G.add_edge('c','m');G.add_edge('d','f');G.add_edge('d','m');G.add_edge('e','f');G.add_edge('e','i');
G.add_edge('e','j');G.add_edge('e','m');G.add_edge('f','g');G.add_edge('f','h');G.add_edge('f','m');G.add_edge('i','j');
G.add_edge('k','l');



aaa=tree_decomposition(G);
aaa.nodes()
aaa.edges()







