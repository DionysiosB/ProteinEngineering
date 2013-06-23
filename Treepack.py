import networkx as nx
import random
from operator import itemgetter
from numpy import *
import sys




######################################################################################################################################################
######################################################################################################################################################

def sort_by_degree(G):
    return sorted(G.degree(with_labels=True).items(),key = itemgetter(1))

######################################################################################################################################################
######################################################################################################################################################

def my_very_simple_dict_reverse_lookup(input_dictionary, input_value):
    for dict_index in input_dictionary:
	if input_dictionary[dict_index]==input_value:
	    return dict_index
	else:
	    print 'Did not find the requested value in the dictionary'
	    return -1


######################################################################################################################################################
######################################################################################################################################################

def my_very_simple_tuple_intersection(tuple1,tuple2):
    return tuple(set(tuple1)& set(tuple2))

######################################################################################################################################################
######################################################################################################################################################




def make_pairs(input_list):
    output_list=[]
    for k in range(len(input_list)):
	for m in range(k+1,len(input_list)):
	    output_list.append((input_list[k],input_list[m]))

    return output_list


######################################################################################################################################################
######################################################################################################################################################

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

######################################################################################################################################################
######################################################################################################################################################
def tree_decomposition(input_graph):

    current_graph=input_graph.copy()
    decomposition_tree_vertices=list()
    counter=0;
    decomposition_tree=nx.Graph()
    tree_connectivity_dictionary=dict()
    for graph_vertex in current_graph.nodes():
	tree_connectivity_dictionary[graph_vertex]=[]


    while current_graph.order()>0:
	nodes_sorted_by_degree=sort_by_degree(current_graph)
	minimum_degree_vertex=nodes_sorted_by_degree[0][0]
	cliques_of_minimum_degree_vertex=nx.cliques_containing_node(current_graph,minimum_degree_vertex)
	number_of_cliques_containing_vertex=len(cliques_of_minimum_degree_vertex)
	minimum_degree_vertex_neighbors=current_graph.neighbors(minimum_degree_vertex)
	new_tree_vertex=[minimum_degree_vertex]
	new_tree_vertex.extend(minimum_degree_vertex_neighbors)
	new_tree_vertex=tuple(new_tree_vertex)
	decomposition_tree.add_node(new_tree_vertex)
	if number_of_cliques_containing_vertex>1:
	    pairs_of_neighbors=make_pairs(minimum_degree_vertex_neighbors)
	    for additional_edge in pairs_of_neighbors:current_graph.add_edge(additional_edge[0],additional_edge[1])
	    toberemoved=[minimum_degree_vertex]
	else:
	    toberemoved=[minimum_degree_vertex]
	    number_of_clique_edges_per_vertex=len(minimum_degree_vertex_neighbors)
	    for temp_vertex in minimum_degree_vertex_neighbors:
		if current_graph.degree(temp_vertex)==number_of_clique_edges_per_vertex:
		    toberemoved.append(temp_vertex)
	for graph_vertex in new_tree_vertex:
	    if graph_vertex in toberemoved:
		current_graph.delete_node(graph_vertex)
		tree_vertices_waiting=tree_connectivity_dictionary[graph_vertex]
		for tree_vertex_waiting in tree_vertices_waiting:
		    decomposition_tree.add_edge(new_tree_vertex,tree_vertex_waiting)
		for tree_vertex_waiting in tree_vertices_waiting:
		    common_graph_nodes_between_tree_vertices=list(my_very_simple_tuple_intersection(new_tree_vertex,tree_vertex_waiting))
		    for graph_vertex in common_graph_nodes_between_tree_vertices:
			tree_connectivity_dictionary[graph_vertex].remove(tree_vertex_waiting)


	    else:
		tree_connectivity_dictionary[graph_vertex].append(new_tree_vertex)

    return decomposition_tree


######################################################################################################################################################
######################################################################################################################################################

def find_tree_leaves(nx_tree_input):

    tree_leaves=list()
    for tree_vertex in nx_tree_input.nodes():
	if nx_tree_input.degree(tree_vertex)==1:tree_leaves.append(tree_vertex)

    return tree_leaves

######################################################################################################################################################
######################################################################################################################################################

def find_optimal_tree_root(nx_tree_input):

    tree_root=nx.center(nx_tree_input)
    return tree_root[0]

######################################################################################################################################################
######################################################################################################################################################

def find_combinations_list(input_dict_of_lists):

    number_of_sets=len(input_dict_of_lists)
    cardinality_dict=dict()
    for k in input_dict_of_lists:
	cardinality_dict[k]=len(input_dict_of_lists[k])
	if cardinality_dict[k]<=0:
	    print 'The elements of the list must be strictly positive integers. Exiting....'
	    return -1

    repetition_dict=dict()
    temp_repetition=0
    for m in cardinality_dict:
	if temp_repetition==0:
	    repetition_dict[m]=1
	    temp_repetition=cardinality_dict[m]
	else:
	    repetition_dict[m]=temp_repetition
	    temp_repetition=temp_repetition*cardinality_dict[m]

    total_number_of_combinations=temp_repetition

    output_combination_list=list()

    for combination_number in range(total_number_of_combinations):
	current_combination_list=list()
	for current_set in input_dict_of_lists:
	    current_combination_list.append( input_dict_of_lists[current_set][(combination_number/repetition_dict[current_set])%cardinality_dict[current_set]])

	output_combination_list.append(current_combination_list)


    return output_combination_list


######################################################################################################################################################
######################################################################################################################################################


def find_tree_structure(nx_tree_input):

    tree_root=find_optimal_tree_root(nx_tree_input)
    tree_leaves=find_tree_leaves(nx_tree_input)

    tree_structure_children_to_parent=dict()
    tree_structure_parent_to_children=dict()

    for current_leaf in tree_leaves:
	current_path=nx.shortest_path(nx_tree_input,tree_root,current_leaf)
	current_path_length=len(current_path)
	for m in range(1,current_path_length):
	    tree_structure_children_to_parent[current_path[m]]=current_path[m-1]

	    if current_path[m-1] not in tree_structure_parent_to_children:tree_structure_parent_to_children[current_path[m-1]]=[current_path[m]]
	    elif current_path[m] not in tree_structure_parent_to_children[current_path[m-1]]:tree_structure_parent_to_children[current_path[m-1]].append(current_path[m])
	    else: continue


    return [tree_structure_children_to_parent,tree_structure_parent_to_children]



######################################################################################################################################################
######################################################################################################################################################


def Dynamic_Programming_for_decomposed_trees(input_tree,input_dictionary,interaction_dictionary):     #Input dictionary= The alternative rotamers for each residue

    current_tree=input_tree.copy()
    master_dictionary=dict()
    for dummy in current_tree.nodes():
	master_dictionary[dummy]=dict()


    tree_root=find_optimal_tree_root(input_tree)
    next_tree_leaves=find_tree_leaves(current_tree)
    current_tree_leaves=find_tree_leaves(current_tree)

    [tree_structure_children_to_parent,tree_structure_parent_to_children]=find_tree_structure(input_tree)


    while len(current_tree_leaves)>0:

	current_tree_leaves=next_tree_leaves[:]
	next_tree_leaves=list()
	if tree_root in current_tree_leaves: current_tree_leaves.remove(tree_root)     #The root HAS to be computed after ALL the other nodes are computed

	for current_node in current_tree_leaves:
	    parent_dict=dict()
	    children_dict=dict()

	    if current_node in tree_structure_parent_to_children:
		for child in tree_structure_parent_to_children[current_node]:
		    children_dict[child]=master_dictionary[child]

	    parent_of_node=tree_structure_children_to_parent[current_node]
	    if parent_of_node not in next_tree_leaves:next_tree_leaves.append(parent_of_node)

	    master_dictionary[current_node]=find_optimal_combination(input_dictionary,interaction_dictionary,current_node,parent_of_node,children_dict)


	#Now, once we are done with all the other nodes, we move on to the tree root

    root_node=tree_root
    parent_of_root=-1
    children_dict=dict()

    for root_child in tree_structure_parent_to_children[root_node]:
	children_dict[root_child]=master_dictionary[root_child]

    master_dictionary[root_node]=find_optimal_combination(input_dictionary,interaction_dictionary,root_node,parent_of_root,children_dict)

    final_dictionary=master_dictionary[root_node]    #This is a dictionary of the form: set:value

    best_combination=final_dictionary.keys()[0]
    minimum_value=final_dictionary[best_combination]

    return [best_combination, minimum_value]

######################################################################################################################################################
######################################################################################################################################################
def find_optimal_combination(input_dictionary,interaction_dictionary,current_node,parent_of_node,children_dict):


    if parent_of_node != -1:
	node_with_parent_intersection=tuple( set(current_node) & set(parent_of_node)  )
	node_not_parent_elements=tuple( set( current_node) - set(parent_of_node))
    else:
	node_with_parent_intersection=tuple()
	node_not_parent_elements=current_node


    if len(children_dict)>0:
	leaf_indicator=0
	children_of_current_node=children_dict.keys()
    else:
	leaf_indicator=1

    iterator_dictionary=dict()
    variable_dictionary=dict()
    output_dictionary=dict()


    if parent_of_node != -1:
	for iterator in node_with_parent_intersection:
	    iterator_dictionary[iterator]=input_dictionary[iterator]

    for variable in node_not_parent_elements:
	variable_dictionary[variable]=input_dictionary[variable]

    all_iterator_combinations=find_combinations_list(iterator_dictionary)
    all_variable_combinations=find_combinations_list(variable_dictionary)


    if len(all_iterator_combinations)>0:
	for current_iterator_combination in all_iterator_combinations:

	    optimal_variable_combination=list()
	    smallest_value=sys.maxint

	    for current_variable_combination in all_variable_combinations:

		current_node_interactions_value=find_total_combination_value(interaction_dictionary,current_iterator_combination,current_variable_combination)

		integrated_value=current_node_interactions_value

		if leaf_indicator==0:
		    provided_set=(set(current_iterator_combination) | set(current_variable_combination)  )
		    integrated_set=provided_set
		    for current_child in children_of_current_node:
			for dummytuple in children_dict[current_child]:
			    dummyset=set(dummytuple)
			    if dummyset.issubset(provided_set):
				integrated_set= ( set(children_dict[current_child][dummytuple][0]) | integrated_set)
				integrated_value+=children_dict[current_child][dummytuple][1]
		else:
		    provided_set=(set(current_iterator_combination) | set(current_variable_combination)	     )
		    integrated_set=provided_set
		if integrated_value < smallest_value:
		    smallest_value=integrated_value
		    optimal_integrated_combination=integrated_set

	    output_dictionary[tuple(current_iterator_combination)]=[tuple(optimal_integrated_combination), smallest_value]

    else:
	current_iterator_combination=[]
	optimal_variable_combination=list()
	smallest_value=sys.maxint

	for current_variable_combination in all_variable_combinations:
	    current_node_interactions_value=find_total_combination_value(interaction_dictionary,current_iterator_combination,current_variable_combination)
	    integrated_value=current_node_interactions_value
	    provided_set=set(current_variable_combination)
	    integrated_set=provided_set

	    if leaf_indicator==0:
		for current_child in children_of_current_node:
		    for dummytuple in children_dict[current_child]:
			dummyset=set(dummytuple)
			if dummyset.issubset(provided_set):
			    integrated_set= ( set(children_dict[current_child][dummytuple][0]) | integrated_set)
			    integrated_value+=children_dict[current_child][dummytuple][1]
	    else:
		provided_set=(set(current_iterator_combination) | set(current_variable_combination)	     )
		integrated_set=provided_set

	    if integrated_value < smallest_value:
		smallest_value=integrated_value
		optimal_integrated_combination=integrated_set

	output_dictionary[tuple(optimal_integrated_combination)]=smallest_value
	
	print 'Optimal combination:', tuple(optimal_integrated_combination)
	print 'Minimum Energy: ', smallest_value


    return output_dictionary



######################################################################################################################################################
######################################################################################################################################################

def find_total_combination_value(interaction_dictionary,list1, list2):

    #Check for common elements in the list

    if len( set(list1) & set(list2))>0:
	print 'There are common elements in the two lists... This is not permitted. Returning -1'
	return -1

    total_list=list1[:]
    total_list.extend(list2)
    number_of_elements=len(total_list)

    output=0

    for k in range(number_of_elements):
	for m in range(k,number_of_elements):
	    if tuple([total_list[k],total_list[m]]) in interaction_dictionary:
		output+=interaction_dictionary[tuple([total_list[k],total_list[m]])];

    return output





######################################################################################################################################################
######################################################################################################################################################

dict1={'a':[6,18], 'b':[19,4], 'c':[7,17],'d':[5,20],'e':[8,26],'f':[16,3],'g':[21,9],'h':[15,10], 'i':[14,2], 'j':[11,23],'k':[22,12], 'l':[13,24],'m':[25,1],  }
dict2=dict()
for k in range(45):
    for m in range(45):
	dict2[tuple([k,m])]=abs(k-m)



G=nx.Graph()
G.add_node('a');G.add_node('b');G.add_node('c');G.add_node('d');G.add_node('e');G.add_node('f');
G.add_node('g');G.add_node('h');G.add_node('i');G.add_node('j');G.add_node('k');G.add_node('l');G.add_node('m');

G.add_edge('a','b');G.add_edge('a','c');G.add_edge('b','d');G.add_edge('c','d');G.add_edge('c','e');G.add_edge('c','k');
G.add_edge('c','l');G.add_edge('c','m');G.add_edge('d','f');G.add_edge('d','m');G.add_edge('e','f');G.add_edge('e','i');
G.add_edge('e','j');G.add_edge('e','m');G.add_edge('f','g');G.add_edge('f','h');G.add_edge('f','m');G.add_edge('i','j');
G.add_edge('k','l');

test=tree_decomposition(G)
test.nodes()
test.edges()

Dynamic_Programming_for_decomposed_trees(test,dict1,dict2)



