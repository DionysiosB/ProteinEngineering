#!/usr/bin/env python
import networkx as nx
import random
from operator import itemgetter
from numpy import *
import sys
#import matplotlib.pyplot as plt
import pickle
import random

############################################################################################################################
#############################	     TREEPACK ALGORITHM	       #############################################################
#############################	    DIONYSIOS BARMPOUTIS       #############################################################
############################################################################################################################


######################################################################################################################################################
######################################################################################################################################################

def interaction_pair_graph_builder(input_dictionary, interaction_dictionary, cutoff=0.001):

    interactions_graph=nx.Graph()

    for interaction_pair in interaction_dictionary:
	if abs(interaction_dictionary[interaction_pair]) < cutoff: continue
	first_element=dict_of_disjoint_lists_reverse_lookup(input_dictionary,interaction_pair[0])
	second_element=dict_of_disjoint_lists_reverse_lookup(input_dictionary,interaction_pair[1])
	interactions_graph.add_edge(first_element,second_element)

    if len(nx.connected_components(interactions_graph))>1:
	print 'WARNING WARNING WARNING Not all nodes in the input dictionary have interactions. Returning BIGGEST COMPONENT ONLY'
	graph_components_list=nx.connected_component_subgraphs(interactions_graph)
	biggest_component=graph_components_list[0]
	for current_component in graph_components_list:
	    if current_component.order()>biggest_component.order():
		biggest_component=current_component
	print 'BIGGEST COMPONENT', biggest_component.nodes()
	return biggest_component

    else:
	return interactions_graph


######################################################################################################################################################
######################################################################################################################################################

def sort_by_degree(G):
    return sorted(G.degree(with_labels=True).items(),key = itemgetter(1))

######################################################################################################################################################
######################################################################################################################################################

def dict_of_disjoint_lists_reverse_lookup(input_dictionary, input_value):
    for dict_index in input_dictionary:
	if input_value in input_dictionary[dict_index]:
	    return dict_index

    print 'Did not find the requested value in the dictionary'
    return -1


######################################################################################################################################################
######################################################################################################################################################

def dict_remove_small_values(input_dict, threshold):

    output_dict=input_dict.copy()

    for dummy_key in input_dict:
	if abs(input_dict[dummy_key])<abs(threshold):
	    del(output_dict[dummy_key])

    return output_dict


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
    current_graph.remove_edges_from(current_graph.selfloop_edges())

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
		current_graph.remove_node(graph_vertex)
		tree_vertices_waiting=tree_connectivity_dictionary[graph_vertex]
		for tree_vertex_waiting in tree_vertices_waiting:
		    decomposition_tree.add_edge(new_tree_vertex,tree_vertex_waiting)


		temp_copy_tree_vertices_waiting=tree_vertices_waiting[:]
		for tree_vertex_waiting in temp_copy_tree_vertices_waiting:
		    common_graph_nodes_between_tree_vertices=my_very_simple_tuple_intersection(new_tree_vertex,tree_vertex_waiting)
		    for candidate in common_graph_nodes_between_tree_vertices:
			if tree_vertex_waiting in tree_connectivity_dictionary[candidate]:tree_connectivity_dictionary[candidate].remove(tree_vertex_waiting)

		del tree_connectivity_dictionary[graph_vertex]
	    else:
		tree_connectivity_dictionary[graph_vertex].append(new_tree_vertex)


    if ((decomposition_tree.number_of_nodes()-decomposition_tree.number_of_edges()) < 1):
	print 'WARNING WARNING WARNING: THE OUTPUT GRAPH IS ****NOT**** A TREE, IT INCLUDES CYCLES'
    elif ((decomposition_tree.number_of_nodes()-decomposition_tree.number_of_edges()) > 1):
	print 'WARNING WARNING WARNING: THE OUTPUT GRAPH IS ****NOT**** A TREE, IT IS DISCONNECTED'
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
		    provided_set=(set(current_iterator_combination) | set(current_variable_combination)	 )
		    integrated_set=provided_set
		    for current_child in children_of_current_node:
			for dummytuple in children_dict[current_child]:
			    dummyset=set(dummytuple)
			    if dummyset.issubset(provided_set):
				integrated_set= ( set(children_dict[current_child][dummytuple][0]) | integrated_set)
				integrated_value=integrated_value+children_dict[current_child][dummytuple][1]

		else:
		    provided_set=(set(current_iterator_combination) | set(current_variable_combination)	    )
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
		provided_set=(set(current_iterator_combination) | set(current_variable_combination))
		integrated_set=provided_set
	    if integrated_value < smallest_value:
		smallest_value=integrated_value
		optimal_integrated_combination=integrated_set

	output_dictionary[tuple(optimal_integrated_combination)]=smallest_value

    return output_dictionary



######################################################################################################################################################
######################################################################################################################################################

def find_combinations_list(input_dict_of_lists):

    number_of_sets=len(input_dict_of_lists)
    cardinality_dict=dict()
    for k in input_dict_of_lists:
	cardinality_dict[k]=len(input_dict_of_lists[k])
	if cardinality_dict[k]<=0:
	    print 'Each list in the dictionary must include at least one element. Exiting....'
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

def find_total_combination_value(interaction_dictionary,iterator_list, variable_list):

    #Check for common elements in the list
    if len( set(iterator_list) & set(variable_list))>0:
	print 'There are common elements in the two lists... This is not permitted. Returning -1'
	return -1

    total_list=iterator_list[:]
    total_list.extend(variable_list)

    iterator_list_elements=len(iterator_list)
    number_of_elements=len(total_list)

    output=0
    for a in range(iterator_list_elements):
	for b in range(a,iterator_list_elements):
	    if tuple([iterator_list[a],iterator_list[b]]) in interaction_dictionary:
		output-=interaction_dictionary[tuple([iterator_list[a],iterator_list[b]])]
	    elif tuple([iterator_list[b],iterator_list[a]]) in interaction_dictionary:
		output-=interaction_dictionary[tuple([iterator_list[b],iterator_list[a]])]

    for k in range(number_of_elements):
	for m in range(k,number_of_elements):
	    if tuple([total_list[k],total_list[m]]) in interaction_dictionary:
		output+=interaction_dictionary[tuple([total_list[k],total_list[m]])]
	    elif tuple([total_list[m],total_list[k]]) in interaction_dictionary:
		output+=interaction_dictionary[tuple([total_list[m],total_list[k]])]


    return output

######################################################################################################################################################
######################################################################################################################################################
def pick_random_initial_variable_combination(variable_dictionary):

    variable_initial_list=list()
    for current_variable in variable_dictionary:
	dummy_list=variable_dictionary[current_variable]
	try:variable_initial_list.append(dummy_list[0])
	except: print 'One or more entries in the given dictionary has no values. Returning -1....';return -1

    return variable_initial_list

######################################################################################################################################################
######################################################################################################################################################
def change_variable_rotamer(interaction_dictionary,variable_dictionary,iterator_rotamers,variable_rotamers,current_value):

    if len(variable_dictionary)==0: print 1/0
    target_residue=random.choice(variable_dictionary.keys())
    target_residue_rotamers=variable_dictionary[target_residue][:]
    old_variable_rotamer=list(set(variable_rotamers) & set(target_residue_rotamers))

    if len(old_variable_rotamer)==0:
	print 'WARNING:The provided residue has no rotamers in the given variable list! Returning -1....'
	return -1[tuple([new_variable_rotamer_list]),new_value]
    elif len(old_variable_rotamer)>1:
	print 'WARNING:The provided variable list has more than one rotamers from the same residue! Returning -1...'
	return -1

    old_variable_rotamer=old_variable_rotamer[0]
    #target_residue_rotamers.remove(old_variable_rotamer)
    new_variable_rotamer=random.choice(target_residue_rotamers)

    new_variable_rotamers=variable_rotamers[:]
    new_variable_rotamers.remove(old_variable_rotamer)
    new_variable_rotamers.append(new_variable_rotamer)

    all_other_rotamers=variable_rotamers[:]
    all_other_rotamers.remove(old_variable_rotamer)
    all_other_rotamers.extend(iterator_rotamers)

    new_value=current_value

    for current_rotamer in all_other_rotamers:
	if tuple([current_rotamer,old_variable_rotamer]) in interaction_dictionary:
	    new_value-=interaction_dictionary[tuple([current_rotamer,old_variable_rotamer])]
	elif tuple([old_variable_rotamer,current_rotamer]) in interaction_dictionary:
	    new_value-=interaction_dictionary[tuple([old_variable_rotamer,current_rotamer])]
	if tuple([current_rotamer,new_variable_rotamer]) in interaction_dictionary:
	    new_value+=interaction_dictionary[tuple([current_rotamer,new_variable_rotamer])]
	elif tuple([new_variable_rotamer,current_rotamer]) in interaction_dictionary:
	    new_value+=interaction_dictionary[tuple([new_variable_rotamer,current_rotamer])]

    print 'tuple([new_variable_rotamers,new_value])', tuple([new_variable_rotamers,new_value])
    return tuple([new_variable_rotamers,new_value])

######################################################################################################################################################
######################################################################################################################################################
def Dynamic_Programming_for_decomposed_trees(input_tree,input_dictionary,interaction_dictionary):

    #Input dictionary= The alternative rotamers for each residue

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

	    master_dictionary[current_node]=find_stochastic_optimal_combination(input_dictionary,interaction_dictionary,current_node,parent_of_node,children_dict)


    #Now, once we are done with all the other nodes, we move on to the tree root

    root_node=tree_root
    parent_of_root=-1
    children_dict=dict()

    for root_child in tree_structure_parent_to_children[root_node]:
	children_dict[root_child]=master_dictionary[root_child]


    master_dictionary[root_node]=find_optimal_combination(input_dictionary,interaction_dictionary,root_node,parent_of_root,children_dict)
    #master_dictionary[root_node]=find_stochastic_optimal_combination(input_dictionary,interaction_dictionary,root_node,parent_of_root,children_dict)
    final_dictionary=master_dictionary[root_node]    #This is a dictionary of the form: set:value

    best_combination=final_dictionary.keys()[0]
    minimum_value=final_dictionary[best_combination]

    print 'Best combination:', best_combination, '    Minimum Value:', minimum_value
    return [best_combination, minimum_value]

######################################################################################################################################################
######################################################################################################################################################
def find_stochastic_optimal_combination(input_dictionary,interaction_dictionary,current_node,parent_of_node,children_dict):

    print '-----------------------------------------------------------------------------------------------------------------------'
    print 'CURRENT NODE:', current_node

    number_of_trials=2

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
    print 'ALL iterator combinations:  ', all_iterator_combinations


    if len(all_iterator_combinations)>0:
	for current_iterator_combination in all_iterator_combinations:

	    print '================================================'
	    print 'current_iterator_combination	 ',current_iterator_combination


	    current_variable_combination=pick_random_initial_variable_combination(variable_dictionary)
	    current_interactions_value=find_total_combination_value(interaction_dictionary,current_iterator_combination,current_variable_combination)
	    optimal_variable_combination=list()
	    smallest_value=sys.maxint

	    for trial_number in range(number_of_trials):

		print 'Trial number ----> ', trial_number
		old_variable_rotamers=current_variable_combination
		print 'Current Variable rotamers , ', current_variable_combination
		old_value=current_interactions_value
		print 'old_value, ', old_value
		new_attempt=change_variable_rotamer(interaction_dictionary,variable_dictionary,current_iterator_combination,old_variable_rotamers,old_value)
		print 'new attempt, ', new_attempt


		current_variable_combination=new_attempt[0]
		current_interactions_value=new_attempt[1]

		integrated_value=current_interactions_value
		print 'Integrated value:  ', integrated_value

		provided_set=(set(current_iterator_combination) | set(current_variable_combination)	)
		print 'Provided set: ', provided_set
		integrated_set=provided_set

		if leaf_indicator==0:
		    print 'Not a leaf, going one level down.....'
		    for current_child in children_of_current_node:
			for dummytuple in children_dict[current_child]:
			    dummyset=set(dummytuple)
			    if dummyset.issubset(provided_set):
				integrated_set= ( set(children_dict[current_child][dummytuple][0]) | integrated_set)
				integrated_value=integrated_value+children_dict[current_child][dummytuple][1]
				print 'Integrated set, and integrated value: ', integrated_set, '---', integrated_value

		if integrated_value < smallest_value:
		    print 'Will CHANGE THE SMALLEST VALUE'
		    smallest_value=integrated_value
		    optimal_integrated_combination=integrated_set

	    output_dictionary[tuple(current_iterator_combination)]=[tuple(optimal_integrated_combination), smallest_value]

    else:
	print '-------------------------------------------------------------------------------REACHED THE ROOT-------------------------'
	TimeWaster=raw_input('ROOT REACHED press ENTER.....')
	current_iterator_combination=[]
	optimal_variable_combination=list()
	smallest_value=sys.maxint
	current_variable_combination=pick_random_initial_variable_combination(variable_dictionary)
	current_interactions_value=find_total_combination_value(interaction_dictionary,current_iterator_combination,current_variable_combination)

	for trial_number in range(number_of_trials):

	    print 'Trial number ----> ', trial_number
	    old_variable_rotamers=current_variable_combination
	    old_value=current_interactions_value
	    print 'Old variable rotamers and old value--->', old_variable_rotamers, '^^^^^', old_value

	    new_attempt=change_variable_rotamer(interaction_dictionary,variable_dictionary,current_iterator_combination,old_variable_rotamers,old_value)
	    current_variable_combination=new_attempt[0]
	    current_interactions_value=new_attempt[1]
	    integrated_value=current_interactions_value
	    print 'new attempt, combination and value: ', current_variable_combination, '****', current_interactions_value

	    provided_set=set(current_variable_combination)
	    integrated_set=provided_set

	    if leaf_indicator==0:
		print 'Not a leaf, going down'
		for current_child in children_of_current_node:
		    for dummytuple in children_dict[current_child]:
			dummyset=set(dummytuple)
			if dummyset.issubset(provided_set):
			    integrated_set= ( set(children_dict[current_child][dummytuple][0]) | integrated_set)
			    integrated_value+=children_dict[current_child][dummytuple][1]
			    print 'After that child, integrated set and value:', integrated_set,'++', integrated_value

	    if integrated_value < smallest_value:
		print 'FOUND SMALLER VALUE, UPDATING ROOT'
		smallest_value=integrated_value
		optimal_integrated_combination=integrated_set

	output_dictionary[tuple(optimal_integrated_combination)]=smallest_value

    return output_dictionary



######################################################################################################################################################
######################################################################################################################################################



######################################################################################################################################################
######################################################################################################################################################
#########################TEST DICT----REMOVE AFTER DEBUGGING####################################################################################
######################################################################################################################################################
######################################################################################################################################################



input_dictionary={'a':[18,6], 'b':[19,4], 'c':[17,7],'d':[5,20],'e':[26,8],'f':[16,3],'g':[9,21],'h':[15,10], 'i':[14,2], 'j':[11,23],'k':[22,12], 'l':[13,24],'m':[25,1],  }

dict2=dict()

dict2[tuple([4,6])]=2
dict2[tuple([4,18])]=14
dict2[tuple([19,6])]=13
dict2[tuple([19,18])]=1

dict2[tuple([4,5])]=1
dict2[tuple([4,20])]=16
dict2[tuple([5,19])]=14
dict2[tuple([19,20])]=1

dict2[tuple([6,7])]=1
dict2[tuple([6,17])]=11
dict2[tuple([7,18])]=11
dict2[tuple([17,18])]=1

dict2[tuple([7,13])]=6
dict2[tuple([7,24])]=17
dict2[tuple([17,13])]=4
dict2[tuple([17,24])]=7

dict2[tuple([7,12])]=5
dict2[tuple([7,22])]=15
dict2[tuple([17,12])]=5
dict2[tuple([17,22])]=5

dict2[tuple([12,13])]=1
dict2[tuple([12,24])]=12
dict2[tuple([13,22])]=9
dict2[tuple([24,22])]=2

dict2[tuple([8,11])]=3
dict2[tuple([8,23])]=15
dict2[tuple([26,11])]=15
dict2[tuple([26,23])]=3

dict2[tuple([8,2])]=6
dict2[tuple([8,14])]=6
dict2[tuple([26,2])]=24
dict2[tuple([26,14])]=12

dict2[tuple([2,11])]=9
dict2[tuple([2,23])]=21
dict2[tuple([14,11])]=3
dict2[tuple([14,23])]=9

dict2[tuple([10,3])]=7
dict2[tuple([10,16])]=6
dict2[tuple([15,3])]=12
dict2[tuple([15,16])]=1

dict2[tuple([3,9])]=6
dict2[tuple([3,21])]=18
dict2[tuple([16,9])]=7
dict2[tuple([16,21])]=5

dict2[tuple([5,3])]=2
dict2[tuple([5,16])]=11
dict2[tuple([20,3])]=17
dict2[tuple([20,16])]=4

dict2[tuple([7,8])]=1
dict2[tuple([7,26])]=19
dict2[tuple([17,8])]=9
dict2[tuple([17,26])]=9

dict2[tuple([3,8])]=5
dict2[tuple([3,26])]=23
dict2[tuple([16,8])]=8
dict2[tuple([16,26])]=10

dict2[tuple([5,7])]=2
dict2[tuple([5,17])]=12
dict2[tuple([20,7])]=13
dict2[tuple([20,17])]=3

dict2[tuple([7,1])]=6
dict2[tuple([7,25])]=18
dict2[tuple([17,1])]=16
dict2[tuple([17,25])]=8

dict2[tuple([3,1])]=2
dict2[tuple([3,25])]=22
dict2[tuple([16,1])]=15
dict2[tuple([16,25])]=9

dict2[tuple([5,1])]=4
dict2[tuple([5,25])]=20
dict2[tuple([20,1])]=19
dict2[tuple([20,25])]=5

dict2[tuple([8,1])]=7
dict2[tuple([8,25])]=17
dict2[tuple([26,1])]=25
dict2[tuple([26,25])]=1

dict2[tuple([8,11])]=3
dict2[tuple([8,23])]=15
dict2[tuple([26,11])]=15
dict2[tuple([26,23])]=3

dict2[tuple([2,11])]=9
dict2[tuple([2,23])]=21
dict2[tuple([14,11])]=3
dict2[tuple([14,23])]=9

interaction_dictionary=dict2
######################################################################################################################################################
######################################################################################################################################################



G=nx.Graph()
G=interaction_pair_graph_builder(input_dictionary, interaction_dictionary)
#nx.draw_spring(G)
#plt.show()
decomp_tree=tree_decomposition(G)
#nx.draw_spring(decomp_tree)
#plt.show()
#TimeWaster=raw_input('After you are done admiring the plots, press ENTER.....')
Dynamic_Programming_for_decomposed_trees(decomp_tree,input_dictionary,interaction_dictionary)

import pickle
resi2rots = pickle.load(open('resi2rots.pkl','rb'))
Edict = pickle.load(open('Edict.pkl','rb'))
#G = pickle.load(open('G.pkl','rb'))
H = nx.Graph()
H = interaction_pair_graph_builder(resi2rots, Edict)
decomp_tree = tree_decomposition(H)
Dynamic_Programming_for_decomposed_trees(decomp_tree, resi2rots, Edict)





