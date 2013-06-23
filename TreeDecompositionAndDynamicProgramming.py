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
	#print current_graph.order()
	nodes_sorted_by_degree=sort_by_degree(current_graph)
	#print 'nodes_sorted_by_degree', nodes_sorted_by_degree
	minimum_degree_vertex=nodes_sorted_by_degree[0][0]
	#print 'Minimum Degree_vertex' , minimum_degree_vertex
	cliques_of_minimum_degree_vertex=nx.cliques_containing_node(current_graph,minimum_degree_vertex)
	#print 'cliques_of_minimum_degree_vertex',cliques_of_minimum_degree_vertex
	number_of_cliques_containing_vertex=len(cliques_of_minimum_degree_vertex)
	#print 'number_of_cliques_containing_vertex', number_of_cliques_containing_vertex
	minimum_degree_vertex_neighbors=current_graph.neighbors(minimum_degree_vertex)
	#print 'minimum_degree_vertex_neighbors', minimum_degree_vertex_neighbors
	new_tree_vertex=[minimum_degree_vertex]
	#print 'new_tree_vertex First element: ',new_tree_vertex
	new_tree_vertex.extend(minimum_degree_vertex_neighbors)
	new_tree_vertex=tuple(new_tree_vertex)
	decomposition_tree.add_node(new_tree_vertex)
	#print 'decomposition_tree_vertices',decomposition_tree.nodes()
	if number_of_cliques_containing_vertex>1:
	    #print 'Not Clique, will remove only one vertex'
	    pairs_of_neighbors=make_pairs(minimum_degree_vertex_neighbors)
	    #print 'pairs_of_neighbors',pairs_of_neighbors
	    for additional_edge in pairs_of_neighbors:current_graph.add_edge(additional_edge[0],additional_edge[1])
	    toberemoved=[minimum_degree_vertex]
	    #print 'toberemoved ', toberemoved
	else:
	    toberemoved=[minimum_degree_vertex]
	    #print 'Clique detected, will try to remove more than one vertex'
	    number_of_clique_edges_per_vertex=len(minimum_degree_vertex_neighbors)
	    #print 'number_of_clique_edges_per_vertex',number_of_clique_edges_per_vertex
	    #print 'Checking all the vertex`s neighbors...'
	    #print 'minimum_degree_vertex_neighbors', minimum_degree_vertex_neighbors
	    for temp_vertex in minimum_degree_vertex_neighbors:
		if current_graph.degree(temp_vertex)==number_of_clique_edges_per_vertex:
		    toberemoved.append(temp_vertex)
		    print 'Will ALSO remove vertex  ', temp_vertex
	for graph_vertex in new_tree_vertex:
	    if graph_vertex in toberemoved:
		current_graph.remove_node(graph_vertex)
		#print 'Removed original graph vertex', graph_vertex
		tree_vertices_waiting=tree_connectivity_dictionary[graph_vertex]
		#print 'For the removed node, tree_vertices_waiting: ' , tree_vertices_waiting
		for tree_vertex_waiting in tree_vertices_waiting:
		    #print 'New Tree vertex:	 ' , new_tree_vertex
		    #print 'Tree Vertex waiting:', tree_vertex_waiting
		    decomposition_tree.add_edge(new_tree_vertex,tree_vertex_waiting)
		    #print 'Connected tree vertices', new_tree_vertex, 'and  ' ,   tree_vertex_waiting
		    #print 'The tree edges are now:   ', decomposition_tree.edges()
		    #print 'THE NUMBER OF TREE EDGES ARE NOW:	', len(decomposition_tree.edges())
		for tree_vertex_waiting in tree_vertices_waiting:
		    common_graph_nodes_between_tree_vertices=list(my_very_simple_tuple_intersection(new_tree_vertex,tree_vertex_waiting))
		    for graph_vertex in common_graph_nodes_between_tree_vertices:
			tree_connectivity_dictionary[graph_vertex].remove(tree_vertex_waiting)
			#print 'Removed from dictionary entry', graph_vertex , 'tree node ', tree_vertex_waiting
			#print 'Now the new dictionary is:  ' , tree_connectivity_dictionary


	    else:
		tree_connectivity_dictionary[graph_vertex].append(new_tree_vertex)
		#print 'New tree_connectivity_dictionary node appended. New tree_connectivity_dictionary ', tree_connectivity_dictionary
	#print 'tree_connectivity_dictionary:  '     , tree_connectivity_dictionary
	#print 'decomposition_tree.nodes:    ', decomposition_tree.nodes()
	#print 'decomposition_tree.edges:    ', decomposition_tree.edges()



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

    print 'CARDINALITY DICT' , cardinality_dict
    repetition_dict=dict()
    temp_repetition=0
    for m in cardinality_dict:
	if temp_repetition==0:
	    repetition_dict[m]=1
	    temp_repetition=cardinality_dict[m]
	    #print 'repetition TEMP', temp_repetition
	    #print 'repetition dict', repetition_dict
	else:
	    repetition_dict[m]=temp_repetition
	    temp_repetition=temp_repetition*cardinality_dict[m]
	    #print 'repetition TEMP', temp_repetition
	    #print 'repetition dict', repetition_dict



    total_number_of_combinations=temp_repetition
    print 'total_number_of_combinations=  ', total_number_of_combinations

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
    print 'MASTER DICTIONARY = ', master_dictionary



    tree_root=find_optimal_tree_root(input_tree)
    print 'Tree root is	 ', tree_root
    next_tree_leaves=find_tree_leaves(current_tree)
    current_tree_leaves=find_tree_leaves(current_tree)


    [tree_structure_children_to_parent,tree_structure_parent_to_children]=find_tree_structure(input_tree)

    print 'tree_structure_children_to_parent   ',tree_structure_children_to_parent
    print 'tree_structure_parent_to_children:  ',tree_structure_parent_to_children



    print' ############################################################################################################################################################'

    while len(current_tree_leaves)>0:

	current_tree_leaves=next_tree_leaves[:]
	next_tree_leaves=list()
	if tree_root in current_tree_leaves: current_tree_leaves.remove(tree_root)     #The root HAS to be computed after ALL the other nodes are computed
	print 'REMOVED TREE ROOT'
	print 'Current_tree_leaves ', current_tree_leaves

	for current_node in current_tree_leaves:
	    print 'Current Node:   ', current_node
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
    print 'FINAL DICTIONARY  ', final_dictionary

    best_combination=final_dictionary.keys()[0]
    minimum_value=final_dictionary[best_combination]



    return [best_combination, minimum_value]

######################################################################################################################################################
######################################################################################################################################################
def find_optimal_combination(input_dictionary,interaction_dictionary,current_node,parent_of_node,children_dict):

    print 'ENTERED find_optimal_combination FUNCTION'
    print 'Input dictionary: ', input_dictionary
    #print 'Interaction dictionary:  ', interaction_dictionary
    print 'Current_node:  ' , current_node
    print 'Parent of node:  ', parent_of_node
    print 'Children_dict: ', children_dict

    if parent_of_node != -1:
	node_with_parent_intersection=tuple( set(current_node) & set(parent_of_node)  )
	print 'node_with_parent_intersection  ', node_with_parent_intersection
	node_not_parent_elements=tuple( set( current_node) - set(parent_of_node))
	print 'node_not_parent_elements	 ', node_not_parent_elements
    else:
	node_with_parent_intersection=tuple()
	print 'node_with_parent_intersection  ', node_with_parent_intersection
	node_not_parent_elements=current_node
	print 'node_not_parent_elements	 ', node_not_parent_elements


    if len(children_dict)>0:
	leaf_indicator=0
	print 'leaf_indicator  ', leaf_indicator
	children_of_current_node=children_dict.keys()
	print 'children_of_current_node:  ', children_of_current_node
    else:
	leaf_indicator=1
	print 'leaf_indicator  ', leaf_indicator


    iterator_dictionary=dict()
    variable_dictionary=dict()
    output_dictionary=dict()


    if parent_of_node != -1:
	for iterator in node_with_parent_intersection:
	    iterator_dictionary[iterator]=input_dictionary[iterator]


    print 'ITERATOR DICT  :	  ', iterator_dictionary

    for variable in node_not_parent_elements:
	variable_dictionary[variable]=input_dictionary[variable]

    print 'VARIABLE DICT  :	  ', variable_dictionary

    all_iterator_combinations=find_combinations_list(iterator_dictionary)
    all_variable_combinations=find_combinations_list(variable_dictionary)

    print 'ITERATOR COMBINATIONS  :	  ', all_iterator_combinations
    print 'all_variable_combinations  :	  ', all_variable_combinations


    if len(all_iterator_combinations)>0:
	for current_iterator_combination in all_iterator_combinations:

	    optimal_variable_combination=list()
	    smallest_value=sys.maxint

	    for current_variable_combination in all_variable_combinations:

		current_node_interactions_value=find_total_combination_value(interaction_dictionary,current_iterator_combination,current_variable_combination)
		print 'node_interactions_value =' ,current_node_interactions_value

		integrated_value=current_node_interactions_value

		if leaf_indicator==0:
		    print 'THIS IS NOT A LEAF, SO IT HAS CHILDREN.....'
		    for current_child in children_of_current_node:

			print 'current_child', current_child

			provided_set=(set(current_iterator_combination) | set(current_variable_combination)  )
			print 'provided_set', provided_set
			integrated_set=provided_set
			print 'Integrated set:	 ', integrated_set

			for dummytuple in children_dict[current_child]:
			    dummyset=set(dummytuple)
			    if dummyset.issubset(provided_set):
				print 'THIS dummyset IS SUBSET	 : ',  dummyset
				integrated_set= ( set(children_dict[current_child][dummytuple][0]) | integrated_set)
				print 'integrated_set  : ', integrated_set
				integrated_value+=children_dict[current_child][dummytuple][1]
				print 'integrated_value', integrated_value
		else:
		    print 'THIS IS A LEAF, NO CHILDREN, NO RECURSIVE FUNCTIONS'
		    provided_set=(set(current_iterator_combination) | set(current_variable_combination)	     )
		    integrated_set=provided_set
		if integrated_value < smallest_value:
		    print 'SMALLEST VALUE HAS TO BE UPDATED'
		    smallest_value=integrated_value
		    print 'smallest_value=  ', smallest_value
		    optimal_integrated_combination=integrated_set
		    print 'optimal_integrated_combination : ', optimal_integrated_combination

	    output_dictionary[tuple(current_iterator_combination)]=[tuple(optimal_integrated_combination), smallest_value]
	    print 'output_dictionary', output_dictionary




    else:
	print 'BBBBBBBBBBBBBBBBBBB'

	current_iterator_combination=[]

	optimal_variable_combination=list()
	smallest_value=sys.maxint

	for current_variable_combination in all_variable_combinations:
	    current_node_interactions_value=find_total_combination_value(interaction_dictionary,current_iterator_combination,current_variable_combination)

	    print 'node_interactions_value =' ,current_node_interactions_value
	    integrated_value=current_node_interactions_value
	    provided_set=set(current_variable_combination)
	    print 'provided_set', provided_set
	    integrated_set=provided_set
	    print 'Integrated set:		 ', integrated_set

	    if leaf_indicator==0:
		print 'THIS IS NOT A LEAF, SO IT HAS CHILDREN.....'
		for current_child in children_of_current_node:

		    print 'current_child', current_child

		    for dummytuple in children_dict[current_child]:
			dummyset=set(dummytuple)
			#print 'dummyset:  ', dummyset
			if dummyset.issubset(provided_set):
			    print 'THIS dummyset IS SUBSET   : ',  dummyset
			    print 'THE NEW INTEGRATED SET WILL BE THE UNION OF THE FOLLOWING SETS:', set(children_dict[current_child][dummytuple][0]), integrated_set
			    integrated_set= ( set(children_dict[current_child][dummytuple][0]) | integrated_set)
			    print 'integrated_set  : ', integrated_set
			    integrated_value+=children_dict[current_child][dummytuple][1]
			    print 'integrated_value', integrated_value
	    else:
		print 'THIS IS A LEAF, NO CHILDREN, NO RECURSIVE FUNCTIONS'
		provided_set=(set(current_iterator_combination) | set(current_variable_combination)	     )
		integrated_set=provided_set

	    if integrated_value < smallest_value:
		print 'SMALLEST VALUE HAS TO BE UPDATED'
		smallest_value=integrated_value
		print 'smallest_value=	    ', smallest_value
		optimal_integrated_combination=integrated_set
		print 'optimal_integrated_combination : ', optimal_integrated_combination

	output_dictionary[tuple(optimal_integrated_combination)]=smallest_value
	print 'output_dictionary', output_dictionary


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
		#print 'k,m	';print k;print	 m ;
		#print 'VALUE ADDED:	 ', interaction_dictionary[tuple([total_list[k],total_list[m]])];
		#print 'OUTPUT:	 ', output

    return output

######################################################################################################################################################
######################################################################################################################################################


test_tree=nx.Graph()
node0=tuple(['a','b'])
node1=tuple(['a','c'])
node2=tuple(['b','f'])
test_tree.add_node(node0)
test_tree.add_node(node1)
test_tree.add_node(node2)
test_tree.add_edge(node0,node1)
test_tree.add_edge(node0,node2)
dict1={'g':[37,41], 'a':[3,5], 'b':[7,11], 'c':[13,17],'d':[19,23], 'f':[29,31] }
dict2=dict()
for k in range(45):
    for m in range(45):
	dict2[tuple([k,m])]=abs(k-m)

print test_tree.nodes()
print test_tree.edges()


[X, Y]=Dynamic_Programming_for_decomposed_trees(test_tree,dict1,dict2)
print X
print Y