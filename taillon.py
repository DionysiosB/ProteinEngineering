import networkx as nx
import random
from operator import itemgetter

#############################################################################################

def sort_by_degree(G):
    return sorted(G.degree(with_labels=True).items(),key = itemgetter(1))

##############################################################################################

def dimacs2nx(filename):
    G = nx.Graph()
    for line in open(filename).readlines():
	l = line.split()
	if l[0]=='c': continue
	if l[0]=='p':
	    N = int(l[2])
	    for n in range(N):
		G.add_node(n)
	if l[0]=='e':
	    G.add_edge(int(l[1]),int(l[2]))
    return G


##############################################################################################

def make_nx_graph(term2energy,hirot2resid):
    termlist = term2energy.keys()

    # Each correction term is a node in a graph
    G = nx.Graph()
    for term in termlist:
	G.add_node(term)

    # We start by drawing edges between nodes that collide (i.e. are not disjoint)
    # For the possibility of overcounting the triplets must have two hirotamers in common, and the union of residues must number 4
    for i,term1 in enumerate(termlist):
	res1 = [hirot2resid[x] for x in term1]
	for j,term2 in enumerate(termlist):
	    if j <= i: continue
	    res2 = [hirot2resid[x] for x in term2]
	    if len(set(term1).intersection(term2)) == 2 and len(set(res1).union(res2)) == 4:
		G.add_edge(term1,term2)

    print 'Made graph with ',len(G.nodes()),'nodes','and',len(G.edges()),'edges'
    return G



##############################################################################################

def make_nx_graph(term2energy,hirot2resid,verbose=False):
    termlist = term2energy.keys()

    # Each correction term is a node in a graph
    G = nx.Graph()
    for term in termlist:
	G.add_node(term)

    # We start by drawing edges between nodes that collide (i.e. are not disjoint)
    # For the possibility of overcounting the triplets must have two hirotamers in common, and the union of residues must number 4
    for i,term1 in enumerate(termlist):
	res1 = [hirot2resid[x] for x in term1]
	for j,term2 in enumerate(termlist):
	    if j <= i: continue
	    res2 = [hirot2resid[x] for x in term2]
	    #if len(set(term1).intersection(term2)) == 2 and len(set(res1).union(res2)) == 4:
	    sharedhirots = set(term1).intersection(term2)
	    allhirots = set(term1).union(term2)
	    allresids = set(res1).union(res2)
	    if len(sharedhirots) > 1 and len(allresids) == len(allhirots):
		G.add_edge(term1,term2)

    if verbose: print 'Made graph with ',len(G.nodes()),'nodes','and',len(G.edges()),'edges'
    return G



##############################################################################################


def heuristic_independent_set_finder(Graph,Weight_vector,first_node,second_node_index):

    first_node_neighborhood=Graph[first_node].keys()

    second_node=first_node_neighborhood[second_node_index]
    second_node_neighborhood=Graph[second_node].keys()

    #print 'second node',second_node
    #print 'first node_neighborhood', first_node_neighborhood
    #print 'second node_neighborhood', second_node_neighborhood

    first_node_neighborhood_set=set(first_node_neighborhood)
    second_node_neighborhood_set=set(second_node_neighborhood)

    common_neighbors=first_node_neighborhood_set.intersection(second_node_neighborhood_set)

    #print 'common neighbors:', common_neighbors

    induced_subgraph=nx.subgraph(Graph,common_neighbors)
    complement_induced_subgraph=nx.complement(induced_subgraph)

    #print 'induced_subgraph.nodes()', induced_subgraph.nodes()
    #print 'induced_subgraph.edges()', induced_subgraph.edges()
    #print 'complement_induced_subgraph.nodes()', complement_induced_subgraph.nodes()
    #print 'complement_induced_subgraph.edges()', complement_induced_subgraph.edges()

    current_complement_induced_subgraph=nx.complement(induced_subgraph)
    current_independent_set=[]

    print current_complement_induced_subgraph.nodes()
    print current_complement_induced_subgraph.edges()

    while (len(current_complement_induced_subgraph.edges())>0):
	current_minimum_degree_node=sort_by_degree(current_complement_induced_subgraph)[0][0]
	current_minimum_degree_node_neighborhood=current_complement_induced_subgraph[current_minimum_degree_node].keys()
	current_independent_set.append(current_minimum_degree_node)
	current_complement_induced_subgraph.remove_nodes_from(current_minimum_degree_node_neighborhood)
	current_complement_induced_subgraph.remove_node(current_minimum_degree_node)

	print current_minimum_degree_node
	print current_minimum_degree_node_neighborhood
	print current_independent_set
	print current_complement_induced_subgraph.nodes()

    current_independent_set.extend(current_complement_induced_subgraph.nodes())
    current_independent_set.extend([first_node,second_node])
    maximum_clique=current_independent_set
    return maximum_clique



################################################################################################

def taillon_algorithm(Graph,weight_vector,number_of_iterations=-1):

    maximum_clique=[]

    weights_dictionary=dict( (node, weight_vector[node]) for node in Graph.nodes() )
    #neighborweights = dict( (term,F[term]) for term in N_v)
    ordered_weights_dictionary=sorted(weights_dictionary.items(),key=itemgetter(1),reverse=True)

    if number_of_iterations==-1:number_of_iterations=Graph.order()

    completed_iterations=0

    while (completed_iterations<number_of_iterations):

	first_node=random.choice(Graph.nodes())
	print 'First node:', first_node
	first_node_neighborhood=Graph[first_node]
	number_of_first_node_neighbors=len(first_node_neighborhood)

	if (number_of_first_node_neighbors>1):

	    second_node_index=random.choice(range(number_of_first_node_neighbors))
	    print 'SECOND node index:', second_node_index

	    current_maximum_clique=heuristic_independent_set_finder(Graph,weight_vector,first_node,second_node_index)
	    if len(current_maximum_clique)>len(maximum_clique): maximum_clique=current_maximum_clique
	    completed_iterations=completed_iterations+1

    return maximum_clique





######################################################################################################







