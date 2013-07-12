import networkx as nx
import random
from operator import itemgetter
from numpy import *

######################################################################################################################################################
def counter(input_list_of_lists):

    number_of_sets=len(input_list_of_lists)
    cardinality_vector=list()
    for k in range(number_of_sets):
	cardinality_vector.append(len(input_list_of_lists[k]))
	if cardinality_vector[k]<=0:
	    print 'The elements of the list must be strictly positive integers. Exiting....'
	    return -1

    repetition_vector=[1]
    for m in range(1,number_of_sets):
	repetition_vector.append(repetition_vector[m-1]*cardinality_vector[m-1])

    total_number_of_combinations=repetition_vector[number_of_sets-1]*cardinality_vector[number_of_sets-1]
    print 'total_number_of_combinations=  ', total_number_of_combinations

    output_combination_list=list()

    for combination_number in range(total_number_of_combinations):
	current_combination_list=list()
	for set_number in range(number_of_sets):
	    current_combination_list.append( input_list_of_lists[set_number][(combination_number/repetition_vector[set_number])%cardinality_vector[set_number]])

	output_combination_list.append(current_combination_list)

    return output_combination_list
######################################################################################################################################################
