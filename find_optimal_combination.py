import networkx as nx
import random
from operator import itemgetter
from numpy import *
import sys


######################################################################################################################################################
######################################################################################################################################################

def counter(input_dict_of_lists):

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
	    print 'repetition TEMP', temp_repetition
	    print 'repetition dict', repetition_dict
	else:
	    repetition_dict[m]=temp_repetition
	    temp_repetition=temp_repetition*cardinality_dict[m]
	    print 'repetition TEMP', temp_repetition
	    print 'repetition dict', repetition_dict



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



def find_optimal_combination(input_dictionary,interaction_dictionary,list_of_iterators, list_of_variables):

    iterator_dictionary=dict()
    variable_dictionary=dict()
    output_dictionary=dict()

    for iterator in list_of_iterators:
	iterator_dictionary[iterator]=input_dictionary[iterator]

    print 'ITERATOR DICT  :   ', iterator_dictionary

    for variable in list_of_variables:
	variable_dictionary[variable]=input_dictionary[variable]

    print 'VARIABLE DICT  :   ', variable_dictionary

    all_iterator_combinations=counter(iterator_dictionary)
    all_variable_combinations=counter(variable_dictionary)

    print 'ITERATOR COMBINATIONS  :   ', all_iterator_combinations
    print 'all_variable_combinations  :	  ', all_variable_combinations

    for current_iterator_combination in all_iterator_combinations:

	optimal_variable_combination=list()
	smallest_value=sys.maxint     #A very large value

	for current_variable_combination in all_variable_combinations:

	    temp_value=find_total_combination_value(interaction_dictionary,current_iterator_combination,current_variable_combination)
	    print '  current_variable_combination     temp_value =  ', current_variable_combination, temp_value

	    if temp_value < smallest_value:
		smallest_value=temp_value
		print 'smallest_value=	', smallest_value
		optimal_variable_combination=current_variable_combination
		print 'optimal_variable_combination :	', optimal_variable_combination

	output_dictionary[tuple(current_iterator_combination)]=optimal_variable_combination
	print output_dictionary


    return output_dictionary




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
		#print 'k,m ';print k;print  m ;
		#print 'VALUE ADDED:	 ', interaction_dictionary[tuple([total_list[k],total_list[m]])];
		print 'OUTPUT:	 ', output

    return output





######################################################################################################################################################


aaa=dict()
aaa['a']=[1,2,3]
aaa['b']=[4,5]
aaa['c']=[7,9]
aaa['d']=[11,13]

list1=['a','b']
list2=['c','d']

ddd=dict()
for k in range(15):
    for m in range(15):
	ddd[tuple([k,m])]=exp(2*cos(k+m))


#find_total_combination_value(ddd,[1,3],[2,5])


find_optimal_combination(aaa,ddd,list1, list2)





