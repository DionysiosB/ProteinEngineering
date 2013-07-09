import random
import itertools

#########################################################################################################################################################################################
#########################################################################################################################################################################################
def simple_array_transpose(input_array):
	
	number_of_columns=len(input_array)
	number_of_rows=len(input_array[0])
	output_array=list()
	
	for row_number in range(number_of_rows):
		
		current_row=list()
		for column_number in range(number_of_columns):
			current_row.append(input_array[column_number][row_number])
		
		output_array.append(current_row)
	
	#print output_array
	return output_array
#########################################################################################################################################################################################
#########################################################################################################################################################################################
def list_modulo_operation(input_list, modulo_number):

	output_list=list()
	for temp in input_list:
		output_list.append(temp % modulo_number)
	
	return output_list


#########################################################################################################################################################################################
#########################################################################################################################################################################################
def make_unique_combination_array_odd_numbers(input_integer,first_matrix_column):
	
	output_list_of_lists=list()
	output_list_of_lists.append(range(input_integer))
	
	for first_row_element in first_matrix_column:
		output_list_of_lists.append(list_modulo_operation(range(first_row_element,first_row_element+input_integer),input_integer))
	
	return output_list_of_lists
#########################################################################################################################################################################################
#########################################################################################################################################################################################
def make_protein_sequences(input_integer,intermediate_array):
	
	output_array=list()
	#intermediate_array=make_unique_combination_array_odd_numbers(input_integer)
	intermediate_array.pop(0)
	#print intermediate_array
	integer_sequence=range(input_integer)
	
	first_column=list()
	for k in integer_sequence:
		first_column.extend(input_integer*[k])
		
	second_column=integer_sequence*input_integer
	
	output_array.append(first_column)
	output_array.append(second_column)
	
	
	first_index=0
	
	for dummy in range(input_integer-1):
		
		current_column=list()
		current_column.extend(integer_sequence)
		for k in range(input_integer-1):
			current_column.extend(intermediate_array[(k+first_index)%(input_integer-1)])
		
		output_array.append(current_column)
		first_index+=1
	
	
	output_array_transpose=simple_array_transpose(output_array)
	#print '------------------------------------------------'
	#print output_array_transpose
	
	return output_array_transpose
	


######################################################################################################################################################
######################################################################################################################################################
def testing_function(list_of_variables,string_length, list_of_lists):
	
	
	for test_list in list_of_lists:
		if len(test_list) != string_length:
			print 'At least one of the given strings has not the right length. Please try again with another set of strings. Exiting....'
			return -1
		
		
	all_possible_pair_combinations=list()
	for first_element in list_of_variables:
		for second_element in list_of_variables:
			all_possible_pair_combinations.append([first_element,second_element])
	
	
	for first_position in range(string_length-1):
		for second_position in range(first_position+1, string_length):
			existing_pairs=list()
			for temp_list in list_of_lists:
				additional_pair=[temp_list[first_position], temp_list[second_position]]
				existing_pairs.append(additional_pair)
			
			for pair_combination in all_possible_pair_combinations:
				if pair_combination not in existing_pairs:
					#print 'FIRST POSITION:', first_position, '   ', 'SECOND POSITION: ', second_position
					#print 'Combination ', pair_combination , 'does not exist in the existing lists. Exiting...'
					return -1
	
	
	print 'All possible combinations have been found.'
	return 0



######################################################################################################################################################
######################################################################################################################################################
###############################################################################################################################################
######################################################################################################################################################
def valid_permutations(n):
	
	input_vector=range(n)
	dummy=itertools.permutations(input_vector,n)
	for current_permutation in dummy:
		flag=1
		for k in range(n):
			if current_permutation[k] == k:flag=0
		if flag==1:
			print current_permutation
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

def modulo_addup_to_list(input_list,number,modulo):
	
	output_list=list()
	for temp in input_list:
		new_number=(temp+number)%modulo
		output_list.append(new_number)
	
	return output_list

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
def list_cyclic_shifts(input_list):
	
	output_list=list()
	list_length=len(input_list)
	for k in range(1,list_length):
		new_list=modulo_addup_to_list(input_list,k,list_length)
		output_list.append(new_list)
	
	return output_list
#######################################################################################################################
#######################################################################################################################
def compatibility_test(list_group, n):
	
	allowed_positions=list()
	for q in range(n):
		allowed_positions.append(range(n))
		allowed_positions[q].remove(q)
	#print allowed_positions
	
	for inspected_list in list_group:
		for k in range(n):
			#print 'k ', k
			if inspected_list[k] in allowed_positions[k]:
				allowed_positions[k].remove(inspected_list[k])
				#print 'allowed_positions', allowed_positions
			else:
				#print 'bye'
				return -1
	
	#print 'Passed first test: ', list_group
	
	dummy_products=itertools.product(range(n),repeat=2)
	allowed_pairs=list()
	for temp in dummy_products:
		if temp[0] != temp[1]:
			allowed_pairs.append(list(temp))
	#print allowed_pairs
	
	for list_index in range(n-1):
		for k in range(1,n):
			#print 'aaa, ',k 
			new_pair=[list_group[list_index][k],list_group[(list_index+1)%(n-1)][k]]
			if new_pair in allowed_pairs:
				allowed_pairs.remove(new_pair)
				#print 'Removed Pair: ', new_pair
			else:
				#print list_group[list_index][k]
				#print list_group[list_index+1][k]
				#print 'bye'
				return -2
				
	
	#print 'Passed Second test: ', list_group
	return 0
	

######################################################################################################################################################
######################################################################################################################################################
def list_differences(input_list,modulo_integer):
	
	list_length=len(input_list)
	temp_list=list()
	output_list=list()
	
	for k in range(list_length-1):
		temp_list.append(input_list[k+1]-input_list[k])
	temp_list.append(input_list[0]-input_list[-1])
	for k in temp_list:
		output_list.append(k%modulo_integer)
	
	return output_list

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
def first_conjecture_search(input_integer):
	
	if input_integer%2==0:
		print 'The conjecture only tries odd numbers. Exiting...'
		return -1
	
	all_cyclic_shifts=list_cyclic_shifts(range(input_integer))
	
	number_of_rows=len(all_cyclic_shifts)
	all_permutations=itertools.permutations(range(number_of_rows),number_of_rows)
	
	for current_permutation in all_permutations:
		
		#print current_permutation
		current_list=list()
		for k in current_permutation:
			current_list.append(all_cyclic_shifts[k])
		#print current_list
		#print '----------------------------------------------'
		
		if current_list[0][0]>1:
			continue
		
		test_result=compatibility_test(current_list,input_integer)
		if test_result==0:
			#print '---------------------------------------------'
			#print current_list
			firsts_list=list()
			for k in current_list:
				firsts_list.append(k[0])
			#print 'Firsts List --------->' , firsts_list
			
			aaa=make_unique_combination_array_odd_numbers(input_integer,firsts_list)
			bbb=make_protein_sequences(input_integer,aaa)
			test_result_new=testing_function(range(input_integer),input_integer+1, bbb)
			if test_result_new==0:
				print 'THIS PERMUTATION WORKS------> ',firsts_list
				print '-----------------------------------------'
	
	

######################################################################################################################################################
######################################################################################################################################################

first_conjecture_search(11)



