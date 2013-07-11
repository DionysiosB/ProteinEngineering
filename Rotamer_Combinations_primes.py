######################################################################################################################################################
######################################	OPTIMAL ROTAMER COMBINATION SEQUENCES  #######################################################################
######################################		DIONYSIOS BARMPOUTIS	       #######################################################################
######################################################################################################################################################


######################################################################################################################################################
######################################################################################################################################################
def testing_function(list_of_variables,string_length, list_of_lists):

    if list_of_lists==-1:
	print 'No sequences passed to the testing function. Exiting...'
	return -1

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
		    #print 'FIRST POSITION:', first_position, '	  ', 'SECOND POSITION: ', second_position
		    #print 'Combination ', pair_combination , 'does not exist in the existing lists. Exiting...'
		    return -1


    print 'All possible combinations have been found.'
    return 1

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
def make_unique_combination_array_odd_numbers(input_integer):

    print "Calculating the seed permutation. This may take some time...."
    first_matrix_column=output_sequence_seed_permutation(input_integer)
    if first_matrix_column==-1:return -1
    print 'Seed permutation found: ', first_matrix_column

    output_list_of_lists=list()
    output_list_of_lists.append(range(input_integer))

    for first_row_element in first_matrix_column:
	output_list_of_lists.append(list_modulo_operation(range(first_row_element,first_row_element+input_integer),input_integer))

    return output_list_of_lists
#########################################################################################################################################################################################
#########################################################################################################################################################################################
def make_protein_sequences(input_integer):

    output_array=list()
    intermediate_array=make_unique_combination_array_odd_numbers(input_integer)
    if intermediate_array==-1:return -1
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
    print '------------------------------------------------'
    print output_array_transpose

    return output_array_transpose

#####################################################################################################################################################
#####################################################################################################################################################
def primetest(number_trial):

    #HARDWIRED PRIME NUMBERS, NEED TO BE CHANGED
    primes_list=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,137,139,149,151]

    if number_trial in primes_list:
	return 1
    else:
	return -1

#####################################################################################################################################################
#####################################################################################################################################################
def output_sequence_seed_permutation(input_integer):

    if primetest(input_integer)==-1:
	print 'The given number is not a prime. Exiting...'
	return -1

    if input_integer==2:
	return [1]
    if input_integer==3:
	return [1,2]


    for test_number in range(2,input_integer):

	initial_element=1
	output_sequence_seed_permutation=[initial_element]
	temp_number=initial_element

	for k in range(input_integer-2):
	    temp_number=(temp_number*test_number)%input_integer
	    if temp_number in output_sequence_seed_permutation:continue
	    else: output_sequence_seed_permutation.append(temp_number)

	if len(output_sequence_seed_permutation)==(input_integer-1):
	    return output_sequence_seed_permutation


    print 'Could not find any seed permutation. Exiting.'
    return -1


#####################################################################################################################################################
#####################################################################################################################################################

input_prime=3
aaa=make_protein_sequences(input_prime)
testing_function(range(input_prime),input_prime+1,aaa)





