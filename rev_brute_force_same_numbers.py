import itertools
import random
import sys

######################################################################################################################################################
##########################		"REVERSE BRUTE FORCE" ALGORITHM FOR ROTAMER SEQUENCES		 #########################################################
##########################				DIONYSIOS BARMPOUTIS		                         #########################################################
######################################################################################################################################################
######################################################################################################################################################


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
		    print 'FIRST POSITION:', first_position, '	 ', 'SECOND POSITION: ', second_position
		    print 'Combination ', pair_combination , 'does not exist in the existing lists. Exiting...'
		    return -1


    print 'Verification complete. All possible combinations have been found.'
    return 0

######################################################################################################################################################
######################################################################################################################################################

def master_dict_builder(number_of_rotamers,number_of_blocks):

    output_dictionary=dict()
    blocks_enumeration_vector=range(number_of_blocks)
    rotamer_enumeration_vector=range(number_of_rotamers)

    block_combinations_iterator=itertools.combinations(blocks_enumeration_vector, 2)
    rotamer_product_iterator=itertools.product(rotamer_enumeration_vector,repeat=2)

    all_rotamer_pairs=list()
    for temp_product in rotamer_product_iterator:
	all_rotamer_pairs.append(temp_product)

    for temp_combination in block_combinations_iterator:

	output_dictionary[temp_combination]=all_rotamer_pairs[:]

    return output_dictionary

######################################################################################################################################################
######################################################################################################################################################

def master_dictionary_partition(input_master_dictionary):

    consecutive_block_dict=dict()
    nonconsecutive_block_dict=dict()

    for temp_pair in input_master_dictionary:

	if temp_pair[1]==temp_pair[0]+1:
	    consecutive_block_dict[temp_pair]=input_master_dictionary[temp_pair]
	else:
	    nonconsecutive_block_dict[temp_pair]=input_master_dictionary[temp_pair]

    output_list=[consecutive_block_dict, nonconsecutive_block_dict]
    return output_list

######################################################################################################################################################
######################################################################################################################################################

def eligible_pairs_list(input_list,first_element):

    output_list=list()
    for dummy_tuple in input_list:
	if dummy_tuple[0]==first_element:
	    output_list.append(dummy_tuple)
    return output_list



######################################################################################################################################################
######################################################################################################################################################
def check_feasibility_sequence_enrichment(block_pair,rotamer_pair,current_list):

    output=1

    if ( (current_list[block_pair[0]] !='X') and (current_list[block_pair[0]] !=rotamer_pair[0]) ):
	output=0
    if ( (current_list[block_pair[1]] !='X') and (current_list[block_pair[1]] !=rotamer_pair[1]) ):
	output=0


    return output


######################################################################################################################################################
######################################################################################################################################################

def sequences_list_generator(number_of_rotamers,number_of_blocks):

    input_master_dictionary=master_dict_builder(number_of_rotamers,number_of_blocks)
    input_master_dictionary_partition=master_dictionary_partition(input_master_dictionary)

    consecutive_block_pairs=input_master_dictionary_partition[0]
    nonconsecutive_block_pairs=input_master_dictionary_partition[1]


    rotamer_enumeration_vector=range(number_of_rotamers)
    rotamer_pair_iterator=itertools.product(rotamer_enumeration_vector,repeat=2)


    output_sequence_list=list()

    for first_pair in rotamer_pair_iterator:
	new_block_sequence=list()
	new_block_sequence.extend(list(first_pair))

	consecutive_block_pairs[(0,1)].remove(first_pair)

	for block_number in range(1,number_of_blocks-1):
	    current_block_pair=(block_number,(block_number+1))
	    current_last_block=new_block_sequence[-1]
	    current_available_rotamer_pairs=consecutive_block_pairs[current_block_pair]
	    eligible_rotamer_pairs=eligible_pairs_list(current_available_rotamer_pairs,current_last_block)
	    next_pair=random.choice(eligible_rotamer_pairs)
	    new_block_sequence.append(next_pair[1])

	    consecutive_block_pairs[current_block_pair].remove(next_pair)

	output_sequence_list.append(new_block_sequence)


	for first_block in range(number_of_blocks-1):
	    for second_block in range(first_block+2,number_of_blocks):
		nonneighboring_block_pair=(first_block,second_block)
		current_rotamer_pair=(new_block_sequence[first_block],new_block_sequence[second_block])
		if current_rotamer_pair in nonconsecutive_block_pairs[nonneighboring_block_pair]:
		    nonconsecutive_block_pairs[nonneighboring_block_pair].remove(current_rotamer_pair)


    finished_flag=0

    while (finished_flag==0):

	finished_flag=1

	new_sequence=['X']*number_of_blocks


	first_block_range=range(number_of_blocks-1)
	random.shuffle(first_block_range)
	for first_block in first_block_range:
	    second_block_range=range(first_block+2,number_of_blocks)
	    random.shuffle(second_block_range)
	    for second_block in second_block_range:
		nonneighboring_block_pair=(first_block,second_block)
		specific_blockpair_remaining_rotamer_pairs=nonconsecutive_block_pairs[nonneighboring_block_pair]
		random.shuffle(specific_blockpair_remaining_rotamer_pairs)
		for remaining_rotamer_pair in specific_blockpair_remaining_rotamer_pairs:
		    if check_feasibility_sequence_enrichment(nonneighboring_block_pair,remaining_rotamer_pair,new_sequence)==1:

			new_sequence[first_block]=remaining_rotamer_pair[0]
			new_sequence[second_block]=remaining_rotamer_pair[1]
			nonconsecutive_block_pairs[nonneighboring_block_pair].remove(remaining_rotamer_pair)

			finished_flag=0


	output_sequence_list.append(new_sequence)
    dummy=output_sequence_list.pop()

    return output_sequence_list

######################################################################################################################################################
######################################################################################################################################################
def monte_carlo_sequences_list_generator(input_rotamers,number_of_blocks, number_of_trials):

    final_sequences_list=[]
    final_sequences_list_length=sys.maxint

    for trial_number in range(number_of_trials):
	#print 'Trial number:', trial_number
	trial_list=sequences_list_generator(input_rotamers,number_of_blocks)
	trial_length=len(trial_list)
	if trial_length <= final_sequences_list_length:
	    final_sequences_list=trial_list
	    final_sequences_list_length=trial_length
	    print 'UPDATED: NEW NUMBER OF SEQUENCES:',final_sequences_list_length

    return final_sequences_list


######################################################################################################################################################
######################################################################################################################################################


length=8
alternatives=3
aaa= monte_carlo_sequences_list_generator(3,8,1500)
testing_function(range(alternatives),length, aaa)
print 'LENGTH=', len(aaa)


