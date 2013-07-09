import itertools
import random
import sys



######################################################################################################################################################
##########################		    "REVERSE BRUTE FORCE" ALGORITHM FOR ROTAMER SEQUENCES    #######################################
##########################				DIONYSIOS BARMPOUTIS		   #######################################
##########################				 MENTOR:CHRIS SNOW		       #######################################
######################################################################################################################################################
######################################################################################################################################################





######################################################################################################################################################
######################################################################################################################################################

def master_dict_builder(input_rotamers,number_of_blocks):

    output_dictionary=dict()
    blocks_enumeration_vector=range(number_of_blocks)
    block_combinations_iterator=itertools.combinations(blocks_enumeration_vector, 2)

    if type(input_rotamers)==int:
	number_of_rotamers=input_rotamers
	rotamer_enumeration_vector=range(number_of_rotamers)
	rotamer_product_iterator=itertools.product(rotamer_enumeration_vector,repeat=2)

	all_rotamer_pairs=list()
	for temp_rotamer_product in rotamer_product_iterator:
	    all_rotamer_pairs.append(temp_rotamer_product)
	for temp_block_combination in block_combinations_iterator:
	    output_dictionary[temp_block_combination]=all_rotamer_pairs[:]

    if type(input_rotamers)==list:

	if len(input_rotamers) != number_of_blocks:
	    print 'Rotamers list given, should be of length equal to the number of blocks'
	    return -1
	if type(input_rotamers[0])==int:
	    input_rotamer_lists=list()
	    for k in range(number_of_blocks):
		input_rotamer_lists.append(range(input_rotamers[k]))
	elif type(input_rotamers[0])==list:
	    input_rotamer_lists=input_rotamers[:]
	else:
	    print 'The type of the elements of the list is not recognized, please try again.Exiting...'
	    return -1
	#print input_rotamer_lists
	#print '-------------'

	for temp_block_combination in block_combinations_iterator:

	    rotamer_product_iterator=itertools.product(input_rotamer_lists[temp_block_combination[0]],input_rotamer_lists[temp_block_combination[1]])
	    current_rotamer_pairs=list()
	    for temp_rotamer_product in rotamer_product_iterator:
		current_rotamer_pairs.append(temp_rotamer_product)
	    output_dictionary[temp_block_combination]=current_rotamer_pairs[:]



    return output_dictionary

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

def sequences_list_generator(input_rotamers,number_of_blocks):

    master_dictionary=master_dict_builder(input_rotamers,number_of_blocks)
    output_sequence_list=list()

    finished_flag=0

    while finished_flag==0:


	new_sequence=['X']*number_of_blocks
	#print 'Just made a new sequence:', new_sequence
	finished_flag=1

	first_block_range=range(number_of_blocks-1)
	random.shuffle(first_block_range)

	for first_block in first_block_range:

	    second_block_range=range(first_block+1,number_of_blocks)
	    random.shuffle(second_block_range)

	    for second_block in second_block_range:
		current_block_pair=(first_block,second_block)

		specific_blockpair_remaining_rotamer_pairs=master_dictionary[current_block_pair]
		random.shuffle(specific_blockpair_remaining_rotamer_pairs)

		for remaining_rotamer_pair in specific_blockpair_remaining_rotamer_pairs:

		    finished_flag=0

		    if check_feasibility_sequence_enrichment(current_block_pair,remaining_rotamer_pair,new_sequence)==1:
			new_sequence[first_block]=remaining_rotamer_pair[0]
			new_sequence[second_block]=remaining_rotamer_pair[1]
			master_dictionary[current_block_pair].remove(remaining_rotamer_pair)
	output_sequence_list.append(new_sequence)
    output_sequence_list.pop()

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


#aaa=monte_carlo_sequences_list_generator([ [1,4,7,9],['a','d','n'],['t'],['q','w','r'],['$','@']   ],5,20)
#aaa=monte_carlo_sequences_list_generator([ ['a','b'],['c','d'],['f','g'] ],3,10)
aaa=monte_carlo_sequences_list_generator(3,8,3000)
for k in aaa:
    print k
print 'LENGTH=', len(aaa)



