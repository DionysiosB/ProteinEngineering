import random

def differences_dfs(input_list,cumulative_list,allowed_numbers,maxint):

    #print '-----------------------------------------------------------------------------------------'
    #print 'Function inputs:'
    #print 'Input_list: ', input_list
    #print 'Cumulative_list:', cumulative_list
    #print 'Allowed_numbers:', allowed_numbers

    if len(input_list)==maxint-2:
	return input_list


    last_element_cumulative_list=cumulative_list[-1]
    #print 'last_element_cumulative_list: ', last_element_cumulative_list
    possible_numbers=list()

    for test_number in allowed_numbers:
	if ((test_number+last_element_cumulative_list) %maxint) !=0:
	    if ((test_number+last_element_cumulative_list) %maxint) not in cumulative_list:
		possible_numbers.append(test_number)

    #print 'possible_numbers:', possible_numbers


    for trial_number in possible_numbers:

	output_list=input_list[:]
	output_list.append(trial_number)

	new_cumulative_list=cumulative_list[:]
	new_cumulative_list.append(((trial_number+last_element_cumulative_list) %maxint))
	#print 'New cumulative list element: ', (trial_number+last_element_cumulative_list)

	new_allowed_numbers=allowed_numbers[:]
	new_allowed_numbers.remove(trial_number)

	#print 'Output list:', output_list
	#print 'New cumulative list:', new_cumulative_list
	#print 'New allowed_numbers, ', new_allowed_numbers


	result= differences_dfs(output_list,new_cumulative_list,new_allowed_numbers,maxint)

	if result != -1:
	    return result


    return -1

#####################################################################################################################################################
#####################################################################################################################################################

test_integer=101
initial_limits=test_integer/2
initial_limits=range(-initial_limits, initial_limits+1)
initial_limits.remove(0)
initial_limits.remove(1)
random.shuffle(initial_limits)
#print initial_limits
test=differences_dfs([1],[1],initial_limits,test_integer)
print test


