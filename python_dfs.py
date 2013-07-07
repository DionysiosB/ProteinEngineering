def dfs_function(input_list,allowed_differences,allowed_numbers,maxint):

    #print '----------------------------------------------------------------'
    #print 'Current input list:', input_list
    #print 'Current allowed_differences:', allowed_differences
    #print 'Current allowed numbers:', allowed_numbers

    current_allowed_numbers=list()

    if len(allowed_differences)==1:
	if len(allowed_numbers)>0:
	    #print "No allowed differences. Exiting"
	    return -1
	else:
	    return input_list

    last_list_element=input_list[-1]
    #print 'Last list element:', last_list_element
    for temp_number in allowed_numbers:
	temp_difference=(temp_number-last_list_element)%maxint
	if temp_difference in allowed_differences:
	    current_allowed_numbers.append(temp_number)

    if len(current_allowed_numbers)==0:
	return -1

    for number_trial in current_allowed_numbers:
	new_difference=(number_trial-last_list_element)%maxint

	output_list=input_list[:]
	output_list.append(number_trial)

	output_allowed_numbers=allowed_numbers[:]
	output_allowed_numbers.remove(number_trial)

	output_allowed_differences=allowed_differences[:]
	output_allowed_differences.remove(new_difference)

	result=dfs_function(output_list,output_allowed_differences,output_allowed_numbers,maxint)

	if result !=-1:
	    return result

    return -1

#####################################################################################################################################################
#####################################################################################################################################################

test_integer=31
test=dfs_function([1],range(1,test_integer),range(2,test_integer),test_integer)
print test







