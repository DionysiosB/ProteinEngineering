import random

def structured_differences(input_integer):

    if ((input_integer<=2)  | (input_integer%2==0)):
	print 'The input integer should be an odd number greater than zero.Exiting...'
	return -1


    limit=input_integer/2

    differences_vector=[1,-1]

    if input_integer==3:
	return differences_vector




    current_limit=1

    while current_limit<limit:

	current_limit+=1
	print 'Current_limit --->', current_limit
	differences_vector=list_pair_placement(differences_vector,current_limit,2*current_limit+1)


    return differences_vector


#########################################################################################################################################################################################
#########################################################################################################################################################################################


def list_cumulative_sum(input_list):

    output_list=list()
    temp=0
    for k in input_list:
	temp+=k
	output_list.append(temp)

    return output_list


#########################################################################################################################################################################################
#########################################################################################################################################################################################
def list_modulo_operation(input_list, modulo_number):

    output_list=list()
    for temp in input_list:
	output_list.append(temp % modulo_number)
    return output_list



#########################################################################################################################################################################################
#########################################################################################################################################################################################
def condition_checker(input_list,modulo_integer):

    temp_list=input_list[:]
    #print temp_list

    list_length=len(input_list)

    for current_start in range(list_length-1):
	current_checking_list=temp_list[current_start:(list_length-1)]
	current_checking_list=list_cumulative_sum(current_checking_list)
	#print current_checking_list
	resulting_list=list_modulo_operation(current_checking_list,modulo_integer)
	if 0 in resulting_list:
	    return -1


    return 1


#########################################################################################################################################################################################
#########################################################################################################################################################################################
def list_pair_placement(input_list,limit_integer,modulo_integer):


    list_length=len(input_list)
    for position_1 in range(len(input_list)):
	for position_2 in range(len(input_list)+2,0,-1):

	    output_list=input_list[:]
	    output_list.insert(position_1,-limit_integer)
	    output_list.insert(position_2,limit_integer)
	    print 'Mock Vector --->', output_list

	    result=condition_checker(output_list,modulo_integer)
	    if result==1:
		print 'Passed test --->', output_list
		return output_list

    for position_1 in range(len(input_list)):
	for position_2 in range(2,len(input_list)+2):

	    output_list=input_list[:]
	    output_list.insert(position_1,limit_integer)
	    output_list.insert(position_2,-limit_integer)
	    #print 'Mock Vector --->', output_list

	    result=condition_checker(output_list,modulo_integer)
	    if result==1:
		print 'Passed test --->', output_list
		return output_list

    print 'Could not find ANY vector... Exit status -1'
    return -1




#########################################################################################################################################################################################
#########################################################################################################################################################################################

def structured_differences_2(input_integer):

    if ((input_integer<=2)	| (input_integer%2==0)):
	print 'The input integer should be an odd number greater than zero.Exiting...'
	return -1


    limit=input_integer/2

    output_list=range(1,limit+1)
    negative_list=range(-1,-limit-1,-1)
    output_list.extend(negative_list)
    random.shuffle(output_list)

    print output_list
    list_length=len(output_list)

    differences_list_working=-1

    counter =0
    while differences_list_working==-1:

	counter+=1
	if counter%10000==0:
	    print counter
	if counter==500000:
	    break

	differences_list_working=1
	temp_list=output_list[:]

	for current_start in range(list_length-1):
	    current_checking_list=temp_list[current_start:(list_length-1)]
	    current_checking_list=list_cumulative_sum(current_checking_list)
	    #print 'Current Checking list: ', current_checking_list
	    resulting_list=list_modulo_operation(current_checking_list,input_integer)
	    #print 'Resulting list: ', resulting_list
	    if 0 in resulting_list:
		differences_list_working=-1
		#print 'Current_start:', current_start
		temp_element=output_list.pop(current_start)
		random_position=random.randint(1,list_length)
		output_list.insert(random_position,temp_element)
		#print 'New output list:', output_list
		break;


    if differences_list_working==1:
	print 'FOUND IT!!!'
    else:
	print 'COULD NOT FIND IT'

    return output_list






#########################################################################################################################################################################################
#########################################################################################################################################################################################



aaa=structured_differences_2(35)
print aaa
#print sorted(aaa)















