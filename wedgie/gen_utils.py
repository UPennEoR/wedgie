def findMiddle(input_list):
    """
    This method returns the central value of a list.
    If there is an even number of values, it returns
    the average of the central two. 
    
    i.e. this should be used for constantly-spaced-value lists.
    
    It's kinda stupid.
    """
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return input_list[int(middle - .5)]
    else:
        return np.mean([input_list[int(middle)], input_list[int(middle-1)]])

def display_baseline_length_dict():
    # This function runs and displays information useful for troubleshooting calfile stuff.
	total_baselines = 0
	total_pairs = 0
	baselines = get_baselines()
	for baseline in baselines:
		baselines[baseline].sort()

		total_baselines += 1

		for pairs in baselines[baseline]:
			total_pairs += 1

	pprint(baselines)
	print "Total Unique Baselines:", total_baselines
	print "Total Pairs:", total_pairs
	print "Bad Antennae:", bad_ants