# This function runs and displays information useful for troubleshooting.
def display_baseline_length_dict():
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

display()