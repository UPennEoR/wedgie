import numpy as np

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
	print("Total Unique Baselines:", total_baselines)
	print("Total Pairs:", total_pairs)
	print("Bad Antennae:", bad_ants)

def RMSE(fit_function, real_function):
    """Finds and return the RMSE between a function and the fit of that function."""
    residuals = fit_function - real_function
    squared_residuals = residuals ** 2
    mean = np.mean(squared_residuals)
    RMSE = np.sqrt(mean)
    return RMSE