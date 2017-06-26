"""
This program contains an importable function get_baselines().  This function takes 
HERA position information from hsa7458_v001.py.  If there is no position data for an 
antenna in hsa7458_v001.py (i.e., the antenna does not yet exist) the z-coordinate 
should be set to -1 in hsa7458_v001.py.  This will make sure that those antennae 
are not used in this program.

If hsa7458_v001.py is not imported, hardcoded position information will 
be used for 19 antennae.

This position information is used to calculate all unique baselines, and 
find all pairs of antennae with each baseline.  The output of this program is a 
Python dictionary. Each key is a unique baseline, and each value is a list of tuples.  
Each tuple contains two numbers that represent antennae.
Example:
{14.6: [(9, 20),(9, 22),(9, 53), ...
...
...

Antennae 9 and 20, 9 and 22, and 9 and 53 are separated by 14.6 meters.
All baselines are stored as floats. All antennae are stored as ints.

When some antennae aren't functioning properly the command line argument --ex_ants or -x can be used.
Example:
python2.7 baseline_lengths.py --ex_ants=80,104,96
python2.7 baseline_lengths.py --x=80,104,96

Both of those commands will exclude antennae 80, 104, and 96 when calculating baselines.

Author: Austin Fox Fortino <fortino@sas.upenn.edu>
Date: June 21, 2017
Last Updated: June 23, 2017
"""
from pprint import pprint

# This block adds the ability for the program to not include certain antennae.
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-x", "--ex_ants", type=str)
args = parser.parse_args()
bad_ants = args.ex_ants

if bad_ants != None: bad_ants = [int(ant) for ant in bad_ants.split(",")]
else: bad_ants = []


def calculate_baseline(antennae, pair):
	"""
	The decimal module is necessary for keeping the number of decimal places small.
	Due to small imprecision, if more than 8 or 9 decimal places are used, 
	many baselines will be calculated that are within ~1 nanometer to ~1 picometer of each other.
	Because HERA's position precision is down to the centimeter, there is no 
	need to worry about smaller imprecision.
	"""
	from decimal import getcontext, Decimal
	getcontext().prec = 6

	dx = antennae[pair[0]]['top_x'] - antennae[pair[1]]['top_x']
	dy = antennae[pair[0]]['top_y'] - antennae[pair[1]]['top_y']
	baseline = float((Decimal(dx)**2 + Decimal(dy)**2).sqrt())
	
	return baseline


def get_baselines():

	# This try-except block is necessary for testing purposes when the hsa7458_v001.py file is unavailable.
	try:
		import hsa7458_v001 as cal
		antennae = cal.prms['antpos_ideal']
		print "Successfully imported hsa7458_v001."
	except Exception:
		print "Unsuccessfully imported hsa7458_v001. Using hardcoded antenna position information."
		print "As of June 22, 2017, this hardcoded information is up to date and comprehensive for 19 antennae."
		antennae = {
			80  :{'top_x':-14.6,	'top_y':-25.28794179,	'top_z':0.0 },	#a0
			104 :{'top_x':0.0,		'top_y':-25.28794179,	'top_z':0.0 },	#a1
			96  :{'top_x':14.6,		'top_y':-25.28794179,	'top_z':0.0 },	#a2

			64  :{'top_x':-21.9,	'top_y':-12.6439709,	'top_z':0.0 }, 	#a3
			53  :{'top_x':-7.3,		'top_y':-12.6439709,	'top_z':0.0 }, 	#a4
			31  :{'top_x':7.3,		'top_y':-12.6439709,	'top_z':0.0 }, 	#a5
			65  :{'top_x':21.9,		'top_y':-12.6439709,	'top_z':0.0 }, 	#a5


			88  :{'top_x':-29.2,	'top_y':0.0,			'top_z':0.0 }, 	#a5
			9   :{'top_x':-14.6,	'top_y':0.0,			'top_z':0.0 }, 	#a5
			20  :{'top_x':0.0,		'top_y':0.0,			'top_z':0.0 }, 	#a5
			89  :{'top_x':14.6,		'top_y':0.0,			'top_z':0.0 }, 	#a5
			43  :{'top_x':29.2,		'top_y':0.0,			'top_z':0.0 }, 	#a5

			105 :{'top_x':-21.9,	'top_y':12.6439709,		'top_z':0.0 }, 	#a3
			22  :{'top_x':-7.3,		'top_y':12.6439709,		'top_z':0.0 }, 	#a4
			81  :{'top_x':7.3,		'top_y':12.6439709,		'top_z':0.0 }, 	#a5
			10  :{'top_x':21.9,		'top_y':12.6439709,		'top_z':0.0 }, 	#a5


			72  :{'top_x':-14.6,	'top_y':25.28794179,	'top_z':0.0 }, 	#a0
			112 :{'top_x':0.0,		'top_y':25.28794179,	'top_z':0.0 }, 	#a1
			97  :{'top_x':14.6,		'top_y':25.28794179,	'top_z':0.0 }, 	#a2
		}

	baselines = {}

	# This block contains the algorithm that determines the baseline and places them in the dictionary.
	for antenna_i in antennae:
		if antennae[antenna_i]['top_z'] == -1: continue
		if antenna_i in bad_ants: continue
		
		for antenna_j in antennae:
			if antennae[antenna_j]['top_z'] == -1: continue
			if antenna_j in bad_ants: continue

			if antenna_i == antenna_j: continue
			elif antenna_i < antenna_j: pair = (antenna_i, antenna_j)
			elif antenna_i > antenna_j: pair = (antenna_j, antenna_i)

			baseline = calculate_baseline(antennae, pair)

			if (baseline not in baselines): baselines[baseline] = [pair]
			elif (pair in baselines[baseline]): continue
			else: baselines[baseline].append(pair)

	return baselines

# This function runs and displays information useful for troubleshooting.
def display():
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