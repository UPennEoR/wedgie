"""
This piece of code takes antenna position information from hsa7458_v001 and calculates all of the unique baselines.
If an antenna does not have position data, the z position in hsa7458_v001 must be set to -1 so that it is not used in this program.

The output of the program is a dictionary.  The keys of the dictionary are unique baseline lengths.
The value of the dictionary is a list of tuples. Each tuples is a pair of antenna numbers. For instance:

14.6: [ (9, 20), (9, 22) ... ]
...
...


14.6 (a float) is the baseline length. (9, 20) and (9, 22) are pairs of antennae that have a baseline length of 14.6.


Author: Austin Fox Fortino <fortino@sas.upenn.edu>
Date: June 21, 2017
"""

from pprint import pprint


def get_baselines():
	import decimal

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

	decimal.getcontext().prec = 6
	lengths = {}

	for antenna_i in antennae:
		if antennae[antenna_i]['top_z'] == -1:
			continue

		for antenna_j in antennae:
			if antennae[antenna_j]['top_z'] == -1:
				continue
			if antenna_i == antenna_j:
				continue
			
			#Lines 77-80 calculate the baseline between antenna_i and antenna_j.
			dx = antennae[antenna_i]['top_x'] - antennae[antenna_j]['top_x']
			dy = antennae[antenna_i]['top_y'] - antennae[antenna_j]['top_y']
			baseline_length = (decimal.Decimal(dx)**2 + decimal.Decimal(dy)**2).sqrt()
			baseline_length = float(baseline_length)

			#Lines 82-85 makes sure that for a pair of antennae, (a, b), a is less than b.
			if antenna_i < antenna_j:
				pair = (antenna_i, antenna_j)
			else:
				pair = (antenna_j, antenna_i)

			#If baseline_length is a new unqiue baseline, then add it to the dictionary along with its pair.
			if (baseline_length not in lengths):
				lengths[baseline_length] = [pair]

			#If (a, b) or (b, a) already exists, do not add it to the dictionary.  Else, add it to the corresponding existing unique baseline.
			if ((antenna_i, antenna_j) in lengths[baseline_length]) or ((antenna_j, antenna_i) in lengths[baseline_length]):
				continue
			else:
				lengths[baseline_length].append(pair)

	return lengths

#Lines 100-103 Pretty print the dictionary with the list of pairs of antennae sorted.  This is for convenience when testing.
lengths = get_baselines()
for baseline in lengths:
	lengths[baseline].sort()
pprint(lengths)

#Lines 105-108 count how many unique baselines there are.
total_baselines = 0
for key in lengths:
	total_baselines += 1
print "Total Baselines:", total_baselines
