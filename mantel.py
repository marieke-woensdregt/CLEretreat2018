#!/usr/bin/env python

from sys import argv
from scipy import mean, random, spatial, stats, std

############ INSERT read_file() ############



############ END OF read_file() ############


############ INSERT pairwise_distances() ############


def pairwise_distances(listofstrings):
	n = len(listofstrings)
	pwdistlist = []
	for i in range(0, n-1):
		for j in range(i+1, n):
			dist = levenshtein_distance(listofstrings[i], listofstrings[j])
			pwdistlist.append(dist)
	return pwdistlist


############ END OF pairwise_distances() ############


############ INSERT monte_carlo() ############



############ END OF monte_carlo() ############

def levenshtein_distance(s1, s2):
	'''
	Takes two stirngs and returns the normalized Levenshtein distance
	'''
	if len(s1) > len(s2):
		s1,s2 = s2,s1
	distances = range(len(s1) + 1)
	for index2, char2 in enumerate(s2):
		newDistances = [index2 + 1]
		for index1, char1 in enumerate(s1):
			if char1 == char2:
				newDistances.append(distances[index1])
			else:
				newDistances.append(1 + min((distances[index1], distances[index1+1], newDistances[-1])))
		distances = newDistances
	return float(distances[-1])/max(len(s1), len(s2))

def shuffle_distances(distances):
	'''
	Takes a list of pairwise distances, converts it to a distance matrix,
	shuffles the matrix, and returns the upper triangle as a vector.
	'''
	matrix = spatial.distance.squareform(distances, 'tomatrix')
	shuffled_vector = []
	n = len(matrix)
	shuffle_order = range(0, n)
	random.shuffle(shuffle_order)
	c = 0
	for i in range(0, n-1):
		for j in range(i+1, n):
			shuffled_vector.append(matrix[shuffle_order[i]][shuffle_order[j]])
			c += 1
	return shuffled_vector

def mantel_test(distances1, distances2, randomizations):
	'''
	Takes two lists of pairwise distances and performs a Mantel test. Returns
	the veridical correlation (r), the mean (m) and standard deviation (sd)
	of the Monte Carlo sample correlations, and a Z-score (z) quantifying the
	significance of the veridical correlation.
	'''
	r = stats.pearsonr(distances1, distances2)[0]
	m, sd = monte_carlo(distances1, distances2, randomizations)
	z = (r-m)/sd
	return r, m, sd, z

def run_mantel(filename):
	'''
	Reads in file and runs Mantel Test on the data
	'''
	strings, meanings = read_file(filename)
	pairwise_dist_strings = pairwise_distances(strings)
	pairwise_dist_meanings = pairwise_distances(meanings)
	r, m, sd, z = mantel_test(pairwise_dist_strings, pairwise_dist_meanings, 10000)
	return r, m, sd, z

if __name__ == '__main__':
	print(run_mantel(argv[1]))
