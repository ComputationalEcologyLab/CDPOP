# -------------------------------------------------------------------------------------------------
# CDPOP_Modules.py
# Author: Erin L Landguth
# Created: June 2010
# Description: This is the function/module file for CDPOP v0.88
# --------------------------------------------------------------------------------------------------

# Import Modules with Except/Try statements

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."

# CDPOP functions
try:
	from CDPOP_PreProcess import *
except ImportError:
	raise ImportError, "CDPOP PreProcess required."
	
# Python specific functions
import os, random, copy, pdb, sys
from collections import Counter

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False


# --------------------------------------------------------------------------
def PrepTextFile(textpath):
	'''
	PrepTextFile() - Prepare the input files
	'''
	
	return textpath.strip('\n').strip('\r')
	
	# End::PrepTextFile()

# --------------------------------------------------------------------------
def logMsg(outf,msg):
	'''
	logMsg() --log file message handler.
	Inputs:
	outf - open file handle
	msg -- string containing formatted message
	--always outputs to log file by default.
	--using msgVerbose, can be set to "Tee" output to stdout as well
	'''
	outf.write(msg+ '\n')
	if msgVerbose:
		print("%s"%(msg))
		
	# End::logMsg()

# ---------------------------------------------------------------------------------------------------	 
def w_choice_general(lst):
	'''
	w_choice_general()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(x[1] for x in lst)
	n=random.uniform(0,wtotal)
	count = 0
	for item, weight in lst:
		if n < weight:
			break
		n = n-weight
		count = count + 1
	return item,count
	
	#End::w_choice_general()
	
# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ---------------------------------------------------------------------------------------------------	 
def ReadGrid(FIDnew,idnew,agenew,xgridnew,ygridnew,\
genesnew,equalsexratio,sexnew,subpopnew,infectionnew,allelst,geneswap,gen,intgenesans):	
	
	'''
	DoReadGrid()
	This function is reads the previous generations
	grid information at the start of the time loop.
	'''
	
	FID = FIDnew
	subpop = subpopnew
	xgrid = xgridnew
	xgridcopy = copy.deepcopy(xgridnew)
	ygrid = ygridnew
	ygridcopy = copy.deepcopy(ygridnew)
	id = idnew
	age = agenew	
	infection = infectionnew
	
	# Store the number of grids
	nogrids = len(FID)
	
	# Get the number of actual filled grids and np.where indexing
	#	different in windows vs. linux: could be version
	if len(np.where(np.asarray(age)=='NA')) == 0:
		filledgrids = nogrids
	else:
		filledgrids = nogrids-len(np.where(np.asarray(age)=='NA')[0])
	
	# If gen is < geneswap, do nothing, else intialize genes here
	if intgenesans != 'known' and gen == geneswap:
	
		# Loop through all grids
		genes = []
		for i in xrange(nogrids):
									
			# And store genes information
			genes.append([])
						
			# For each loci:
			for j in xrange(len(allelst[0])):			
				# Take a random draw from the w_choice function at jth locus
				if len(allelst) > 1:
					rand1 = w_choice_general(allelst[int(subpop[i])-1][j])[0]
					rand2 = w_choice_general(allelst[int(subpop[i])-1][j])[0]
				else:					
					rand1 = w_choice_general(allelst[0][j])[0]
					rand2 = w_choice_general(allelst[0][j])[0]

				# Store genes loci spot
				genes[i].append([])
				
				# Append assinment onto indall array - run through each condition for assignment of 1s or 2s or 0s
				# 	1s = heterozygous at that locus
				#	2s = homozygous at that locus
				#	0s = absence of allele
				for k in xrange(len(allelst[0][j])):
					
					# Somebody not in this spot
					if age[i] == 'NA':
						# And to genes list
						genes[i][j].append('NA')
						
					# Else if somebody is in spot, assign genes
					else:					
						# Assignment of 2, the rest 0
						if rand1 == rand2: 
							if k < rand1 or k > rand1:
								tempindall = 0
							elif k == rand1:
								tempindall = 2
								
						# Assignment of 1s, the rest 0
						if rand1 != rand2:
							if k < min(rand1,rand2) or k > max(rand1,rand2):
								tempindall = 0
							elif k == rand1 or k == rand2:
								tempindall = 1
							else:
								tempindall = 0
								
						# And to genes list
						genes[i][j].append(tempindall)
			
	else:
		genes = genesnew
	
	# If sex is not equal ratio
	if equalsexratio == 'N' or equalsexratio == 'AtBirth':
		sex = sexnew
	
	# If equal sex ratio is Y, then split up sex into equal parts
	if equalsexratio == 'WrightFisher':
		sex = []
		# Get unique number of subpops, make equal sex ratio for each subpopulation
		nosubpops = len(np.unique(subpop))
		unique_subpops = np.unique(subpop)
		
		# Count up the unique number of subgrids appending to subgrids
		subgridtotal = []
		# Create list of lists storage spots for number of subgrids
		for i in xrange(nosubpops):
			subgridtotal.append([])
		for i in xrange(len(subpop)):
			# Loop through unique subpops
			for j in xrange(nosubpops):
				# If subpop exits append to subgrid spot
				if subpop[i] == unique_subpops[j]:
					subgridtotal[int(unique_subpops[j])-1].append(1)
		# And then sum them up
		for i in xrange(nosubpops):
			subgridtotal[i] = sum(subgridtotal[i])
			# If the subpopulation number is not even then sys exit
			if np.mod(subgridtotal[i],2) == 1:
				print("You have equal sex ratio turned on and this population is not even.")
				sys.exit(-1)
			
		# Then loop through each subpopulation
		for i in xrange(nosubpops):			
			# Then create half males and females and shuffle
			sextemp = np.append(np.zeros(subgridtotal[i]/2,"int"),np.ones(subgridtotal[i]/2,"int"))
			np.random.shuffle(sextemp)
			# Loop through these individuals and append to LIST sex
			for j in xrange(len(sextemp)):			
				# The add them together and shuffle and append to sex
				sex.append(str(sextemp[j]))
		
		# Delete extra stuff
		del(sextemp)
	
	#Return this functions variables
	tupReadGrid = FID,sex,id,age,xgrid,xgridcopy,ygrid,\
	ygridcopy,genes,nogrids,subpop,infection,filledgrids
	return tupReadGrid
	
	#End::DoReadGrid()

# ---------------------------------------------------------------------------------------------------	 
def GetMetrics(Population,nogrids,loci,alleles,genes,gen,Ho,\
Alleles,He,subpop,p1,p2,q1,q2,Population_age,age,agemort,geneswap,\
allelst,F=None,FST=None,FIS=None,FIT=None):
	'''
	GetMetrics()
	This function summarizes the genotypes and
	produces genetic metrics.
	Ho - Observed heterozygosity per generation
	He - Expected heterozygoisty per generation
	Alleles - Total number of unique alleles in genotype*individuals
	'''
	
	# Track population age numbers
	Population_age.append([])
	for i in xrange(len(agemort)):
		Population_age[gen].append([])
	countage = Counter(np.asarray(np.asarray(age,dtype='|S10')[np.where(np.asarray(age,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	for i in xrange(1,len(agemort)+1):
		Population_age[gen][i-1].append(countage[i])
	
	# Sum up population age tracker
	for i in xrange(len(agemort)):
		Population_age[gen][i] = sum(Population_age[gen][i])
		
	# List for total, left, and right
	unique_alleles = Alleles
	
	# Get unique number of subpops
	nosubpops = len(np.unique(subpop))
	unique_subpops = np.unique(subpop)
	
	# --------------
	# Eventually clean this up for subpops that are not seq 1 to N
	#subgridtest = [ (i,subpop.count(i)) for i in set(subpop) ]
	#subgridcount = Counter(subpop)
	subgridtotal = []	
	# Create list of lists storage spots for number of subgrids
	for i in xrange(nosubpops):
		subgridtotal.append([])
	for i in xrange(len(subpop)):
		# Loop through unique subpops
		for j in xrange(nosubpops):
			# If subpop exits append to subgrid spot
			if subpop[i] == unique_subpops[j] and age[i] != 'NA':
				subgridtotal[int(unique_subpops[j])-1].append(1)

	# And then sum them up - and store numbers
	Population.append([])
	for i in xrange(nosubpops):
		subgridtotal[i] = sum(subgridtotal[i])
		# Add information to Population tracker
		Population[gen].append(subgridtotal[i])
	# And then get the number of filled grids
	filledgrids = sum(subgridtotal)
	# Add Population total
	Population[gen].insert(0,filledgrids)
	
	# Only complete if greater than startGenes
	if gen >= geneswap:
	
		# Get allele location as seqence from alleles array
		allele_numbers = []
		for i in xrange(loci):
			for j in xrange(alleles[i]):
				allele_numbers.append(j)
		allele_numbers = np.asarray(allele_numbers)
		
		# Remove the 'NA' gene values
		tempgenes = []
		for i in xrange(len(genes)):
			if age[i] != 'NA':
				tempgenes.append(genes[i])
		
		# Cast genes as an numpy array as byte type
		genes_array_woNA = np.asarray(tempgenes,dtype='float')
		genes_array_wNA = np.asarray(genes)
		
		# The total number of alleles
		total_alleles = len(allele_numbers)
		
		# Create a list to store the subpopulation grid number information
		subgrids = []
		all_freq_sub = []
		ho_count_sub = []
		ho_sub = []
		all_freq_sq_sub = []
		homozygosity_sub = []
		he_sub = []
		sumsubpopsHo = []
		sumsubpopsHe = []
		alleles_sub = []
		
		# Create list of lists storage spots for number of subgrids
		for i in xrange(nosubpops):
			subgrids.append([])
			all_freq_sub.append([])
			ho_count_sub.append([])
			ho_sub.append([])
			all_freq_sq_sub.append([])
			homozygosity_sub.append([])
			he_sub.append([])
			alleles_sub.append([])
			
		# Count up the unique number of subgrids appending to subgrids
		for i in xrange(len(subpop)):
			# Loop through unique subpops
			for j in xrange(nosubpops):
				# If subpop exits append to subgrid spot
				if subpop[i] == unique_subpops[j] and age[i] != 'NA':
					subgrids[int(unique_subpops[j])-1].append(i)
			
		# Get allele frequency for total
		if filledgrids != 0:
			all_freq_tot = np.asarray(np.sum(genes_array_woNA,axis=0),dtype = 'float').reshape(total_alleles)
			all_freq_tot = all_freq_tot/(2*filledgrids)
		else:
			all_freq_tot = np.zeros(total_alleles,float)
			
		# Get allele frequency for subpopulations
		for i in xrange(nosubpops):
			if subgridtotal[i] != 0:
				all_freq_sub[i].append(np.asarray(np.sum(np.asarray(genes_array_wNA[subgrids[i],:,:],dtype='float'),axis=0),dtype = 'float').reshape(total_alleles))
				all_freq_sub[i] = all_freq_sub[i][0]/(2*subgridtotal[i])
			else:
				all_freq_sub[i].append(np.zeros(total_alleles,float))
				all_freq_sub[i] = all_freq_sub[i][0]
		
		# Create an array to fill up with allele frequencies - only for total
		all_freq_list = np.zeros((total_alleles,2))		
		all_freq_list[:,0] = allele_numbers
		all_freq_list[:,1] = all_freq_tot	
		
		#Calculate the number of homogenous alleles for total
		ho_count_tot = (np.array(genes_array_woNA==2)).sum()
		# Calculate the number of homogenous alleles in each subpop
		for i in xrange(nosubpops):
			ho_count_sub[i].append((np.array(np.asarray(genes_array_wNA[subgrids[i],:,:],dtype='float')==2)).sum())
		
		# Calculate the observed het for total
		if filledgrids != 0:
			ho_tot = (float(filledgrids*loci - ho_count_tot)/(loci*filledgrids))
		else:
			ho_tot = 0.0		
		# Calculate the observed het in each subpop
		for i in xrange(nosubpops):
			if subgridtotal[i] != 0:
				ho_sub[i].append((float(subgridtotal[i]*loci - ho_count_sub[i][0])/(loci*subgridtotal[i])))
			else:
				ho_sub[i].append(0.0)
			
		# Append Ho information (Observed Het)
		Ho.append([ho_tot])
		for i in xrange(nosubpops):
			Ho[gen].append(ho_sub[i][0])
		
		# Get the sqare of the allele frequency for total
		all_freq_sq_tot = all_freq_tot**2
		
		# Get the square of the allele frequency for subpops
		for i in xrange(nosubpops):
			all_freq_sq_sub[i].append(all_freq_sub[i]**2)
		
		# Calculate the homozygosity for total populations
		homozygosity_tot = sum(all_freq_sq_tot)/loci
		
		# Calculate the homozygosity for subpopulations
		for i in xrange(nosubpops):
			homozygosity_sub[i].append(sum(all_freq_sq_sub[i][0])/loci)
		
		# Store He for [Total]
		if filledgrids != 0:
			he_tot = (1. - homozygosity_tot)
		else:
			he_tot = 0.0
		
		# Store He for subpopulations
		for i in xrange(nosubpops):
			if subgridtotal[i] != 0:
				he_sub[i].append(1. - homozygosity_sub[i][0])
			else:
				he_sub[i].append(0.0)
			
		# Append He information (Expected Het)
		He.append([he_tot])
		for i in xrange(nosubpops):
			He[gen].append(he_sub[i][0])
		
		# Get total number of alleles
		alleles_tot = np.array(all_freq_tot>0.).sum()
		
		# Get the total number of alleles in each subpop
		for i in xrange(nosubpops):
			alleles_sub[i].append(np.array(all_freq_sub[i]>0.).sum())
		
		# Append allele total information
		unique_alleles.append([alleles_tot])
		for i in xrange(nosubpops):
			unique_alleles[gen].append(alleles_sub[i][0])
		
		# Get allele frequency totals for selection section
		p1.append(all_freq_tot[0])
		p2.append(all_freq_tot[1])
		q1.append(all_freq_tot[2])
		q2.append(all_freq_tot[3])

	# Return and store numbers if less than geneswap time
	else:
		p1.append(np.nan)
		p2.append(np.nan)
		q1.append(np.nan)
		q2.append(np.nan)
		unique_alleles.append([np.nan])
		He.append([np.nan])
		Ho.append([np.nan])
		for i in xrange(nosubpops+1):			
			unique_alleles[gen].append(np.nan)
			He[gen].append(np.nan)
			Ho[gen].append(np.nan)
	
	# Return variables from this function
	tupGetMetrics = unique_alleles,filledgrids,subgridtotal
	return tupGetMetrics
	
	#End::GetMetrics()
	
# ---------------------------------------------------------------------------------------------------	 
def InheritGenes(gen,AllelesMutated,offspringno,\
offspring,genes,loci,muterate,mtdna,mutationans,geneswap):
	'''
	InheritGenes()
	Pass along gentic information to survived offspring from parents
	Input: offspring, genes 
	Output: [femaleid,maleid,cdmatidofmother,cdmatidoffather,sex],[
	genetic information]		
	'''		
	
	# Check for generation to start swapping genes
	if gen >= geneswap:
	
		# Storage for tracking how many alleles mutated
		noallelesmutated = []
			
		# If there are offspring
		if int(offspringno) != int(0):
		
			# Begin loop through offspring
			for i in xrange(offspringno):
				
				# Add spot in offspring array for individual i's genes
				offspring[i].append([])
				
				# Temp storage for i's mother's genes
				mothergenes=genes[offspring[i][0]]
				# Temp storage for i's father's genes
				fathergenes=genes[offspring[i][1]]
				
				# Loop through each locus
				for jspot in xrange(loci):
					
					# Temporary index storage
					tempindfather = []
					tempindmother = []
										
					# Loop through each allele-mother and father have same len of alleles
					for kspot in xrange(len(fathergenes[jspot])):
												
						# Check the homogeneous 2 case - in father genes
						if int(fathergenes[jspot][kspot])==2:
							
							# Get a random number for allele mutation
							mutationrandno = rand()
							
							# Check if random number is less than or equal to muterate
							if mutationrandno <= muterate:
							
								# If backward and forward mutation
								if mutationans == 'random':
									# Randomly choose another allele								
									randallelespot = int(len(fathergenes[jspot])*rand())
									while randallelespot == kspot:
										randallelespot = int(len(fathergenes[jspot])*rand())
									tempindfather.append(randallelespot)									
									# Count a mutation
									noallelesmutated.append(1)
									
								# If just forward mutation
								elif mutationans == 'forward':
									if kspot != len(fathergenes[jspot])-1:
										tempindfather.append(kspot+1)									
										# Count a mutation
										noallelesmutated.append(1)
									else:
										# THen no mutation
										tempindfather.append(kspot)

								# If just forward mutation
								elif mutationans == 'backward':
									if kspot != 0:
										tempindfather.append(kspot-1)									
										# Count a mutation
										noallelesmutated.append(1)
									else:
										# THen no mutation
										tempindfather.append(kspot)
										
								# If forward and backward mutation
								elif mutationans == 'forwardbackward':
									# Then random forward or backward step
									randstep = rand()
									# To go left, but it can't be the first allele
									if randstep < 0.5 and kspot != 0:
										tempindfather.append(kspot-1)
										# Count a mutation
										noallelesmutated.append(1)
									# To go right, but it can't be the last allele
									elif randstep >= 0.5 and kspot != len(fathergenes[jspot])-1:
										tempindfather.append(kspot+1)
										# Count a mutation
										noallelesmutated.append(1)
									else:
										# THen no mutation
										tempindfather.append(kspot)
								
								# If forward mutation in A and backward mutation for b (A -> a, b -> B)
								elif mutationans == 'forwardAbackwardBrandomN':
									if jspot == 0 and kspot == 0:
										tempindfather.append(kspot+1)									
										# Count a mutation
										noallelesmutated.append(1)
									elif jspot == 1 and kspot == 1:
										tempindfather.append(kspot-1)
										# Count a mutation
										noallelesmutated.append(1)
									elif jspot != 0 and jspot != 1:
										# Randomly choose another allele								
										randallelespot = int(len(fathergenes[jspot])*rand())
										while randallelespot == kspot:
											randallelespot = int(len(fathergenes[jspot])*rand())
										tempindfather.append(randallelespot)									
										# Count a mutation
										noallelesmutated.append(1)
									else:
										# THen no mutation
										tempindfather.append(kspot)	
								
								# No other mutation models matched
								else:
									print('The mutation model does not exist.')
									sys.exit(-1)								
							
							# and if random number is not less than or equal to muterate
							else:
								tempindfather.append(kspot)
						
						# Check the homogeneous 2 case - in mother genes	
						if int(mothergenes[jspot][kspot])==2:
							
							# Get a random number for allele mutation
							mutationrandno = rand()
							
							# Check if random number is less than or equal to muterate
							if mutationrandno <= muterate:
							
								# If backward and forward mutation
								if mutationans == 'random':
									# Randomly choose another allele
									randallelespot = int(len(mothergenes[jspot])*rand())
									while randallelespot == kspot:
										randallelespot = int(len(mothergenes[jspot])*rand())
									tempindmother.append(randallelespot)
									# Count a mutation
									noallelesmutated.append(1)
									
								# If just forward mutation
								elif mutationans == 'forward':
									if kspot != len(mothergenes[jspot])-1:
										tempindmother.append(kspot+1)									
										# Count a mutation
										noallelesmutated.append(1)
									else:
										# THen no mutation
										tempindmother.append(kspot)

								# If just backward mutation
								elif mutationans == 'backward':
									if kspot != 0:
										tempindmother.append(kspot-1)									
										# Count a mutation
										noallelesmutated.append(1)
									else:
										# THen no mutation
										tempindmother.append(kspot)
										
								# If forward or backward step mutation
								elif mutationans == 'forwardbackward':
									# Random draw for forward or backward
									randstep = rand()
									# TO go left, but it can't be the first allele
									if randstep < 0.5 and kspot != 0:
										tempindmother.append(kspot-1)									
										# Count a mutation
										noallelesmutated.append(1)
									# To go right, but it can't be the last allele
									elif randstep >= 0.5 and kspot != len(fathergenes[jspot])-1:
										tempindmother.append(kspot+1)									
										# Count a mutation
										noallelesmutated.append(1)
									else:
										# THen no mutation
										tempindmother.append(kspot)
								
								# If forward mutation in A and backward mutation in b
								elif mutationans == 'forwardAbackwardBrandomN':
									if jspot == 0 and kspot == 0:
										tempindmother.append(kspot+1)									
										# Count a mutation
										noallelesmutated.append(1)
									elif jspot == 1 and kspot == 1:
										tempindmother.append(kspot-1)
										# Count a mutation
										noallelesmutated.append(1)
									if jspot != 0 and jspot != 1:
										# Randomly choose another allele
										randallelespot = int(len(mothergenes[jspot])*rand())
										while randallelespot == kspot:
											randallelespot = int(len(mothergenes[jspot])*rand())
										tempindmother.append(randallelespot)
										# Count a mutation
										noallelesmutated.append(1)
									else:
										# THen no mutation
										tempindmother.append(kspot)
								
								# No other mutation models matched
								else:
									print('The mutation model does not exist.')
									sys.exit(-1)					

							# and if random number is not less than or equal to muterate
							else:
								tempindmother.append(kspot)
							
						# Check the hetero 1 case for father genes
						if int(fathergenes[jspot][kspot])==1:
							tempindfather.append(kspot)
							
						# Check the hetero 1 case for mother genes
						if int(mothergenes[jspot][kspot])==1:
							tempindmother.append(kspot)
					
					# Check if the tempindex has a length of 2 (which means it was not homo at
					#	at this locus), then randomly select one of them, and then check for mutation
					# Check from father genes
					if len(tempindfather) == 2:
					
						# Then randomly select one of the homo alleles
						temprandnofather = int(2*rand())
						
						# Delete from list
						del(tempindfather[temprandnofather])
						thealleleselected = tempindfather[0]
						
						# Get a random number for allele mutation
						mutationrandno = rand()
						
						# Check if random number is less than or equal to muterate
						if mutationrandno <= muterate:
							
							# If backward and forward mutation
							if mutationans == 'random':
								# Randomly choose another allele
								randallelespot = int(len(fathergenes[jspot])*rand())
								while randallelespot == thealleleselected:
									randallelespot = int(len(fathergenes[jspot])*rand())
								# and then reassign this spot
								tempindfather[0] = randallelespot							
								# Count a mutation
								noallelesmutated.append(1)
								
							# Else if just forward mutation
							elif mutationans == 'forward':
								if thealleleselected != len(fathergenes[jspot])-1:
									# and then reassign this spot
									tempindfather[0] = thealleleselected+1									
									# Count a mutation
									noallelesmutated.append(1)
																
							# Else if just backward mutation
							elif mutationans == 'backward':
								if thealleleselected != 0:
									# and then reassign this spot
									tempindfather[0] = thealleleselected-1									
									# Count a mutation
									noallelesmutated.append(1)
														
							# Else if forward backward mutation
							elif mutationans == 'forwardbackward':
								# Random draw for forward or backward
								randstep = rand()
								# TO go left, but it can't be the first allele
								if randstep < 0.5 and thealleleselected != 0:
									# and then reassign this spot
									tempindfather[0] = thealleleselected-1									
									# Count a mutation
									noallelesmutated.append(1)
								# To go right, but it can't be the last allele
								elif randstep >= 0.5 and thealleleselected != len(fathergenes[jspot])-1:
									# and then reassign this spot
									tempindfather[0] = thealleleselected+1									
									# Count a mutation
									noallelesmutated.append(1)
							
							# Else if forward mutation for A and backward mutation for b
							elif mutationans == 'forwardAbackwardBrandomN':
								if jspot == 0 and thealleleselected == 0:
									# and then reassign this spot
									tempindfather[0] = thealleleselected+1									
									# Count a mutation
									noallelesmutated.append(1)									
								elif jspot == 1 and thealleleselected == 1:
									tempindfather.append(kspot-1)
									# Count a mutation
									noallelesmutated.append(1)
								elif jspot != 0 and jspot != 1:
									# Randomly choose another allele
									randallelespot = int(len(fathergenes[jspot])*rand())
									while randallelespot == thealleleselected:
										randallelespot = int(len(fathergenes[jspot])*rand())
									# and then reassign this spot
									tempindfather[0] = randallelespot							
									# Count a mutation
									noallelesmutated.append(1)
								else:
									# THen no mutation
									tempindfather[0] = thealleleselected	
										
							# No other mutation models matched
							else:
								print('The mutation model does not exist.')
								sys.exit(-1)		
					
					# Check from mother genes
					if len(tempindmother) == 2:
					
						# THen randomly select on of the homo alleles
						temprandnomother = int(2*rand())
						# Delete from list
						del(tempindmother[temprandnomother])
						thealleleselected = tempindmother[0]
						
						# Get a random number for allele mutation
						mutationrandno = rand()
						
						# Check if random number is less than or equal to muterate
						if mutationrandno <= muterate:
						
							# If backward and forward mutation
							if mutationans == 'random':
								# Randomly choose another allele
								randallelespot = int(len(mothergenes[jspot])*rand())
								while randallelespot == thealleleselected:
									randallelespot = int(len(mothergenes[jspot])*rand())
								# and then reassign this spot
								tempindmother[0] = randallelespot
								# Count a mutation
								noallelesmutated.append(1)
								
							# Else if just forward mutation
							elif mutationans == 'forward':
								if thealleleselected != len(mothergenes[jspot])-1:
									# and then reassign this spot
									tempindmother[0] = thealleleselected+1									
									# Count a mutation
									noallelesmutated.append(1)
																
							# Else if just backward mutation
							elif mutationans == 'backward':
								if thealleleselected != 0:
									# and then reassign this spot
									tempindmother[0] = thealleleselected-1									
									# Count a mutation
									noallelesmutated.append(1)
																
							# Else if forward backward mutation
							elif mutationans == 'forwardbackward':
								# Random draw for forward or backward
								randstep = rand()
								# TO go left, but it can't be the first allele
								if randstep < 0.5 and thealleleselected != 0:
									# and then reassign this spot
									tempindmother[0] = thealleleselected-1									
									# Count a mutation
									noallelesmutated.append(1)
								# To go right, but it can't be the last allele
								elif randstep >= 0.5 and thealleleselected != len(mothergenes[jspot])-1:
									# and then reassign this spot
									tempindmother[0] = thealleleselected+1									
									# Count a mutation
									noallelesmutated.append(1)
							
							# Else if forward mutation for A and backward mutation for b
							elif mutationans == 'forwardAbackwardBrandomN':
								if jspot == 0 and thealleleselected == 0:
									# and then reassign this spot
									tempindmother[0] = thealleleselected+1								
									# Count a mutation
									noallelesmutated.append(1)								
								elif jspot == 1 and thealleleselected == 1:
									tempindmother.append(kspot-1)
									# Count a mutation
									noallelesmutated.append(1)
								elif jspot != 0 and jspot != 1:
									# Randomly choose another allele
									randallelespot = int(len(mothergenes[jspot])*rand())
									while randallelespot == thealleleselected:
										randallelespot = int(len(mothergenes[jspot])*rand())
									# and then reassign this spot
									tempindmother[0] = randallelespot
									# Count a mutation
									noallelesmutated.append(1)
								else:
									# THen no mutation
									tempindmother[0] = thealleleselected
							
							# No other mutation models matched
							else:
								print('The mutation model does not exist.')
								sys.exit(-1)		
					
					# Now write to offspring genes array the selected alleles in locus j
					for kspot in xrange(len(fathergenes[jspot])):
						
						# Hetero case 1 AB
						if tempindfather[0] == kspot and tempindmother[0] != kspot:
							offspring[i][6].append(1)
						# Homo case AA or BB
						elif tempindfather[0] == kspot and tempindmother[0] == kspot:
							offspring[i][6].append(2)
						# Hetero case 2 BA
						elif tempindmother[0] == kspot and tempindfather[0] != kspot:
							offspring[i][6].append(1)
						# Or nothing there at all
						elif tempindmother[0] != kspot and tempindfather[0] != kspot:
							offspring[i][6].append(0)
				
				# If mtdna is turned on, then erase the last loci and force it to be mothergenes
				if mtdna == 'Y':
					
					# Force last locus to be mothergenes
					for imtdna in xrange(len(mothergenes[loci-1])):
						offspring[i][6][len(offspring[i][6])-(len(mothergenes[loci-1]))+imtdna] = mothergenes[loci-1][imtdna]
					
			# Now store the total number of alleles that mutated
			AllelesMutated.append(sum(noallelesmutated))
			
			# Delete temp variables to free up space
			del(mothergenes)
			del(fathergenes)
			del(tempindmother)
			del(tempindfather)
			
		# If there are no offspring
		elif int(offspringno) == int(0):
		
			# Store the total number of alleles that mutated
			AllelesMutated.append(0.0)				
	
	# If not generation to swap genes
	else:
		
		# Store the total number of alleles that mutated
		AllelesMutated.append(0.0)	
		
	# Return variables from this argument
	return offspring
	
	# End::InheritGenes()

# ---------------------------------------------------------------------------------------------------
def ConstantMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Deaths,alleles,FID,agemort,infection,geneswap,mature):
	
	'''
	Constant mortality applied to each age class
	'''	
	
	# Get total number of deaths for each age class, be careful of NAs and strings
	uniqueages = Counter(np.asarray(np.asarray(age,dtype='|S10')[np.where(np.asarray(age,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	agedeaths = []	
	extra_agedeaths = []
	for i in xrange(1,(len(agemort)+1)):
		agedeaths.append(round(agemort[i-1]*uniqueages[i]))
	
	# Switch for ages over age classes, apply last mortality to age classes
	if max(uniqueages) > len(agemort):
		print('Warning: age classes exceeding specified class in Agevars.csv file.')
		for j in xrange(len(agemort)+1,(max(uniqueages)+1)):
			extra_agedeaths.append(round(agemort[-1]*uniqueages[j]))
	
	# Grab locations that are open
	openindex = np.where(np.asarray(sex) == 'NA')[0]
	
	# Grab locations that are not open
	filledindex = np.where(np.asarray(sex) != 'NA')[0]
	
	# Then take a sample from the possible age class indices to delete from
	deleteoldindex = []	
	for i in xrange(1,(len(agedeaths)+1)):			
		deleteoldindex.append(random.sample(np.where(np.asarray(age,dtype = 'str') == str(i))[0],int(agedeaths[i-1])))
		
	# In case there are extra age deaths
	if len(extra_agedeaths) != 0:
		count = len(agemort)+1
		for j in xrange(len(extra_agedeaths)):
			deleteoldindex.append(random.sample(np.where(np.asarray(age) == count)[0],int(extra_agedeaths[j])))
			count = count + 1
	
	# Flatten and turn into array
	deleteoldindex = np.asarray([item for sublist in deleteoldindex for item in sublist])
		
	# Then add indices together
	deleteallindex = np.append(openindex,deleteoldindex)
	deleteallindex = np.asarray(deleteallindex,dtype = 'int')
	
	# Store freegrid locations
	freegrid = np.asarray(FID)[deleteallindex]
	freegrid = list(freegrid)	
	# Delete all of the old generations information.
	sex = np.delete(sex,deleteallindex)
	id = np.delete(id,deleteallindex)
	age = np.delete(age,deleteallindex)
	xgrid = np.delete(xgrid,deleteallindex)
	ygrid = np.delete(ygrid,deleteallindex)
	mature = np.delete(mature,deleteallindex)
	# Keep genes around if within burnin faze
	if gen >= geneswap:
		genes = np.delete(genes,deleteallindex,axis=0)
	FID = np.delete(FID,deleteallindex)
	infection = np.delete(np.asarray(infection),deleteallindex)
	
	# Just one more shuffle to mix up empty versus killed off grids
	shuffle(freegrid)
	
	# Store total number of Deaths information for output 
	# Check that agedeaths equals age distribution file
	Deaths.append(agedeaths)
	if len(extra_agedeaths) != 0:
		Deaths[gen][-1] = Deaths[gen][-1] + sum(extra_agedeaths)
	
	# Return variables from this argument
	tupMort = freegrid,id,sex,age,xgrid,ygrid,genes,FID,infection,mature
	return tupMort
	# End::DoConstantMortality
	
# ----------------------------------------------------------------------------------------------	 
def DDMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Deaths,alleles,FID,infection,geneswap,K_env,popmodel,agemort,mature):
	'''
	DensityDependentMortality()
	Density dependent survival applied to each population.		
	'''
	
	# Get total number of deaths for each age class, be careful of NAs and strings
	uniqueages = Counter(np.asarray(np.asarray(age,dtype='|S10')[np.where(np.asarray(age,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	Nt = len(np.where(np.asarray(age,dtype='|S10') != 'NA')[0])
	agedeaths = []	
	extra_agedeaths = []
	
	for i in xrange(1,(len(agemort)+1)): # Only for ages 1 +		
		# Current Ni,t
		Ni = uniqueages[i]		
		# If last age class, no survivors
		if i == len(agemort):			
			mortreturn = Ni
		else:
			# Next age class
			Niplus1 = uniqueages[i+1]
			
			# Rickers option
			if popmodel == 'Rickers' or popmodel == 'rickers':
				print('Not operational.')
				sys.exit(-1)
			# Richards
			elif popmodel == 'Richards' or popmodel == 'richards':
				print('Not operational.')
				sys.exit(-1)
			# Logistic
			elif popmodel == 'logistic' or popmodel == 'Logistic':
				# Ni+1,t * (1 - Nt/K) * (si-1 * Ni,t - Ni+1,t)
				Ntplus1 = Niplus1 + (1. - Nt/float(K_env))*((1.-agemort[i-1]) * Ni - float(Niplus1))
				mortreturn = Ni - int(Ntplus1)
				if mortreturn < 0:
					mortreturn = 0
			else:
				print('Enter rickers, richards, or logistic for density dependent option')
				sys.exit(-1)
		
		agedeaths.append(mortreturn)
	
	# Switch for ages over age classes, apply last mortality to age classes
	if max(uniqueages) > len(agemort):
		print('Warning: age classes exceeding specified class in Agevars.csv file, apply 0 survival.')
		for j in xrange(len(agemort)+1,(max(uniqueages)+1)):
			extra_agedeaths.append(uniqueages[j])
			
	# Grab locations that are open
	openindex = np.where(np.asarray(sex) == 'NA')[0]
	
	# Grab locations that are not open
	filledindex = np.where(np.asarray(sex) != 'NA')[0]
	
	# Then take a sample from the possible age class indices to delete from
	deleteoldindex = []
	for i in xrange(1,(len(agedeaths)+1)): # index from 1 to ...
		# NA switch
		if len(openindex) == 0:
			deleteoldindex.append(random.sample(np.where(np.asarray(age) == i)[0],int(agedeaths[i-1]))) # index into age deaths - 1
		else:
			deleteoldindex.append(random.sample(np.where(np.asarray(age) == str(i))[0],int(agedeaths[i-1])))	
	# In case there are extra age deaths
	if len(extra_agedeaths) != 0:
		count = len(agemort)+1
		for j in xrange(len(extra_agedeaths)):
			deleteoldindex.append(random.sample(np.where(np.asarray(age) == count)[0],int(extra_agedeaths[j])))
			count = count + 1
	
	# Flatten and turn into array
	deleteoldindex = np.asarray([item for sublist in deleteoldindex for item in sublist],dtype = 'int')
		
	# Then add indices together
	deleteallindex = np.append(openindex,deleteoldindex)
	deleteallindex = np.asarray(deleteallindex,dtype='int')
	
	# Store freegrid locations
	freegrid = np.asarray(FID)[deleteallindex]
	freegrid = list(freegrid)	
	# Delete all of the old generations information.
	sex = np.delete(sex,deleteallindex)
	id = np.delete(id,deleteallindex)
	age = np.delete(age,deleteallindex)
	xgrid = np.delete(xgrid,deleteallindex)
	ygrid = np.delete(ygrid,deleteallindex)
	mature = np.delete(mature,deleteallindex)
	# Keep genes around if within burnin faze
	if gen >= geneswap:
		genes = np.delete(genes,deleteallindex,axis=0)
	FID = np.delete(FID,deleteallindex)
	infection = np.delete(np.asarray(infection),deleteallindex)
	
	# Just one more shuffle to mix up empty versus killed off grids
	shuffle(freegrid)
	
	# Store total number of Deaths information for output 
	# Check that agedeaths equals age distribution file
	Deaths.append(agedeaths)
	if len(extra_agedeaths) != 0:
		Deaths[gen][-1] = Deaths[gen][-1] + sum(extra_agedeaths)
	
	# Return variables from this argument
	tupMort = freegrid,id,sex,age,xgrid,ygrid,genes,FID,infection,mature
	return tupMort
	# End::DDMortality
	
# ---------------------------------------------------------------------------------------------------	 
def AdultSelection(tupMort,fitvals,mature,SelectionDeaths):
	'''
	AdultSelection()
	Mortality of old generation
	using selection values. Only for mature individuals
	'''
	# Unpack tuple
	freegrid = tupMort[0]
	id = tupMort[1]
	sex = tupMort[2]
	age = tupMort[3]
	xgrid = tupMort[4]
	ygrid = tupMort[5]
	genes = tupMort[6]	
	FID = tupMort[7]
	infection = tupMort[8]
	mature = tupMort[9]
	
	deleteallindex = []
	# Loop through mature individuals
	for i in xrange(len(mature)):
		if mature[i] == 1:
			# Find it's location
			usefitvals = fitvals[FID[i]]
			# Check genotype and match fitvals
			if genes[i][0][0] == 2: # AA
				diffmort = float(usefitvals[0])/100.
			elif genes[i][0][0] == 1 and genes[i][0][1] == 1: # Aa
				diffmort = float(usefitvals[1])/100.
			elif genes[i][0][1] == 2: #aa
				diffmort = float(usefitvals[2])/100.
			else: # Another genotype
				diffmort = 0.0
			
			# Then flip the coin to see if  individual survives its location
			randcheck = rand()			
			# If individual did not survive: keep track of delete ones
			if randcheck < diffmort:
				SelectionDeaths.append(1) # Record
				# Keep storage spots to delete
				deleteallindex.append(i)				
				# Add to freegrid
				freegrid.append(FID[i])
	
	# Delete spots
	deleteallindex = np.asarray(deleteallindex,dtype='int')			
	sex = np.delete(sex,deleteallindex)
	id = np.delete(id,deleteallindex)
	age = np.delete(age,deleteallindex)
	xgrid = np.delete(xgrid,deleteallindex)
	ygrid = np.delete(ygrid,deleteallindex)
	mature = np.delete(mature,deleteallindex)
	genes = np.delete(genes,deleteallindex,axis=0)
	FID = np.delete(FID,deleteallindex)
	infection = np.delete(np.asarray(infection),deleteallindex)			
	
	# Return variables from this argument
	tupMort = freegrid,id,sex,age,xgrid,ygrid,genes,FID,infection
	return tupMort
	# End::AdultSelection
	
# ---------------------------------------------------------------------------------------------------	 
def DoMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Deaths,alleles,FID,agemort,infection,geneswap,popmodel,K_env,fitvals,mature,cdevolveans,Opt3SelectionDeaths,burningen):
	'''
	DoMortality()
	Mortality of old generation
	Input: Adult mortality% 
	Output: Old files minus the killed off individuals:
	freegrid = [xycdmatid location of free grid spot 
	in random order]
	Constant or Density dependent functions here.
	'''
	
	Opt3SelectionDeaths.append([]) # Spot added for generation
	# Switch for model choice
	if popmodel == 'exp':
	
		tupMort = ConstantMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Deaths,alleles,FID,agemort,infection,geneswap,mature)
		
	elif popmodel == 'logistic' or popmodel == 'rickers' or popmodel == 'richards':
		
		tupMort = DDMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Deaths,alleles,FID,infection,geneswap,K_env,popmodel,agemort,mature)	
		
	else:
		print('This population model for survival does not exist')
		sys.exit(-1)
	
	# Then apply fitness values to adults if there and mature
	if cdevolveans == '3' and gen >= burningen:
		tupMort = AdultSelection(tupMort,fitvals,mature,Opt3SelectionDeaths[gen])
	Opt3SelectionDeaths[gen] = sum(Opt3SelectionDeaths[gen]) 
		
	return tupMort
	
	# End::DoAdultMortality()