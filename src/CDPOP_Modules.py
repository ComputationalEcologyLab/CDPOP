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
Alleles,He,subpop,p1,p2,q1,q2,Population_age,age,Magemort,geneswap,\
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
	for i in xrange(len(Magemort)):
		Population_age[gen].append([])
	countage = Counter(np.asarray(np.asarray(age,dtype='|S10')[np.where(np.asarray(age,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	for i in xrange(1,len(Magemort)+1):
		Population_age[gen][i-1].append(countage[i])
	
	# Sum up population age tracker
	for i in xrange(len(Magemort)):
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
	Output: [femaleid,maleid,cdmatidofmother,cdmatidoffather,sex,infection,TWINID],[
	genetic information]		
	'''		
	
	# Check for generation to start swapping genes
	if gen >= geneswap:
	
		# Storage for tracking how many alleles mutated
		noallelesmutated = []			
		# If there are offspring
		if int(offspringno) != int(0):
		
			# Sort by twins. 
			offspring.sort(key=lambda x: x[6]) 
			tempoff = np.asarray(copy.deepcopy(offspring))
			countTwins = Counter(tempoff[:,6])
	
			isTwin = False # For checking for twins
			twingenes = [] # Initialize twin genes
			# Begin loop through offspring
			for i in xrange(offspringno):				
				
				# If twins genes were copied already from previous offspring
				if len(twingenes) == 0: 
					# Temp storage for i's mother's  and fathers genes
					mothergenes=genes[int(offspring[i][0])]
					fathergenes=genes[int(offspring[i][1])]								
					fathergenes = sum(fathergenes,[])
					mothergenes = sum(mothergenes,[])
					fathergenes = np.asarray(fathergenes,dtype=int)
					mothergenes = np.asarray(mothergenes,dtype=int)
					# Temp genes storage for offspring
					tempgenes = np.zeros(len(fathergenes),dtype =int)
					# Allele indices
					alleles = np.asarray(range(len(mothergenes))) 
					for iloci in xrange(loci): # Loop through loci
						# Allele indices to sample from - index into tempgenes
						possiblealleles = alleles[(iloci*len(mothergenes)/loci):(iloci*len(mothergenes)/loci+len(mothergenes)/loci)]
						
						# Father and mother locations							
						F2 = np.where(fathergenes[possiblealleles] == 2)[0] # location of 2s
						F1 = np.where(fathergenes[possiblealleles] == 1)[0]
						M2 = np.where(mothergenes[possiblealleles] == 2)[0]
						M1 = np.where(mothergenes[possiblealleles] == 1)[0]
						tempALL = np.concatenate((F2,F2,F1,M2,M2,M1),axis=0) # Repeat 2s twice
						if len(tempALL) < 4:
							print('Error in twinning. Email Erin.')
							sys.exit(-1)
							
						# Sample 2 alleles - these are indeices
						sampleAlleles = random.sample(tempALL,2)
						# Fill in alleles corresponding to sampled spots
						for iall in xrange(2):
							tempgenes[possiblealleles[sampleAlleles[iall]]] = tempgenes[possiblealleles[sampleAlleles[iall]]] + 1
							
					# mtDNA is turned on
					# ------------------
					if mtdna == 'Y':
						# Force last locus to be mothergenes - possible alleles are from the last loop above
						tempgenes[possiblealleles] = mothergenes[possiblealleles]				
						
					# Check for twin copy
					# -------------------
					if len(offspring[i][6].split('T')) > 1:
						# Check if there are indeed more than 1 of these 
						if countTwins[offspring[i][6]] > 1:
							twingenes = copy.deepcopy(tempgenes)
										
				# Then flip the twin check back for next offsprings
				else: # This is a twin
					tempgenes = copy.deepcopy(twingenes)
					isTwin = True
					# Check
					if offspring[i][6] != offspring[i-1][6]:
						print('Twinning algorithm is not working. Email Erin.')
						sys.exit(-1)
														
				# Then check for mutations at each allele
				# ---------------------------------------
				if muterate != 0.0:
					for iloci in xrange(loci): # Loop through loci
						mutationrandnos = rand(2) # Get a random number for checking
						# Allele indices to sample from - index into tempgenes
						possiblealleles = alleles[(iloci*len(mothergenes)/loci):(iloci*len(mothergenes)/loci+len(mothergenes)/loci)]
						# Get the current location of alleles - index into tempgenes 
						thisloci = possiblealleles[np.where(tempgenes[possiblealleles] != 0)[0]]
						# Check case for homo
						if len(thisloci) == 1:
							# Copy the spot
							thisloci = np.concatenate((thisloci,thisloci),axis=0)
						
						# Loop through alleles
						for iall in xrange(2): 			
							
							# Check if random number is less than muterate
							if mutationrandnos[iall] < muterate:
								
								# First remove this allele from tempgenes
								tempgenes[thisloci[iall]] = tempgenes[thisloci[iall]] - 1
																
								# If random kth allele model
								if mutationans == 'random':
									# Randomly choose another allele, but not what allele it was									
									movealleleTO = random.sample(possiblealleles[np.where(thisloci[iall] != possiblealleles)[0]],1)[0]
									# Index into tempgenes and add 1
									tempgenes[movealleleTO] = tempgenes[movealleleTO] + 1
																	
									# Count a mutation
									noallelesmutated.append(1)
									
								# If just forward mutation
								elif mutationans == 'forward':
									# Move allele forward unless it is the last one
									if thisloci[iall] != possiblealleles[-1]:
										tempgenes[thisloci[iall]+1] = tempgenes[thisloci[iall]+1] + 1
																		
										# Count a mutation
										noallelesmutated.append(1)

								# If just forward mutation
								elif mutationans == 'backward':
									# Move allele backward unless it is the first one
									if thisloci[iall] != possiblealleles[0]:
										tempgenes[thisloci[iall]-1] = tempgenes[thisloci[iall]-1] + 1
																		
										# Count a mutation
										noallelesmutated.append(1)
										
								# If forward and backward mutation
								elif mutationans == 'forwardbackward':
									# Then random forward or backward step
									randstep = rand()
									# To go left, but it can't be the first allele
									if randstep < 0.5 and thisloci[iall] != possiblealleles[0]:
										tempgenes[thisloci[iall]-1] = tempgenes[thisloci[iall]-1] + 1
										# Count a mutation
										noallelesmutated.append(1)
									# To go right, but it can't be the last allele
									elif randstep >= 0.5 and thisloci[iall] != possiblealleles[-1]:
										tempgenes[thisloci[iall]+1] = tempgenes[thisloci[iall]+1] + 1
										# Count a mutation
										noallelesmutated.append(1)
								
								# If forward mutation in A and backward mutation for b (A -> a, b -> B)
								elif mutationans == 'forwardAbackwardBrandomN':
									print('Currently not operating. Email Erin.')
									sys.exit(-1)
									if iloci == 0 and thisloci[iall] == possiblealleles[0]:
										tempgenes[thisloci[iall]+1] = tempgenes[thisloci[iall]+1] + 1									
										# Count a mutation
										noallelesmutated.append(1)
									elif iloci == 1 and thisloci[iall] == possiblealleles[1]:
										tempgenes[thisloci[iall]-1] = tempgenes[thisloci[iall]-1] + 1
										# Count a mutation
										noallelesmutated.append(1)
									elif iloci != 0 and iloci != 1:
										# Randomly choose another allele								
										movealleleTO = random.sample(possiblealleles[np.where(thisloci[iall] != possiblealleles)[0]],1)[0]
										# Index into tempgenes and add 1
										tempgenes[movealleleTO] = tempgenes[movealleleTO] + 1
																		
										# Count a mutation
										noallelesmutated.append(1)
								
								# No other mutation models matched
								else:
									print('The mutation model does not exist.')
									sys.exit(-1)	
				
				# Add to offspring list
				# ---------------------
				# Add spot in offspring array for individual i's genes
				offspring[i].append(tempgenes.tolist())
				
				# Copy the saved twingenes for the next offspring
				if isTwin:
					twingenes = [] # Erase the copied twin genes
					isTwin = False # Turn off 
					
			# Now store the total number of alleles that mutated
			AllelesMutated.append(sum(noallelesmutated))
			
			# Delete temp variables to free up space
			del(mothergenes)
			del(fathergenes)
			#del(tempgenes)
			
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
def ConstantMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,alleles,FID,Magemort,Fagemort,infection,geneswap,mature):
	
	'''
	Constant mortality applied to each age class
	'''	
	
	# Get total number of deaths for each age class, be careful of NAs and strings
	uniqueages = Counter(np.asarray(np.asarray(age,dtype='|S10')[np.where(np.asarray(age,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	# Split up for sex
	females = np.where(np.asarray(sex,dtype='|S4')=='0')[0]
	males = np.where(np.asarray(sex,dtype='|S4')=='1')[0]
	Fages = np.asarray(age)[females]
	Mages = np.asarray(age)[males]
	Funiqueages = Counter(np.asarray(np.asarray(Fages,dtype='|S10')[np.where(np.asarray(Fages,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	Muniqueages = Counter(np.asarray(np.asarray(Mages,dtype='|S10')[np.where(np.asarray(Mages,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	Magedeaths = []	
	Fagedeaths = []
	extra_Magedeaths = []
	extra_Fagedeaths = []
	for i in xrange(1,(len(Magemort)+1)):
		Magedeaths.append(round(Magemort[i-1]*Muniqueages[i]))
		Fagedeaths.append(round(Fagemort[i-1]*Funiqueages[i]))
	
	# Switch for ages over age classes, apply last mortality to age classes
	if max(Muniqueages) > len(Magemort):
		print('Warning: age classes exceeding specified class in Agevars.csv file.')
		for j in xrange(len(Magemort)+1,(max(Muniqueages)+1)):
			extra_Magedeaths.append(round(Magemort[-1]*Muniqueages[j]))
	if max(Funiqueages) > len(Fagemort):
		print('Warning: age classes exceeding specified class in Agevars.csv file.')
		for j in xrange(len(Fagemort)+1,(max(Funiqueages)+1)):
			extra_Fagedeaths.append(round(Fagemort[-1]*Funiqueages[j]))
	
	# Grab locations that are open
	openindex = np.where(np.asarray(sex) == 'NA')[0]
	
	# Grab locations that are filled
	filledindex = np.where(np.asarray(sex) != 'NA')[0]
	
	# Then take a sample from the possible age class indices to delete from
	Mdeleteoldindex = []	
	Fdeleteoldindex = []
	for i in xrange(1,(len(Magedeaths)+1)):			
		Mdeleteoldindex.append(random.sample(np.where(np.asarray(Mages,dtype = 'str') == str(i))[0],int(Magedeaths[i-1])))
		Fdeleteoldindex.append(random.sample(np.where(np.asarray(Fages,dtype = 'str') == str(i))[0],int(Fagedeaths[i-1])))
		
	# In case there are extra age deaths
	if len(extra_Magedeaths) != 0:
		count = len(Magemort)+1
		for j in xrange(len(extra_Magedeaths)):
			Mdeleteoldindex.append(random.sample(np.where(np.asarray(Mages) == count)[0],int(extra_Magedeaths[j])))
			count = count + 1
	if len(extra_Fagedeaths) != 0:
		count = len(Fagemort)+1
		for j in xrange(len(extra_Fagedeaths)):
			Fdeleteoldindex.append(random.sample(np.where(np.asarray(Fages) == count)[0],int(extra_Fagedeaths[j])))
			count = count + 1
	
	# Flatten and turn into array and get original index locations and add indices together
	if sum(Mdeleteoldindex,[]) != []:
		Mdeleteoldindex = males[np.asarray(sum(Mdeleteoldindex,[]))]
		deleteallindex = np.append(openindex,Mdeleteoldindex)
	else:
		Mdeleteoldindex = np.asarray([])
		deleteallindex = openindex
	if sum(Fdeleteoldindex,[]) != []:
		Fdeleteoldindex = females[np.asarray(sum(Fdeleteoldindex,[]))]
		deleteallindex = np.append(deleteallindex,Fdeleteoldindex)
	else:
		Fdeleteoldindex = np.asarray([])
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
		genes = np.asarray(genes,dtype=int)
		genes = genes.tolist()
	FID = np.delete(FID,deleteallindex)
	infection = np.delete(np.asarray(infection),deleteallindex)
	
	# Just one more shuffle to mix up empty versus killed off grids
	shuffle(freegrid)
	
	# Store total number of Deaths information for output	
	Track_MDeaths.append(Magedeaths)
	Track_FDeaths.append(Fagedeaths)
	# Check that agedeaths equals age distribution file
	if len(extra_Magedeaths) != 0:
		Track_MDeaths[gen][-1] = Track_MDeaths[gen][-1] + sum(extra_Magedeaths)
	if len(extra_Fagedeaths) != 0:
		Track_FDeaths[gen][-1] = Track_FDeaths[gen][-1] + sum(extra_Fagedeaths)
	
	# Return variables from this argument
	tupMort = freegrid,id,sex,age,xgrid,ygrid,genes,FID,infection,mature
	return tupMort
	# End::DoConstantMortality
	
# ----------------------------------------------------------------------------------------------	 
def DDMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,alleles,FID,infection,geneswap,K_env,popmodel,Magemort,Fagemort,mature):
	'''
	DensityDependentMortality()
	Density dependent survival applied to each population.		
	'''
	
	# Get total number of deaths for each age class, be careful of NAs and strings
	uniqueages = Counter(np.asarray(np.asarray(age,dtype='|S10')[np.where(np.asarray(age,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	Nt = len(np.where(np.asarray(age,dtype='|S10') != 'NA')[0])
	agedeaths = []	
	extra_agedeaths = []
	
	if Magemort != Fagemort:
		print('Warning: Logistic growth is specified and the average age specific mortality of males and females will be used.')
	agemort = (np.asarray(Magemort) + np.asarray(Fagemort)) / 2.
	
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
		genes = np.asarray(genes,dtype=int)
		genes = genes.tolist()
	FID = np.delete(FID,deleteallindex)
	infection = np.delete(np.asarray(infection),deleteallindex)
	
	# Just one more shuffle to mix up empty versus killed off grids
	shuffle(freegrid)

	# Store total number of Deaths information for output 
	# Check that agedeaths equals age distribution file
	Track_MDeaths.append((np.asarray(agedeaths)/2.).tolist())
	Track_FDeaths.append((np.asarray(agedeaths)/2.).tolist())
	if len(extra_agedeaths) != 0:
		Track_MDeaths[gen][-1] = Track_MDeaths[gen][-1] + sum(extra_agedeaths)/2.
		Track_FDeaths[gen][-1] = Track_FDeaths[gen][-1] + sum(extra_agedeaths)/2.
	
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
def DoMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,alleles,FID,Magemort,Fagemort,infection,geneswap,popmodel,K_env,fitvals,mature,cdevolveans,Opt3SelectionDeaths,burningen):
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
	
		tupMort = ConstantMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,alleles,FID,Magemort,Fagemort,infection,geneswap,mature)
		
	elif popmodel == 'logistic' or popmodel == 'rickers' or popmodel == 'richards':
		
		tupMort = DDMortality(filledgrids,nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,alleles,FID,infection,geneswap,K_env,popmodel,Magemort,Fagemort,mature)	
		
	else:
		print('This population model for survival does not exist')
		sys.exit(-1)
	
	# Then apply fitness values to adults if there and mature
	if cdevolveans == '3' and gen >= burningen:
		tupMort = AdultSelection(tupMort,fitvals,mature,Opt3SelectionDeaths[gen])
	Opt3SelectionDeaths[gen] = sum(Opt3SelectionDeaths[gen]) 
		
	return tupMort
	
	# End::DoAdultMortality()