# -------------------------------------------------------------------------------------------------
# CDPOP_Offspring.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for offspring processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	from numpy.random import *
	import numpy as np
except ImportError:
	raise ImportError, "Numpy required."
import pdb,sys
from sets import Set
import random
from collections import Counter

# ---------------------------------------------------------------------------------------------------	
def DoOffspringSex(Bearpairs,Femalepercent,CDpairs,equalsexratio):
	'''
	DoOffspringSex()
	This function assigns the sex of each offspring.
	'''	
	# Create empty variable for storing individual offspring information
	offspring=[]
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
		
		# If equal sex ratio is N
		if equalsexratio == 'N':
		
			# And then loop through each offspring from that mate pair
			for j in xrange(Bearpairs[i][2]):
				
				# Select sex of the jth offspring - select a random number
				randsex = int(100*rand())
				
				# If that random number is less the Femalepercent, assign it to be a female
				if randsex < Femalepercent:
					offsex = 0
				
				# If the random number is greater than the Femalepercent, assign it to be a male
				else:
					offsex = 1

				# And then append all information onto a list storage variable offspring
				offspring.append([Bearpairs[i][0],Bearpairs[i][1],CDpairs[i][0],CDpairs[i][1],offsex])
		
		# If equal sex ratio is Y
		elif equalsexratio == 'AtBirth':
		
			# Error statement for if not even number
			if np.mod(Bearpairs[i][2],2) != 0:
				print('Equal sex ratio does not work for odd litter size. Use offno = 3 and even lambda value.')
				sys.exit(-1)
			
			offsex = np.append(np.zeros(Bearpairs[i][2]/2,"int"),np.ones(Bearpairs[i][2]/2,"int"))
			np.random.shuffle(offsex)
			
			# And then loop through each offspring from that mate pair
			for j in xrange(Bearpairs[i][2]):
			
				# And then append all information onto a list storage variable offspring
				offspring.append([Bearpairs[i][0],Bearpairs[i][1],CDpairs[i][0],CDpairs[i][1],offsex[j]])
	
	# Variables returned
	return offspring
	
	# End::DoOffspringSex()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringRandom(Bearpairs,CDpairs,lmbdavals,age):
	'''
	DoOffspringRandom()
	This function chooses a random number of 
	offspring for a mated pair.
	'''	
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
		
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
			CDpairs[i].append(0)
		
		# If females did mate up, then assign random drawn number
		else:
			
			# Get age of female
			Fage = int(age[Bearpairs[i][0]])
			if Fage >= len(lmbdavals)+1:
				Fage = len(lmbdavals)-1
			lmbda = lmbdavals[Fage-1]
			
			# Randomly choose a number between 0 and 4
			randkidno = int(int((lmbda))*rand())
			
			# Append Offspring number to end of Pairs [F,M,#offspring]
			Bearpairs[i].append(randkidno)
			CDpairs[i].append(randkidno)
	
	# Variables returned
	tupDoOffRand = Bearpairs,CDpairs
	return tupDoOffRand
	# End::DoOffspringRandom()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringPoisson(Bearpairs,CDpairs,lmbdavals,age):
	'''
	DoOffspringPoisson()
	This function chooses a number of offspring 
	from a Poisson distribution for a mated pair.
	'''		
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
			CDpairs[i].append(0)
		
		# If female did mate up, then assign a poisson draw with mean lmbda
		else:
			
			# Get age of female
			Fage = int(age[Bearpairs[i][0]])
			if Fage >= len(lmbdavals)+1:
				Fage = len(lmbdavals)-1
			lmbda = lmbdavals[Fage-1]
			
			# Poisson number with mean, lambda
			poissonkidno = poisson(lmbda)
			
			# Append Offspring number to end of Pairs [F,M,#offspring]
			Bearpairs[i].append(poissonkidno)
			CDpairs[i].append(poissonkidno)
			
	# Variables returned
	tupDoOffPois = Bearpairs,CDpairs
	return tupDoOffPois
	# End::DoOffspringPoisson()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringConstant(Bearpairs,CDpairs,lmbdavals,age):
	'''
	DoOffspringConstant()
	This function chooses a constant number of 
	offspring for each mated pair.
	'''	
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
			CDpairs[i].append(0)
		
		# If females did mate up, then assign lmbda constant
		else:
			
			# Get age of female
			Fage = int(age[Bearpairs[i][0]])
			if Fage >= len(lmbdavals)+1:
				Fage = len(lmbdavals)-1
			lmbda = lmbdavals[Fage-1]
			
			# Assign a 1 [F,M,#offspring]
			Bearpairs[i].append(int(lmbda))
			CDpairs[i].append(int(lmbda))
			
	# Variables returned
	tupDoOffConst = Bearpairs,CDpairs
	return tupDoOffConst
	# End::DoOffspringConstant()
	
# ---------------------------------------------------------------------------------------------------	
def DoOffspringEqual(offspring,lmbdavals,age):
	'''
	DoOffspringEqual()
	This function chooses lmbda number of offspring for each female.
	'''	
	# import up front not working...temporary.
	from sets import Set
	
	# If offspring list is not empty
	if len(offspring) != 0:		
				
		# Get array and list of females
		tempOff = np.asarray(offspring)
		tempFemales = list(tempOff[:,0])
		
		# Get unique Female locations
		uniqueSet = Set(item for item in tempFemales)

		# Loop through unique Set of females and create sample list
		returnOff = []
		for item in uniqueSet:
		
			# Get age of unique female
			Fage = int(age[item])
			if Fage >= len(lmbdavals)+1:
				Fage = len(lmbdavals)-1
			lmbda = lmbdavals[Fage-1] 
			# Get where this female is in offspring list
			tempFloc = np.where(tempFemales == item)[0]
			# Then random sample lmbda of them
			sampOffloc = random.sample(tempFloc,int(lmbda))
			for i in xrange(len(sampOffloc)):
				returnOff.append(offspring[sampOffloc[i]])
	else:
		returnOff = offspring
		
	return returnOff
	# End::DoOffspringEqual()

# ---------------------------------------------------------------------------------------------------
def DoOffspringNormal(Bearpairs,CDpairs,lmbdavals,sigmavals,age):
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
			CDpairs[i].append(0)
		
		# If female did mate up, then assign a normal draw
		else:
			
			# Get age of female
			Fage = int(age[Bearpairs[i][0]])
			if Fage >= len(lmbdavals)+1:
				Fage = len(lmbdavals)-1
			lmbda = lmbdavals[Fage-1]
			sigma = sigmavals[Fage-1]
			
			# Poisson number with mean, lambda
			normalkidno = int(np.random.normal(lmbda,sigma))
			if normalkidno < 0:
				normalkidno == 0
			
			# Append Offspring number to end of Pairs [F,M,#offspring]
			Bearpairs[i].append(normalkidno)
			CDpairs[i].append(normalkidno)
			
	# Variables returned
	tupDoOff = Bearpairs,CDpairs
	return tupDoOff
	# End::DoOffspringNormal()
	
# ---------------------------------------------------------------------------------------------------	 
def DoCDInfect(offspring,infection,transmissionprob):
	'''
	DoCDInfect()
	This function determins whether the offspring gets infected
	or not from the parents infection.  Vertical transmission.
	'''	
	# Get infection status of offspring.
	for ithinfect in xrange(len(offspring)):
	
		# If parent has infection
		if infection[offspring[ithinfect][0]] == 1 or\
		infection[offspring[ithinfect][1]] == 1:
		
			# Get a random number
			randinfection = rand()
			
			# If random flip is less than transmissionprob
			if randinfection < transmissionprob:
			
				# Then append infection status to offspring 
				offspring[ithinfect].append(1)
				
			# If offspring does not get infected
			else:
				offspring[ithinfect].append(0)
			
		# If offspring does not get infected.
		else:
		
			# Then append infection status to offspring
			offspring[ithinfect].append(0)

	# Return variables
	return offspring
		
	# End::DoCDInfect()	
	
# ----------------------------------------------------------------------------------------------	 
def DoOffMortality(offspring,Mnewmortperc,Fnewmortperc,equalsexratio,offno,age,lmbda,sex,Track_MOffDeaths,Track_FOffDeaths):
	'''
	DoOffMortality()
	Mortality functions for age 0 class.
	'''
	
	# s0 given - global survival for age0s
	s0_males = 1. - Mnewmortperc
	s0_females = 1. - Fnewmortperc
	
	# Special case for stable age distribution - adjust survival number
	if offno == '6':
	
		if (Mnewmortperc != Fnewmortperc):
			print('Warning: Age 0 mortality is different for males and females, and option 6 offno is being used. Using average survival.')
		s0 = 1. - ((Mnewmortperc + Fnewmortperc) / 2.)
		# Get unique ages
		countage = Counter(np.asarray(np.asarray(age,dtype='|S10')[np.where(np.asarray(age,dtype='|S10') != 'NA')[0]],dtype=np.int8))
		# Get female locations and age of each
		Fsex = np.where(np.asarray(sex) == '0')[0]
		Fsexage = np.asarray(age)[Fsex]
		Fsexage = np.asarray(Fsexage,dtype = int)		
		
		tempoffspring = [] # storage for survived offspring as loop through age classes
		# Loop through age classes - careful of indexing
		for i in xrange(len(lmbda)):
			# For indexing, get age here
			ageindex = i+1
			
			# Get total individuals in each age class
			Nt = countage[ageindex]
			# Get total females in each age class
			Ft = len(np.where(Fsexage == ageindex)[0])			
			# Get number of offspring produced by each age class
			Ft_offno = Ft * lmbda[i]
						
			if Ft == 0:
				break			
			
			# Get corrected survival number
			correctedS0 = s0*Nt/Ft
			
			# Find the location of this female age class - this is their unique ID/location
			Fspots = Fsex[np.where(Fsexage == ageindex)[0]]
		
			# Get mother locations in offspring list - this is all of the mother's ID that had offspring
			motheroff = np.asarray(offspring)[:,0]
			
			# For each Fspot, return offspring location index
			offspot = [] 
			for j in xrange(len(Fspots)):
				offspot.append(list(np.where(motheroff == Fspots[j])[0]))
			offspot = sum(offspot,[])
			offspot = np.asarray(offspot) # This is the age specific females that had the offspring in the offspring list.
			if len(offspot) == 0:
				tempoffspot = offspot
			else:
				tempoffspot = np.asarray(offspring)[offspot]			
			tempoffspot = tempoffspot.tolist()
			
			# Shuffle the offspring list
			shuffle(tempoffspot)
			
			# Get the number of survivors/deaths
			#offsurvivors = int(round((correctedS0)*len(offspot)))
			offsurvivors = int(round(correctedS0))			
						
			if equalsexratio == 'AtBirth':
				print('Not operating currently with special case of offno=6 for equalsexratio=AtBirth.')
				sys.exit(-1)
				offdeaths = int(len(offspot) - offsurvivors)
				# Sort by sex
				tempoffspot.sort(key=lambda x: x[4])
				# The remove first half of offsurvive and last half
				temp = tempoffspot[(offdeaths/2):]
				temp = temp[0:len(temp)-(offdeaths/2)]
				tempoffspring.append(temp)
			else:		
				# Grab the survived offspring location
				tempoffspring.append(tempoffspot[0:offsurvivors])
			
		# Flatten tempoffspring list
		flattened = [val for sublist in tempoffspring for val in sublist]
		tempoffspring = sum(tempoffspring,[])
		
		# Tracking
		Track_MOffDeaths.append((len(offspring)-len(tempoffspring))/2.)
		Track_FOffDeaths.append((len(offspring)-len(tempoffspring))/2.)
		
	else:
		
		# Get number of survivors for female and male
		offsex = np.asarray(offspring)[:,4]
		Foffsex = np.where(offsex == 0)[0] # Index location
		Moffsex = np.where(offsex == 1)[0] # Index location
		Foffsurvivors = int(round((s0_females)*len(Foffsex)))
		Moffsurvivors = int(round((s0_males)*len(Moffsex)))
		Foffdeaths = int(round((Fnewmortperc)*len(Foffsex)))
		Moffdeaths = int(round((Mnewmortperc)*len(Moffsex)))
					
		# Shuffle the offspring list
		shuffle(offspring)
		
		if equalsexratio == 'AtBirth':
			if (Foffdeaths + Moffdeaths) > 0:
				print('Warning: Equal sex ratio AtBirth is specified with different age 0 mortality values. Using total age 0 deaths.')
			offdeaths = Foffdeaths + Moffdeaths
			# Sort by sex
			offspring.sort(key=lambda x: x[4])
			# The remove first half of offsurvive and last half
			tempoffspring = offspring[(offdeaths/2):]
			tempoffspring = tempoffspring[0:len(tempoffspring)-(offdeaths/2)]
		else:
			# Shuffle Female and male index list
			shuffle(Foffsex)
			shuffle(Moffsex)
			
			# Grab the survived female and male index list 
			tempFoff = Foffsex[0:Foffsurvivors]
			tempMoff = Moffsex[0:Moffsurvivors]
			tempoff = np.concatenate((tempFoff,tempMoff),axis=0)
			# Grab the survived offspring location
			tempoffspring = np.asarray(offspring)[tempoff]
			# Back to list
			tempoffspring = tempoffspring.tolist()
			
		# Tracking
		Track_MOffDeaths.append(Moffdeaths)
		Track_FOffDeaths.append(Foffdeaths)
	
	return tempoffspring
	
	# End::DoOffMortality()
	
# ---------------------------------------------------------------------------------------------------	 
def DoOffspring(offno,lmbda,Bearpairs,CDpairs,\
Femalepercent,Births,infection,transmissionprob,equalsexratio,Mnewmortperc,Fnewmortperc,Track_MOffDeaths,Track_FOffDeaths,sigma,age,sex):
	'''
	DoOffspring()
	Choose numberof Offspring for each mated pair
	Input: selection choice for offspring number distrubution draw
	offno, Bearpairs, lmbda.
	Output: Bear Pairs + # of offspring [Female,Male,#Offspring]
	Also add # offspring to CDpairs
	'''
	
	# Function 1 is a uniform random draw between 0 and lmdba number	
	if (offno=='1'):
		
		tupDoOffRand = DoOffspringRandom(Bearpairs,CDpairs,lmbda,age)
		
		Bearpairs = tupDoOffRand[0]
		CDpairs = tupDoOffRand[1]
	
	# Function 2 is a Poisson draw
	elif (offno=='2'):
	
		tupDoOffPois = DoOffspringPoisson(Bearpairs,CDpairs,lmbda,age)
		
		Bearpairs = tupDoOffPois[0]
		CDpairs = tupDoOffPois[1]
		
	# Function 3 and 4 is a constant of lmbda offspring per each pairing
	elif (offno=='3' or offno=='4' or offno == '6'):
	
		tupDoOffConst = DoOffspringConstant(Bearpairs,CDpairs,\
		lmbda,age)
		
		Bearpairs = tupDoOffConst[0]
		CDpairs = tupDoOffConst[1]
	
	# Function 5 is a normal draw
	elif (offno=='5'):
		
		tupDoOff = DoOffspringNormal(Bearpairs,CDpairs,lmbda,sigma,age)
		Bearpairs = tupDoOff[0]
		CDpairs = tupDoOff[1]
	
	# Other functions
	else:
		print('This offspring birth rate option (offno) does not exist.')
		sys.exit(-1)
	
	# Assign sex of each offspring
	offspring = DoOffspringSex(Bearpairs,Femalepercent,CDpairs,equalsexratio)
	
	# If function 4 is choosen, then weed the offspring here
	if (offno == '4'):
		
		offspring = DoOffspringEqual(offspring,lmbda,age)
	
	# Store number of Births
	tempBirths = len(offspring)
	Births.append(tempBirths)
	
	# Apply Mortality to the egg class here
	offspring = DoOffMortality(offspring,Mnewmortperc,Fnewmortperc,equalsexratio,offno,age,lmbda,sex,Track_MOffDeaths,Track_FOffDeaths)
	
	# Assign infection to each offspring
	offspring = DoCDInfect(offspring,infection,transmissionprob)
		
	# Return variables from this argument
	tupDoOff = offspring,len(offspring)
	return tupDoOff
		
	# End::DoOffspring()