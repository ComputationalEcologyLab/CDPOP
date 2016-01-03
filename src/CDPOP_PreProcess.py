# -------------------------------------------------------------------------------------------------
# CDPOP_PreProcess.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file pre processing.
# --------------------------------------------------------------------------------------------------

# Import Modules with Except/Try statements

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."

# Scipy function KDTree
try:
	from scipy.spatial import KDTree
	scipyAvail = True
except ImportError:
	raise ImportError, "Scipy required."
	scipyAvail = False
	
# CDPOP functions
try:
	from CDPOP_Mate import *
	from CDPOP_Offspring import *
	from CDPOP_Disperse import *
except ImportError:
	raise ImportError, "CDPOP Modules missing."
	
# General python imports
import os,sys

# ----------------------------------------------------------------------------------
def loadFile(filename, header_lines=0, delimiter=None, cdpop_inputvars=False): ###
	'''
	Used to load file hearders according to current UNICOR standards
	as of 5/19/2011
	'''
	try:
		inputfile = open(filename)
	except (IOError,OSError) as e:
		print("Load file: %s the file (%s) is not available!"%(e,filename))
		sys.exit(-1)
	header_dict = {}
	data_list = []
	index_list = [] ###
	
	if delimiter != None:
		lines = [ln.rstrip().split(delimiter) for ln in inputfile]
	else:
		lines = [ln.rstrip().split() for ln in inputfile]
	# Close file
	inputfile.close()
		
	for i,line in enumerate(lines):
		if i < header_lines:
			if len(line) <= 1:
				print("Only one header value in line, skipping line...")
				continue
			#else:	
			elif len(line) == 2: ###
				header_dict[line[0]] = line[1]
			### working with a file where the first line is all header keys
			### and the following lines are their data
			elif cdpop_inputvars: ###
				for j in range(len(line)): ###
					header_dict[line[j]] = [] ###
					index_list.append(line[j]) ###
		else:
			#This is a fix to remove empty entries from from a line if they
			#exist, this is to fix a problem with reading in cdmatrices
			for i in range(line.count('')):
				line.remove('')
			data_list.append(line)
			if cdpop_inputvars: ###
				#tempTuple = ()
				for j in range(len(line)): ###
					# remove the following lines should the tuple representing bar-delimited values break anything -TJJ
					if line[j].find('|') != -1:
						tempList = line[j].split('|')
						line[j] = tuple(tempList)
					#---
					
					header_dict[index_list[j]].append(line[j]) ###
	
	if not cdpop_inputvars:
		return header_dict, data_list
	else:
		n_jobs = len(lines) - header_lines
		return header_dict, index_list, n_jobs

# ---------------------------------------------------------------------------------------------------	 
def GetMaxCDValue(threshold,cdmatrix):	
	'''
	GetMaxCDValue()
	This function calculates the maximum movement thresholds.
	'''	
	
	# movement threshold if max specified
	if str(threshold).endswith('max'):
		# If max
		if len(threshold.strip('max')) == 0:
			threshold = np.amax(cdmatrix)
		else:
			threshold = (float(threshold.strip('max'))/100.)\
			*np.amax(cdmatrix)
	else:
		threshold = float(threshold)
	
	return threshold
# ---------------------------------------------------------------------------------------------------	 
def ReadFitnessSurface(fitsurface):	
	'''
	ReadFitnessSurface()
	This function reads in the ascii fitness surface.
	'''	
	# Open file for reading
	inputfile = open(fitsurface,'r')

	# Read lines from the file
	lines = inputfile.readlines()

	#Close the file
	inputfile.close()

	# Create an empty matrix to append to
	values = []

	# Split up each line in file and append to empty matrix, x
	for i in lines:
		i = i.strip('\n').strip('\r').strip(' ')
		thisline = i.split(' ')
		values.append(thisline)
	
	# Grab some information from fitvalues: number of columns
	lenfirstline = len(values[0])
	ncols = int(values[0][lenfirstline-1])
		
	# Grab some information from values: number of rows
	lensecondline = len(values[1])
	nrows = int(values[1][lensecondline-1])

	# Grab some information from values: x coord value in lower left
	lenthirdline = len(values[2])
	xllcorner = float(values[2][lenthirdline-1])

	# Grab some information from values: y coord value in lower left
	lenfourthline = len(values[3])
	yllcorner = float(values[3][lenfourthline-1])

	# Grab some information from values: cell size
	lenfifthline =len(values[4])
	cellsize = float(values[4][lenfifthline-1])
	
	# Grab some information from values: Nodataval
	lensixthline =len(values[5])
	Nodataval = float(values[5][lensixthline-1])
	
	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)
	
	# Error statement
	if len(values)-6 != nrows or len(values[6]) != ncols:
		print('Spatial selection surface file is not in correct format. Return error number of dimensions.')
		sys.exit(-1)
	# And turn rasterfile into a list with float values without header information
	fitnessvalues = []
	for i in xrange(len(values)-6):
		fitnessvalues.append([])
		for j in xrange(len(values[6])):
			fitnessvalues[i].append(float(values[6+i][j]))
	
	# Return tuple
	tupReadFit = fitnessvalues,ncols,nrows,xll,yll,cellsize,Nodataval
	return tupReadFit
	# End::ReadFitnessSurface()
	
# ---------------------------------------------------------------------------------------------------	
def CreateAlleleList(loci,alleles,xgenes):
	'''
	CreateAlleleList()
	This function creates a list for the allele storage.
	'''	
	
	# Store all information in a list [loci][allele#,probability]
	allelst = []
	for i in xrange(loci):
		allelst.append([])
		for k in xrange(alleles[i]):
			allspot = alleles[i]*i+1+k
			allelst[i].append([int(k),float(xgenes[allspot][1])])
	
	# Return variables
	return allelst
	
	# End::CreateAlleleList()
	
# ---------------------------------------------------------------------------------------------------	 
def InitializeGenes(intgenesans,allefreqfilename,loci,alleles,datadir):	
	
	allelst = []
	
	# If genetic structure intialized by a file...
	if intgenesans == 'file':
		
		# Loop through allele frequency files given
		for ifile in xrange(len(allefreqfilename)):
			
			fileans = allefreqfilename[ifile]
		
			# If genetic structure intialized by a file...
			if fileans != 'random':
		
				# Check statements
				if os.path.exists(datadir+fileans):
					# Open file for reading
					inputfile = open(datadir+fileans,'rU')
				else:
					print("CDPOP InitializeGenes() error: open failed, could not open %s"%(fileans))
					sys.exit(-1)
						
				# Read lines from the file
				lines = inputfile.readlines()
				
				#Close the file
				inputfile.close()
				
				# Create an empty matrix to append to
				xgenes = []
			
				# Split up each line in file and append to empty matrix, x
				for i in lines:
					thisline = i.strip('\n').strip('\r').strip(' ').split(',')
					xgenes.append(thisline)
				
				# Error check here
				if (len(xgenes)-1) != sum(alleles):
					print('Allele frequency file is not the specified number of loci and alleles as in in input file.')
					sys.exit(-1)
				
				# Delete lines from earlier
				del(lines)
			
				# Call CreateAlleleList()
				allelst.append(CreateAlleleList(loci,alleles,xgenes))
			
			# If genetic structure is to be initialize by random
			elif fileans == 'random':
				
				# Create even distribution
				xgenes = []
				xgenes.append(['Allele List','Frequency'])
				for iloci in xrange(loci):
					for iall in xrange(alleles[iloci]):
						xgenes.append(['L'+str(iloci)+'A'+str(iall),str(1.0/alleles[iloci])])
				
				# Call CreateAlleleList()
				allelst.append(CreateAlleleList(loci,alleles,xgenes))	
	
	# If just random for entire population
	elif intgenesans == 'random':
		
		# Create even distribution
		xgenes = []
		xgenes.append(['Allele List','Frequency'])
		for iloci in xrange(loci):
			for iall in xrange(alleles[iloci]):
				xgenes.append(['L'+str(iloci)+'A'+str(iall),str(1.0/alleles[iloci])])
		
		# Call CreateAlleleList()
		allelst.append(CreateAlleleList(loci,alleles,xgenes))
	
		# Delete x variable
		del(xgenes)
	
	# Error statement
	else:
		allelst = []
			
	# Return variables
	return allelst
		
	# End::InitializeGenes()
	
# ---------------------------------------------------------------------------------------------------	 
def InitializeAge(agefilename,nogrids,datadir):
	'''
	InitializeAge()
	This function initializes the age of each population
	with an age distribution list.
	'''
	
	# Loop through age files given
	agelst = []
	for ifile in xrange(len(agefilename)):
		agelst.append([])
		
		fileans = agefilename[ifile]			
				
		# Check statements
		if os.path.exists(datadir+fileans):
			# Open file for reading
			inputfile = open(datadir+fileans,'rU')
		else:
			print("CDPOP InitializeAge() error: open failed, could not open %s"%(datadir+fileans))
			sys.exit(-1)
		
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Error statement on agevars file length
		if len(lines) < 3:
			print('Agevars.csv file must have at least 2 age classes: Age 0 (eggs) and Age 1 (Adults).')
			sys.exit(-1)
		
		# Create an empty matrix to append to
		xage = []
		
		# Split up each line in file and append to empty matrix, x
		for i in lines:
			thisline = i.strip('\n').strip('\r').strip(' ').split('\r')
			for j in xrange(len(thisline)):
				xage.append(thisline[j].split(','))
			
		# Delete lines from earlier
		del(lines)
		
		# Store all information in a list [age,probability]
		ageclass = []
		ageno = []
		agemort = []
		eggs_mean = []
		eggs_sigma = []	
		Mmature = []
		Fmature = []
		for i in xrange(len(xage)-2): # i+2, only get Age 1+
			ageclass.append(int(xage[i+2][0]))
			ageno.append(float(xage[i+2][1]))
			if len(xage[i+2][2].split('|')) == 1: 
				agemort.append(float(xage[i+2][2])/100.)
			else:
				agemort.append(xage[i+2][2].split('|'))
			if len(xage[i+2][3].split('|')) == 1: 
				eggs_mean.append(float(xage[i+2][3]))
			else:
				eggs_mean.append(xage[i+2][3].split('|'))
			if len(xage[i+2][4].split('|')) == 1: 
				eggs_sigma.append(float(xage[i+2][4]))
			else:
				eggs_sigma.append(xage[i+2][4].split('|'))
			Mmature.append(float(xage[i+2][5]))
			Fmature.append(float(xage[i+2][6]))
			
		# Get age distribution list
		for i in xrange(len(ageno)):
			agelst[ifile].append([ageclass[i],ageno[i]/sum(ageno)])
		
		# Get Age 0 mortality
		if len(xage[1][2].split('|')) == 1: 
			newmortperc = float(xage[1][2])/100.
		else:
			newmortperc = xage[1][2].split('|')
			
		# Error checks here: if number of classes does not equal mortality age classes
		if len(agelst[ifile]) != len(agemort):
			print('Agedistribution data not fully entered correctly.')
			sys.exit(-1)
		# Error check that age 0s are initialized
		if xage[1][1] != '0':
			print('Warning: Age 0 was initialized at start of model, enter 0 for distribution in Agevars.csv file.')
			sys.exit(-1)
		# Deletes
		del(xage)
	
	# Return variables
	return agelst,agemort,eggs_mean,eggs_sigma,newmortperc,Mmature,Fmature
	
	# End::InitializeAge()

# ---------------------------------------------------------------------------------------------------	 
def ReadXY(xyfilename):
	'''
	ReadMateXYCDMatrix()
	This function reads in the xy values for the cost distance matrix.
	'''	
	
	
	# Open file for reading
	inputfile = open(xyfilename+'.csv','rU')
	
	# Read lines from the file
	lines = inputfile.readlines()
	
	#Close the file
	inputfile.close()
	
	# Create an empty matrix to append to
	xy = []
	
	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.rstrip('\n').rstrip('\r').split(',')
		xy.append(thisline)
		
	# Delete lines from earlier
	del(lines)
	
	# Return variables
	return xy
	
	# End::ReadXY()

# ---------------------------------------------------------------------------------------------------	 
def ReadCDMatrix(cdmatrixfilename,function,threshold,A,B,C,subpop):
	'''
	ReadMateCDMatrix()
	This function reads in the mating cost distance matrix.
	'''	
	
	# Check statements
	if os.path.exists(cdmatrixfilename+'.csv'):
		# Open file for reading
		inputfile = open(cdmatrixfilename+'.csv','rU')
	else:
		print("CDPOP ReadCDMatrix() error: open failed, could not open %s"%(cdmatrixfilename+'.csv'))
		sys.exit(-1)
	
	# Read lines from the file
	lines = inputfile.readlines()
	
	# Close the file
	inputfile.close()
	
	# Create an empty matrix to append to 
	bigCD = []
	
	# Split up each line in file and append to empty matrix, x
	for spot in lines:
		thisline = spot.strip('\n').split(',')
		bigCD.append(thisline[0:len(lines)])
	bigCD = np.asarray(bigCD,dtype='float')
	
	# Delete lines from earlier
	del(lines)
		
	# Store number of files
	nofiles = len(bigCD)
	
	# Calculate max and min of bigCD matrix
	minbigCD = np.amin(bigCD)
	maxbigCD = np.amax(bigCD)
			
	# Get maximum cdvalue to use for movethreshold if specified
	threshold = GetMaxCDValue(threshold,bigCD)
	
	# Create a matrix of to be filled 
	cdmatrix = []
	
	# Fill up matrix with float value of array x
	for j in xrange(nofiles):
		cdmatrix.append([])
		for k in xrange(nofiles):
			
			# For the linear function
			if function == '1':
				scale_min = 0.
				scale_max = threshold
				# Set = to 0 if cdvalue is greater than movethreshold
				if float(bigCD[j][k]) > threshold:
					cdmatrix[j].append(0.0)
				# If threshold is 0 (philopatry) set to 1 - can't dived by 0
				elif float(bigCD[j][k]) <= threshold and threshold == 0.0:
					cdmatrix[j].append(1.0)
				# Else calculated function value and if not philopatry
				elif float(bigCD[j][k]) <= threshold and threshold != 0.0:
					cdmatrix[j].append(-(1./threshold)*float(bigCD[j][k]) + 1)
				else:
					print('Something off in linear function values.')
					sys.exit(-1)
						
			# For the inverse square function
			elif function == '2':
				
				# This function gets rescale: calculate here
				if threshold == 0:
					scale_min = 0.
				else:
					scale_min = 1./(pow(threshold,2))
				scale_max = 1.
				
				# Set = to 0 if cdvalue is greater than movethreshold
				if float(bigCD[j][k]) > threshold:
					cdmatrix[j].append(0.0)
				# If threshold is 0 (philopatry) set to 1 - can't dived by 0
				elif float(bigCD[j][k]) <= threshold and threshold == 0.0:
					cdmatrix[j].append(1.0)
				# If cd mat is 0. 
				elif float(bigCD[j][k]) <= threshold and threshold != 0.0 and float(bigCD[j][k]) == 0.0 or (minbigCD == maxbigCD or int(maxbigCD) == 0):
					cdmatrix[j].append(1.0)
				# Else calculated function value
				elif float(bigCD[j][k]) <= threshold and threshold != 0.0 and float(bigCD[j][k]) != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
					invsq_val = 1./(pow(float(bigCD[j][k]),2))
					invsq_val = (invsq_val - scale_min) / (scale_max - scale_min)
					cdmatrix[j].append(invsq_val)# Else something else.
				else:
					print('Something off in inv squ function values.')
					sys.exit(-1)
					
							
			# Nearest neighbor function here
			elif function == '3':
				print('Nearest neighbor function is not currently implemented.')
				print('You can use Linear function with neighbor threshold for approximation. Email Erin.')
				sys.exit(-1)
				
			# Random function here
			elif function == '4':
				scale_min = 0.
				scale_max = threshold
				# Set = to 0 if cdvalue is greater than movethreshold
				if float(bigCD[j][k]) > threshold:
					cdmatrix[j].append(0.0)
				# Else calculated function value
				else:
					cdmatrix[j].append(1.0)
			
			# For the negative binomial function
			elif function == '5':
			
				# This function gets rescale: calculate here
				scale_min = A*pow(10,-B*float(threshold))
				scale_max = A*pow(10,-B*float(minbigCD))
			
				# Set = to 0 if cdvalue is greater than movethreshold
				if float(bigCD[j][k]) > threshold:
					cdmatrix[j].append(0.0)
				# Rescaled value divide by zero check cases
				elif float(bigCD[j][k]) <= threshold and threshold == 0.0 and (minbigCD == maxbigCD or int(maxbigCD) == 0):
					cdmatrix[j].append(1.0)
				# Else calculated function value
				elif float(bigCD[j][k]) <= threshold and threshold != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
					negexp = A*pow(10,-B*float(bigCD[j][k]))
					negexp = (negexp - scale_min) / (scale_max - scale_min)
					cdmatrix[j].append(negexp)
				# Else something else.
				else:
					print('Something off in neg exp function values.')
					sys.exit(-1)
					
			# For in a subpopulation only
			elif function == '6':
				
				scale_min = 0.
				scale_max = 1.
				# Check if within the same subpopulation
				if subpop[j] == subpop[k]:
					cdmatrix[j].append(1.0)
				else:
					cdmatrix[j].append(0.0)
			
			# For Gaussian function 
			elif function == '7':
			
				# This function gets rescale: calculate here
				scale_min = A*np.exp(-((float(threshold)-B)**2)/(2*C**2))
				scale_max = A*np.exp(-((float(minbigCD)-B)**2)/(2*C**2))
			
				# Set = to 0 if cdvalue is greater than movethreshold
				if float(bigCD[j][k]) > threshold:
					cdmatrix[j].append(0.0)
				# Rescaled value divide by zero check cases
				elif float(bigCD[j][k]) <= threshold and threshold == 0.0 and (minbigCD == maxbigCD or int(maxbigCD) == 0):
					cdmatrix[j].append(1.0)
				# Else calculated function value
				elif float(bigCD[j][k]) <= threshold and threshold != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
					gauss_val = A*np.exp(-((float(bigCD[j][k])-B)**2)/(2*C**2))
					gauss_val = (gauss_val - scale_min) / (scale_max - scale_min)
					cdmatrix[j].append(gauss_val)
				# Else something else.
				else:
					print('Something off in gauss function values.')
					sys.exit(-1)
					
			# For cost distance matrix only function 
			elif function == '8':
				
				scale_min = minbigCD
				scale_max = threshold
				
				# Set = to 0 if cdvalue is greater than movethreshold
				if float(bigCD[j][k]) > threshold:
					cdmatrix[j].append(0.0) 
				# Rescaled value divide by zero check cases - philopatry
				elif (float(bigCD[j][k]) <= threshold and threshold == 0.0) and (minbigCD == maxbigCD or int(maxbigCD) == 0 or threshold == minbigCD):
					cdmatrix[j].append(1.0)
				# If cd mat is 0. 
				elif (float(bigCD[j][k]) <= threshold and threshold != 0.0 and float(bigCD[j][k]) == 0.0) or (minbigCD == maxbigCD or int(maxbigCD) == 0 or threshold == minbigCD):
					cdmatrix[j].append(1.0)
				# Else calculated function value
				elif (float(bigCD[j][k]) <= threshold and threshold != 0.0 and float(bigCD[j][k]) != 0.0) and (minbigCD != maxbigCD and int(maxbigCD) != 0 and threshold != minbigCD):
					cd_val = (float(bigCD[j][k]) - scale_min) / (scale_max - scale_min)
					cdmatrix[j].append(1. - cd_val)
				# Else something else.
				else:
					print('Something off in 8 function values.')
					sys.exit(-1)
					
			# error
			else:
				print('This movement function option does not exist.')
				sys.exit(-1)
	
	# Delete variables
	del(bigCD)
	
	# Return variables
	tupReadMat = cdmatrix,threshold,scale_min,scale_max
	return tupReadMat
	
	# End::ReadMateCDMatrix

# ---------------------------------------------------------------------------------------------------	 
def DoCDClimate(datadir,icdtime,cdclimgentime,matecdmatfile,dispcdmatfile,matemoveno,Fdispmoveno,Mdispmoveno,matemovethresh,Fdispmovethresh,Mdispmovethresh,matemoveparA,matemoveparB,matemoveparC,FdispmoveparA,FdispmoveparB,FdispmoveparC,MdispmoveparA,MdispmoveparB,MdispmoveparC,subpop,agemort,offno,lmbda,sigma,K_envvals,newmortperc,fitvals):
	'''
	DoCDClimate()
	Reads in cost distance matrices and converts to probabilities.
	'''
	
	# -------------------------------
	# Extract cdclimate values here
	# -------------------------------
	# Store cdmat file information - header file (loadFile()) passes tuple or string if only 1
	if not isinstance(cdclimgentime, (list,tuple)):
		# Error checks
		if cdclimgentime != ['0']:
			print('If not using CDClimate option, set begin time loop with cdclimgentime at 0.')
			sys.exit(-1)
	
	# Get age mortality here
	tempmort = []
	templmbda = []
	tempsigma = []
	for i in xrange(len(agemort)):
		if not isinstance(agemort[i], float):
			tempmort.append(float(agemort[i][icdtime])/100.)
		else:
			tempmort.append(agemort[i])
		if not isinstance(lmbda[i], float):
			templmbda.append(float(lmbda[i][icdtime]))
		else:
			templmbda.append(lmbda[i])
		if not isinstance(sigma[i], float):
			tempsigma.append(float(sigma[i][icdtime]))
		else:
			tempsigma.append(sigma[i])
	if not isinstance(agemort[i], float):
		tempnewmort = float(newmortperc[icdtime])/100.
	else:
		tempnewmort = newmortperc
	# Get birth rate values here, r and K_env
	if isinstance(offno, (list,tuple)):
		tempoffno = offno[icdtime]
	else:
		tempoffno = offno
	if isinstance(K_envvals, (list,tuple)):
		tempK_env = int(K_envvals[icdtime])
	else:
		tempK_env = int(K_envvals)
		
	# Check for each instance
	if isinstance(matecdmatfile, (list,tuple)):
		matecdmatfile = datadir+matecdmatfile[icdtime]
	else:	
		matecdmatfile = datadir+matecdmatfile
	if isinstance(dispcdmatfile, (list,tuple)):
		dispcdmatfile = datadir+dispcdmatfile[icdtime]
	else:
		dispcdmatfile = datadir+dispcdmatfile
	if isinstance(matemoveno, (list,tuple)):
		matemoveno = matemoveno[icdtime]		
	if isinstance(Fdispmoveno, (list,tuple)):
		Fdispmoveno = Fdispmoveno[icdtime]
	if isinstance(Mdispmoveno, (list,tuple)):
		Mdispmoveno = Mdispmoveno[icdtime]
	if isinstance(matemovethresh, (list,tuple)):
		matemovethresh = matemovethresh[icdtime]
	if isinstance(Fdispmovethresh, (list,tuple)):
		Fdispmovethresh = Fdispmovethresh[icdtime]
	if isinstance(Mdispmovethresh, (list,tuple)):
		Mdispmovethresh = Mdispmovethresh[icdtime]
	if isinstance(matemoveparA, (list,tuple)):
		matemoveparA = float(matemoveparA[icdtime])
	else:
		matemoveparA = float(matemoveparA)
	if isinstance(matemoveparB, (list,tuple)):
		matemoveparB = float(matemoveparB[icdtime])
	else:
		matemoveparB = float(matemoveparB)		
	if isinstance(matemoveparC, (list,tuple)):
		matemoveparC = float(matemoveparC[icdtime])
	else:
		matemoveparC = float(matemoveparC)
	if isinstance(FdispmoveparA, (list,tuple)):			
		FdispmoveparA = float(FdispmoveparA[icdtime])
	else:
		FdispmoveparA = float(FdispmoveparA)
	if isinstance(FdispmoveparB, (list,tuple)):
		FdispmoveparB = float(FdispmoveparB[icdtime])
	else:
		FdispmoveparB = float(FdispmoveparB)		
	if isinstance(FdispmoveparC, (list,tuple)):
		FdispmoveparC = float(FdispmoveparC[icdtime])
	else:
		FdispmoveparC = float(FdispmoveparC)		
	if isinstance(MdispmoveparA, (list,tuple)):
		MdispmoveparA = float(MdispmoveparA[icdtime])
	else:
		MdispmoveparA = float(MdispmoveparA)
	if isinstance(MdispmoveparB, (list,tuple)):
		MdispmoveparB = float(MdispmoveparB[icdtime])
	else:
		MdispmoveparB = float(MdispmoveparB)		
	if isinstance(MdispmoveparC, (list,tuple)):
		MdispmoveparC = float(MdispmoveparC[icdtime])
	else:
		MdispmoveparC = float(MdispmoveparC)
	# Patch based parameters	
	tempfitvals = []
	for isub in xrange(len(subpop)):
		if len(fitvals) > 0:
			tempfitvals.append([])
			for i in xrange(len(fitvals[isub])):
				if len(fitvals[isub][i].split('|')) > 1:
					tempfitvals[isub].append(fitvals[isub][i].split('|')[icdtime])
				else:
					tempfitvals[isub].append(fitvals[isub][i])
						
	# ---------------------------------------------------------
	# Read in cdmatrix.csv and convert to a probability matrix
	# ---------------------------------------------------------
	
	# If mate and disp are the same, then only read in once.
	if (matecdmatfile == dispcdmatfile) \
	and (Fdispmoveno == Mdispmoveno == matemoveno) \
	and (Fdispmovethresh == Mdispmovethresh == matemovethresh):
		tupReadMat = ReadCDMatrix(matecdmatfile,matemoveno,\
		matemovethresh,matemoveparA,matemoveparB,matemoveparC,subpop)
		
		# Unpack tuple
		matecdmatrix = np.asarray(tupReadMat[0])
		matemovethresh = tupReadMat[1]
		mate_ScaleMin = tupReadMat[2]
		mate_ScaleMax = tupReadMat[3]
		
		# Then Set disp = mate
		Fdispcdmatrix = matecdmatrix
		Mdispcdmatrix = matecdmatrix
		Fdispmovethresh = matemovethresh
		Mdispmovethresh = matemovethresh
		Fdisp_ScaleMin = mate_ScaleMin
		Fdisp_ScaleMax = mate_ScaleMax
		Mdisp_ScaleMin = mate_ScaleMin
		Mdisp_ScaleMax = mate_ScaleMax		

	# Else if anything is different	
	else: 
		# ---------------------------------------
		# Read in cdmatrix.csv - For Mating
		# ---------------------------------------	
		tupReadMat = ReadCDMatrix(matecdmatfile,matemoveno,\
		matemovethresh,matemoveparA,matemoveparB,matemoveparC,subpop)
		matecdmatrix = np.asarray(tupReadMat[0])
		matemovethresh = tupReadMat[1]
		mate_ScaleMin = tupReadMat[2]
		mate_ScaleMax = tupReadMat[3]
	
		# ------------------------------------------------
		# Read in cdmatrix.csv - For Female Dispersal 
		# ------------------------------------------------	
		tupReadMat = ReadCDMatrix(dispcdmatfile,Fdispmoveno,\
		Fdispmovethresh,FdispmoveparA,FdispmoveparB,FdispmoveparC,subpop)
		Fdispcdmatrix = np.asarray(tupReadMat[0])
		Fdispmovethresh = tupReadMat[1]
		Fdisp_ScaleMin = tupReadMat[2]
		Fdisp_ScaleMax = tupReadMat[3]

		# ----------------------------------------------
		# Read in cdmatrix.csv - For Male Dispersal
		# ----------------------------------------------	
		tupReadMat = ReadCDMatrix(dispcdmatfile,Mdispmoveno,\
		Mdispmovethresh,MdispmoveparA,MdispmoveparB,MdispmoveparC,subpop)
		Mdispcdmatrix = np.asarray(tupReadMat[0])
		Mdispmovethresh = tupReadMat[1]
		Mdisp_ScaleMin = tupReadMat[2]
		Mdisp_ScaleMax = tupReadMat[3]
	
	# Return this functions variables
	tupClimate = matecdmatrix,Fdispcdmatrix,Mdispcdmatrix,matemovethresh,\
	Fdispmovethresh,Mdispmovethresh,Fdisp_ScaleMin,Fdisp_ScaleMax,Mdisp_ScaleMin,Mdisp_ScaleMax,mate_ScaleMin,mate_ScaleMax,tempmort,tempoffno,templmbda,tempsigma,tempK_env,tempnewmort,tempfitvals,matemoveno,Fdispmoveno,Mdispmoveno	
	return tupClimate
	
	#End::DoCDClimate()
	
# ---------------------------------------------------------------------------- 
def InitializeInfect(cdinfect,Infected,nogrids,sex):
	'''
	InitializeInfect()
	This function initializes the infection status.
	'''
	
	# Create empty list to append to 
	infection = []
	tempinf = []
	
	# If cdinfect...
	if cdinfect == 'Y':
		
		# Assign initial individual a random infection status
		for ithinfect in xrange(nogrids):
		
			# If a NA value
			if sex[ithinfect] == 'NA':
				# Append a NA
				infection.append('NA')
			else:
				# Append a random number 0 or 1
				infection.append(int(2*rand()))
				tempinf.append(infection[ithinfect])
	
	# If not cdinfect...
	else: 
		
		# Assign initial individual a random infection status
		for ithinfect in xrange(nogrids):
		
			# If a NA value
			if sex[ithinfect] == 'NA':
				# Append a NA
				infection.append('NA')
			else:			
				# Append 0
				infection.append(0)
				tempinf.append(infection[ithinfect])
			
	# Append to Infected
	Infected.append(sum(tempinf))
	del(tempinf)
	
	tupInf = infection,Infected
	return tupInf
	
	# End::InitializeInfect()

# ---------------------------------------------------------------------------- 
def DoGridOut_cdpop0(ithmcrundir,gen,loci,alleles,nogrids,subpop,xgrid,ygrid,\
id,sex,age,agelst,genes,intgenesans,infection,allelst,subpopmigration,subpopemigration,geneswap):
	'''
	DoGridOut_cdpop0()
	Output grid0.csv in cdpop format
	'''			
	# Create file to write matrix to
	outputfile = open(ithmcrundir+'grid'+str(0)+'.csv','w')
		
	# Write out the titles
	title = ['Subpopulation','XCOORD','YCOORD','ID','sex','age','infection','DisperseCDist']
			
	# Write out the title from xy points
	for i in xrange(len(title)):
		outputfile.write(title[i]+',')
	
	# Write out the genes title informations
	# Loop through loci
	for i in xrange(loci-1):
		
		# Loop for allele length
		for j in xrange(alleles[i]):
			outputfile.write('L'+str(i)+'A'+str(j)+',')
	
	# To get a return character on the end of the title
	for i in xrange(alleles[loci-1]-1):
		outputfile.write('L'+str(loci-1)+'A'+str(i)+',')
	outputfile.write('L'+str(loci-1)+'A'+str(alleles[loci-1]-1))
	# Get return character
	outputfile.write('\n')
	
	# Write out all of the information.		
	for i in xrange(nogrids):
		outputfile.write(subpop[i]+',')
		outputfile.write(str(float(xgrid[i]))+',')
		outputfile.write(str(float(ygrid[i]))+',')
		if sex[i] == 'NA':
			outputfile.write('OPEN,') #id
			outputfile.write('NA,') # sex
			outputfile.write('NA,') # age
			if intgenesans != 'known':
				age.append('NA')
			outputfile.write('NA,') # infection
			outputfile.write('NA,') # dispersalDist
		else:
			outputfile.write(id[i]+',')
			outputfile.write(str(sex[i])+',')
			if intgenesans == 'known':
				outputfile.write(age[i]+',')			
			else:
				if len(agelst) > 1:
					agetemp = w_choice_general(agelst[int(subpop[i])-1])[0]
				else:
					agetemp = w_choice_general(agelst[0])[0]
				outputfile.write(str(agetemp)+',')
				age.append(agetemp)
			outputfile.write(str(infection[i])+',')
			outputfile.write('Initial,')
		
		# if known genes 
		if intgenesans == 'known' and geneswap == 0:
			
			# Write out gene info
			for jk in xrange(loci-1):
				for kl in xrange(alleles[jk]):
					if sex[i] == 'NA':
						outputfile.write('NA,')
					else:
						outputfile.write(str(genes[i][jk][kl])+',')
			# To get return character on end
			for jk in xrange(alleles[loci-1]-1):
				if sex[i] == 'NA':
					outputfile.write('NA,')
				else:
					outputfile.write(str(genes[i][loci-1][jk])+',')
			if sex[i] == 'NA':
				outputfile.write('NA\n')
			else:
				outputfile.write(str(genes[i][loci-1][alleles[loci-1]-1])+'\n')
			
		# NA values for when geneswap > gen
		elif intgenesans != 'known' and geneswap != 0:
			
			# Write out gene info
			for jk in xrange(loci-1):
				for kl in xrange(alleles[jk]):
					if sex[i] == 'NA':
						outputfile.write('NA,')
					else:
						outputfile.write(str(genes[i][jk][kl])+',')
			# To get return character on end
			for jk in xrange(alleles[loci-1]-1):
				if sex[i] == 'NA':
					outputfile.write('NA,')
				else:
					outputfile.write(str(genes[i][loci-1][jk])+',')
			if sex[i] == 'NA':
				outputfile.write('NA\n')
			else:
				outputfile.write(str(genes[i][loci-1][alleles[loci-1]-1])+'\n')
			
		# if file genes or random genes 			
		elif intgenesans != 'known' and geneswap == 0:
				
			#Store empty array to be appended to for gene info
			indall = []
						
			# And store genes information
			genes.append([])
						
			# For each loci:
			for j in xrange(loci):
			
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
				for k in xrange(alleles[j]):
					
					# Somebody not in this spot
					if sex[i] == 'NA':
						# THen append tempindall to indall
						indall.append('NA')
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
								
						# THen append tempindall to indall
						indall.append(tempindall)						
						# And to genes list
						genes[i][j].append(tempindall)
		
			# Add indall information to outputfile text
			for j in xrange(len(indall)-1):
				outputfile.write(str(indall[j])+',')				
			# To get return character on the end
			outputfile.write(str(indall[len(indall)-1])+'\n')
			
	# Here add unique number of subpop spot in list (make list of lists)
	nosubpops = len(np.unique(subpop))
	unique_subpops = np.unique(subpop)
	# THen add spot in empty
	subpopmigration.append([])
	subpopemigration.append([])
	for i in xrange(nosubpops):
		subpopmigration[0].append([0])
		subpopemigration[0].append([0])
	
	# Close file
	outputfile.close()
	
	# Return variables
	return genes,subpopmigration,subpopemigration
	
	# End::DoGridOut_cdpop0()	
	
# ---------------------------------------------------------------------------------------------------	 
def DoPreProcess(outdir,ibatch,ithmcrun,\
xyfilename,agefilename,equalsexratio,\
loci,intgenesans,allefreqfilename,alleles,gen,logfHndl,\
cdevolveans,cdinfect,Infected,
subpopmigration,subpopemigration,datadir,geneswap):
	'''
	DoPreProcess()
	This function does all the pre-processing work before
	CDPOP begins its time loops.
	'''
	# ----------------------------
	# Create directory
	# ----------------------------		
	ithmcrundir = outdir+'batchrun'+\
	str(ibatch)+'mcrun'+str(ithmcrun)+'/'
	os.mkdir(ithmcrundir)
	
	# ------------------------------------------------------------------
	# Read in xy points file and store info in list
	# ------------------------------------------------------------------
	xy = ReadXY(xyfilename)
	
	# Error statement for 5 column data
	if len(xy[1]) != 17 and intgenesans!='known':
		print('XY input file must be 17 columns, see example input files.')
		sys.exit(-1)
	
	# Store all information in lists by variable name
	FID = []
	subpop = []
	xgrid = []
	ygrid=[]
	id = []
	sex = []	
	age = []	
	genes = []
	infection = []
	fitvals = []
	for i in xrange(len(xy)-1):
		FID.append(i)
		subpop.append(xy[i+1][0])
		if xy[i+1][1] == 'NA' or xy[i+1][2] == 'NA':
			print('You need to specify the (x,y) locations for this location even if it is a NA value.')
			sys.exit(-1)
		if xy[i+1][3] == 'NA' and equalsexratio == 'WrightFisher':
			print('You can not force equal sex ratio (WrightFisher) with fluctuating population size.')
			print('Specify N for Equal Sex Ratio option.')
			sys.exit(-1)
		xgrid.append(float(xy[i+1][1]))
		ygrid.append(float(xy[i+1][2]))
		id.append(xy[i+1][3])
		if intgenesans != 'known':
			if cdevolveans == '1' or cdevolveans == '3':
				fitvals.append([xy[i+1][5],xy[i+1][6],xy[i+1][7]])
			elif cdevolveans == '2':
				fitvals.append([xy[i+1][8],xy[i+1][9],xy[i+1][10],xy[i+1][11],xy[i+1][12],xy[i+1][13],xy[i+1][14],xy[i+1][15],xy[i+1][16]])
		else:
			if cdevolveans == '1' or cdevolveans == '3':
				fitvals.append([xy[i+1][8],xy[i+1][9],xy[i+1][10]])
			elif cdevolveans == '2':
				fitvals.append([xy[i+1][11],xy[i+1][12],xy[i+1][13],xy[i+1][14],xy[i+1][15],xy[i+1][16],xy[i+1][17],xy[i+1][18],xy[i+1][19]])
		
		# Only grab sex information from the file is equal sex ratio is N
		if equalsexratio == 'N' or equalsexratio == 'AtBirth':
			if xy[i+1][3] == 'NA' and xy[i+1][4] != 'NA':
				print('Must specify NA for the sex value of all empty starting locations.')
				sys.exit(-1)
			sex.append(xy[i+1][4])
				
			# Change id to 'OPEN'
			if sex[i] == 'NA':
				id[i] = 'OPEN'
			
			# If sex was F or M...change to 0 and 1
			if sex[i] == 'F':
				sex[i] = '0'
			elif sex[i] == 'M':
				sex[i] = '1'
		
		# Known genes check and switch to randomo
		if intgenesans == 'known' and geneswap != 0:
			print('Warning: Known gene file can not be used with gene start time > 0 Will continue but with random gene assignment.')
			intgenesans = 'random'
		
		# Get known information from file 
		elif intgenesans == 'known' and geneswap == 0:
			if len(xy[1]) != 20+sum(alleles):
				print('Specified known genetic initializtion file. Not correct format. See example files for number of fields needed.')
				sys.exit(-1)
			# Age storage here for known file
			age.append(xy[i+1][5])
			# Error check for age in known file must be 1 or greater
			if age[i] != 'NA':
				if int(age[i]) < 1:
					print('Known file must initize with age 1+.')
					sys.exit(-1)	
			
			# genes[individual][locus][allele]
			genes.append([])			
			# Error check here to make sure gene file matches specified loci and alleles
			if sum(alleles) != len(xy[i+1][20:len(xy[i+1])]):
				print('Known genes file does not match loci and alleles given.')
				sys.exit(-1)							
			for j in xrange(loci):
				genes[i].append(xy[i+1][int(20+sum(alleles[0:j])):int(20+sum(alleles[0:j+1]))])	
		
		# For other cases when geneswap is great than 0
		if intgenesans != 'known' and geneswap != 0:
			# genes[individual][locus][allele]
			genes.append([])
			for j in xrange(loci):
				genes[i].append([])
				for iall in xrange(alleles[j]):
					genes[i][j].append('NA')
			
	# Store the number of grids
	nogrids = len(xy)-1
	
	# Delete x variable
	del(xy)
	
	# --------------------------
	# Error Checks
	# --------------------------	
	
	# For now, subpops need to be ordered 1 to N and not skipping, no 0s
	if len(np.where(np.unique(subpop)=='0')[0]) != 0:
		print('Subpopulation identification field can not have 0 values.')
		sys.exit(-1)
	tempcheck = []
	for i in xrange(len(np.unique(subpop))):
		tempcheck.append(int(np.unique(subpop)[i]))
	tempcheck = np.sort(tempcheck)
	if len(tempcheck) > 1:
		for i in xrange(len(tempcheck)-1):
			if tempcheck[i+1]-tempcheck[i] > 1:
				print('Subpopulation identification field must be labeled sequentially or a single value.')
				sys.exit(-1)	
		
	# If equal sex ratio is Y, then split up sex into equal parts
	if equalsexratio == 'WrightFisher':
	
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
			
			# Then create half males and females and shuffle
			sextemp = np.append(np.zeros(subgridtotal[i]/2,"int"),np.ones(subgridtotal[i]/2,"int"))
			np.random.shuffle(sextemp)
			# Loop through these individuals and append to LIST sex
			for j in xrange(len(sextemp)):			
				# The add them together and shuffle and append to sex
				sex.append(str(sextemp[j]))
				
		# Delete extra stuff
		del(sextemp)
		
	# -------------------------------------------
	# Initialize age structure
	# ------------------------------------------- 
	agelst,agemort,egg_lmbdavals,egg_sigmavals,newmortperc,Mmature,Fmature = InitializeAge(agefilename,nogrids,datadir)
	
	# ------------------------------------------
	# Initialize infection
	# ------------------------------------------
	tupInf = InitializeInfect(cdinfect,Infected,nogrids,sex)
	infection = tupInf[0]
	Infected = tupInf[1]
	
	# ----------------------------------------------
	# Initialize genetic structure
	# ----------------------------------------------
	allelst = InitializeGenes(intgenesans,allefreqfilename,loci,alleles,datadir)
	
	# --------------------------------------------------------------------
	# Create output file grid0.csv and write to it and return genes
	# -------------------------------------------------------------------- 
	genes,subpopmigration,subpopemigration = DoGridOut_cdpop0(ithmcrundir,0,loci,alleles,\
	nogrids,subpop,xgrid,ygrid,id,sex,age,agelst,genes,intgenesans,\
	infection,allelst,subpopmigration,subpopemigration,geneswap)
			
	# Return this functions variables
	tupPreProcess = ithmcrundir,FID,id,sex,age,xgrid,ygrid,genes,\
	nogrids,subpop,fitvals,infection,Infected,subpopmigration,subpopemigration,agemort,egg_lmbdavals,egg_sigmavals,allelst,newmortperc,Mmature,Fmature,intgenesans
	return tupPreProcess
	
	#End::DoPreProcess()
	
# ---------------------------------------------------------------------------------------------------	 		
def DoUserInput(fileans):
	
	'''
	DoUserInput()
	This function reads in the user input and 
	stores the variables.
	'''
	
	# Open file for reading
	inputfile = open(fileans,'r')

	# Read lines from the file
	lines = inputfile.readlines()

	#Close the file
	inputfile.close()

	# Create an empty matrix to append to
	inputvariables = []

	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.split(',')
		inputvariables.append(thisline)
		
	# Delete lines
	del(lines)

	return inputvariables
	
	#End::DoUserInput()