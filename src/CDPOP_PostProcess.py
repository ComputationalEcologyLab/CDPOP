# -------------------------------------------------------------------------------------------------
# CDPOP_PostProcess.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file post processing.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
except ImportError:
	raise ImportError, "Numpy required."
import pdb,random,os,sys,glob
# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

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
def DoDaGeneticDistance(loci,nogrids,alleles,\
ithmcrundir,nthfile,logfHndl):
	'''
	DoDaGeneticDistance()
	This function outputs the genetic distance matrix 
	following Nei's algorithm.
	'''	
	
	print('Nei Da is not functional yet for versions greater than 1.2. Email Erin.')
	sys.exit(-1)
	
	# Insert -1 at beginning of nthfile, to get grid0.csv
	nthfile.insert(0,-1)
	
	# Get the specified nthfile's grid to calculated matrix
	for i in xrange(len(nthfile)):
						
		# Open file for reading
		inputfile = open(ithmcrundir+'grid'+str(nthfile[i]+1)+'.csv','r')
		
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		x = []
		
		# Split up each line in file and append to empty matrix, x
		for l in lines:
			thisline = l.split(',')
			x.append(thisline)
				
		# Store genetic information: genes[individual], but need them as float values
		genes = []
		tempgenes = []
		for k in xrange(len(x)-1):
			# Get list from read in file
			tempgenes.append(x[k+1][8:int(8+sum(alleles))])
			# Create spot in genes
			genes.append([])
			for j in xrange(sum(alleles)):
				# Make each list spot an integer
				genes[k].append(float(tempgenes[k][j]))
		
		# Create a matrix of zeros to be filled
		gendmatrix = np.zeros((nogrids,nogrids),float)
		
		# Loop through each individual k
		for k in xrange(nogrids):
			# Compare individual k to every other inidividual j
			for j in xrange(nogrids):
				# Create a tempvariable to be written over for each comparison
				tempsqrt=[]
				# Loop through each locus
				for locus in xrange(loci):						
					# Loop through each allele value
					for alle in xrange(alleles[locus]):
						# Find the allele frequencies between k and j checking the 4 conditions
						if genes[k][alle+sum(alleles[0:locus])]==2.0:
							if genes[j][alle+sum(alleles[0:locus])]==2.0:
								tempsqrt.append(np.sqrt(1.0))
							elif genes[j][alle+sum(alleles[0:locus])]==1.0:
								tempsqrt.append(np.sqrt(1.0*0.5))
						elif genes[k][alle+sum(alleles[0:locus])]==1.0:
							if genes[j][alle+sum(alleles[0:locus])]==1.0:
								tempsqrt.append(np.sqrt(0.5*0.5))
							elif genes[j][alle+sum(alleles[0:locus])]==2.0:
								tempsqrt.append(np.sqrt(0.5*1.0))
				# Write the Da value to gendmatrix
				gendmatrix[k][j] = 1-float(sum(tempsqrt))/(loci)
		
		# Transpose cdmatrix to get vertical display (just being picky)
		gendmatrix = np.transpose(gendmatrix)
		
		# Create file to write matrix to
		outputfile = open(ithmcrundir+'Gdmatrix'+str(nthfile[i]+1)+'.csv','w')
		
		# Sequence each row in the matrix
		for seqrow in gendmatrix:
		
			# Grab each element in each row and write element to outputfile
			for ele in xrange(len(seqrow)):
				outputfile.write(str(seqrow[ele]))
				# Add comma
				outputfile.write(',')
			
			# Return line
			outputfile.write('\n')
		
		# Close file
		outputfile.close()
		
		# Logging message
		stringout = 'Gdmatrix'+str(nthfile[i])+'.csv has been created'
		logMsg(logfHndl,stringout)
	
	# Delete the -1
	del(nthfile[0])
	
	# End::DoDaGeneticDistance()

# ---------------------------------------------------------------------------------------------------	
def DoDpsGeneticDistance(loci,nogrids,alleles,\
ithmcrundir,nthfile,logfHndl):
	'''
	DoDpsGeneticDistance()
	This function outputs the genetic distance matrix 
	following proportion of shared alleles algorithm.
	'''	
	# List the grid files in directory
	csvfileList = glob.glob(ithmcrundir+'grid*.csv')
		
	# Get the specified nthfile's grid to calculated matrix
	for i in xrange(len(csvfileList)):
						
		# Open file for reading
		inputfile = open(csvfileList[i],'r')
			
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		x = []
		
		# Split up each line in file and append to empty matrix, x
		for l in lines:
			thisline = l.split(',')
			x.append(thisline)
				
		# Store genetic information: genes[individual], but need them as float values
		genes = []
		tempgenes = []
		for k in xrange(len(x)-1):
			# Get list from read in file
			tempgenes.append(x[k+1][8:int(8+sum(alleles))])
			# Create spot in genes
			genes.append([])
			for j in xrange(sum(alleles)):
				# Make each list spot an integer
				genes[k].append(tempgenes[k][j].strip('\n'))
		
		# Create a matrix of zeros to be filled
		gendmatrix = []
		
		# Loop through each individual k
		tempcount = 0		
		for k in xrange(nogrids):
			
			# Break if NA
			if genes[k][0] == 'NA':
				continue
			else:
				gendmatrix.append([])
				tempcount = tempcount+1
				# Compare individual k to every other inidividual j
				for j in xrange(nogrids):
					
					# Break if NA
					if genes[j][0] == 'NA':
						continue
					else:
						# Create a tempvariable to be written over for each comparison
						tempmin=[]
						# Loop through each allele value
						for alle in xrange(sum(alleles)):
							# Find the shared alleles between k and j checking the 4 conditions
							if genes[k][alle]=='2':
								if genes[j][alle]=='2':
									tempmin.append(2)
								elif genes[j][alle]=='1':
									tempmin.append(1)
							elif genes[k][alle]=='1':
								if genes[j][alle]=='1':
									tempmin.append(1)
								elif genes[j][alle]=='2':
									tempmin.append(1)
						# Write the Dps value to gendmatrix
						gendmatrix[tempcount-1].append(1-float(sum(tempmin))/(2*loci))
		
		# Strip directory/filename of grid and add 'Gdmatrix.csv'
		gdpathname = csvfileList[i].replace('grid','Gdmatrix')
	
		# Create file to write matrix to
		outputfile = open(gdpathname,'w')
		
		# Sequence each row in the matrix
		for seqrow in gendmatrix:
		
			# Grab each element in each row and write element to outputfile
			for ele in xrange(len(seqrow)):
				outputfile.write(str(seqrow[ele]))
				# Add comma
				outputfile.write(',')
			
			# Return line
			outputfile.write('\n')
		
		# Close file
		outputfile.close()
		
		# Logging message
		stringout = gdpathname+' has been created'
		logMsg(logfHndl,stringout)
		print('The genetic distance matrix '+gdpathname+' has been created.')
		
	# End::DoDpsGeneticDistance()
		
# ---------------------------------------------------------------------------------------------------	
def DoBrayCurtisGeneticDistance(loci,nogrids,alleles,\
ithmcrundir,nthfile,logfHndl):
	'''
	DoBrayCurtisGeneticDistance()
	This function outputs the genetic distance matrix 
	following the Bray-Curtis algorithm.
	'''	
	print('Bray-Curtis is not functional yet for versions greater than 1.2. Email Erin.')
	sys.exit(-1)
	A = loci*2
	B = loci*2
	
	# Insert -1 at beginning of nthfile, to get grid0.csv
	nthfile.insert(0,-1)
	
	# Get the specified nthfile's grid to calculated matrix
	for i in xrange(len(nthfile)):
						
		# Open file for reading
		inputfile = open(ithmcrundir+'grid'+str(nthfile[i]+1)+'.csv','r')
			
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		x = []
		
		# Split up each line in file and append to empty matrix, x
		for l in lines:
			thisline = l.split(',')
			x.append(thisline)
				
		# Store genetic information: genes[individual], but need them as float values
		genes = []
		tempgenes = []
		for k in xrange(len(x)-1):
			# Get list from read in file
			tempgenes.append(x[k+1][8:int(8+sum(alleles))])
			# Create spot in genes
			genes.append([])
			for j in xrange(sum(alleles)):
				# Make each list spot an integer
				genes[k].append(float(tempgenes[k][j]))	
		
		# Create a matrix of zeros to be filled
		gendmatrix = np.zeros((nogrids,nogrids),float)

		# Loop through the rows
		for igrids in range(nogrids):
			
			# Loop through the columns
			for j in range(nogrids):
				
				# Create a tempvariable to be written over
				tempW=[]
				
				# Loop through each allele value
				for k in range(sum(alleles)):
					
					# Find the minimum value between the i and jth entry
					tempW.append(min(int(genes[igrids][k]),int(genes[j][k])))
					
				# Write the Curtis-Bray value to gendmatrix
				gendmatrix[igrids][j] = 1-(float((2*sum(tempW)))/(A+B))

		# Create file to write matrix to
		outputfile = open(ithmcrundir+'Gdmatrix'+str(nthfile[i]+1)+'.csv','w')
		
		# Sequence each row in the matrix
		for seqrow in gendmatrix:

			# Grab each element in each row and write element to outputfile
			for ele in xrange(len(seqrow)):
				outputfile.write(str(seqrow[ele]))
				# Add comma
				outputfile.write(',')
			
			# Return line
			outputfile.write('\n')

		# Close file
		outputfile.close()

		# Logging message
		stringout = 'Gdmatrix'+str(nthfile[i])+'.csv has been created'
		logMsg(logfHndl,stringout)
	
	# Delete the -1
	del(nthfile[0])
	
	# End::DoBrayCurtisGeneticDistance()
	
# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_cdpop(ithmcrundir,gen,loci,alleles,nogrids,subpopnew,xgridnew,\
ygridnew,idnew,sexnew,agenew,genesnew,logfHndl,infection,AllDispDistCD):
	'''
	DoGridOut_cdpop()
	Output grid.csv in cdpopformat	
	'''	
	
	# Create file to write info to
	outputfile = open(ithmcrundir+'grid'+str(gen+1)+'.csv','w')
	
	# Write out the titles - Add Titles from xypoints
	title = ['Subpopulation','XCOORD','YCOORD','ID','sex','age','infection','DisperseCDist']
	
	# Write out the title from xy points
	for i in xrange(len(title)):
		# Write out FID number
		outputfile.write(title[i]+',')
			
	# Write out the loci title info
	# Loop for loci length
	for i in xrange(loci-1):
		# Loop for allele length
		for j in xrange(alleles[i]):
			outputfile.write('L'+str(i)+'A'+str(j)+',')
	# To get a return character on the end of the title
	for i in xrange(alleles[loci-1]-1):
		outputfile.write('L'+str(loci-1)+'A'+str(i)+',')
	outputfile.write('L'+str(loci-1)+'A'+str(alleles[loci-1]-1)+'\n')

	# Loop through each grid spot and output
	for i in xrange(nogrids):
		
		outputfile.write(subpopnew[i]+',')
		outputfile.write(str(float(xgridnew[i]))+',')
		outputfile.write(str(float(ygridnew[i]))+',')
		outputfile.write(idnew[i]+',')
		outputfile.write(sexnew[i]+',')
		outputfile.write(str(agenew[i])+',')
		outputfile.write(str(infection[i])+',')
		outputfile.write(str(AllDispDistCD[i])+',')
		# Write out gene info
		for jk in xrange(loci-1):
			for kl in xrange(alleles[jk]):
				outputfile.write(str(genesnew[i][jk][kl])+',')
		# TO get return character on end
		for jk in xrange(alleles[loci-1]-1):
			outputfile.write(str(genesnew[i][loci-1][jk])+',')
		outputfile.write(str(genesnew[i][loci-1][alleles[loci-1]-1])+'\n')
										
	# Logging message
	stringout = 'The file grid'+str(gen+1)+'.csv has been created'
	logMsg(logfHndl,stringout)		
	
	# Close file
	outputfile.close()
	# End::DoGridOut_cdpop()
	
# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_general(loci,alleles,ithmcrundir,logfHndl):
	'''
	DoGridOut_general()
	Output grid.csv in general genotype format	
	'''	
		
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in xrange(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each grid
	for igrid in xrange(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		infection_cdpop = []
		subpop_cdpop = []
		for ispot in xrange(nogrids):
			subpop_cdpop.append(x[ispot+1][0])			
			x_cdpop.append(float(x[ispot+1][1]))
			y_cdpop.append(float(x[ispot+1][2]))
			id_cdpop.append(x[ispot+1][3])
			sex_cdpop.append(x[ispot+1][4])
			age_cdpop.append(x[ispot+1][5])
			infection_cdpop.append(x[ispot+1][6])
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in xrange(nogrids):
			genes_cdpop.append([])
			for jspot in xrange(loci):
				genes_cdpop[ispot].append(x[ispot+1][int(8+sum(alleles[0:jspot])):int(8+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(nogrids):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in xrange(loci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in xrange(alleles[ithloci]):
				
					# Check if allele spot is 1
					if genes_cdpop[ithind][ithloci][ithallele] == '1':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif genes_cdpop[ithind][ithloci][ithallele] == '2':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if NA in that spot
					elif genes_cdpop[ithind][ithloci][ithallele] == 'NA':
						
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append('NA')
						GenFormgenes[ithind][ithloci].append('NA')
						
					# Error check
					elif genes_cdpop[ithind][ithloci][ithallele] != '0':
						print('Something wrong in gene general format. Email Erin.')
						sys.exit(-1)					
			
		# Create file to write matrix to
		outputfilename = filename.split('grid')
		outputfile = open(outputfilename[0]+'/generalgrid'+outputfilename[1],'w')
				
		# Write out the titles that match general grid format
		title = ['Subpopulation','X','Y','ID','sex','age','infection']
			
		# Write out the title
		for ititle in range(len(title)):
			outputfile.write(title[ititle])
			outputfile.write(',')
			
		# Write out the loci title 
		for i in range(loci-1):
			outputfile.write('Locus'+str(i+1)+'a')
			outputfile.write(',')
			outputfile.write('Locus'+str(i+1)+'b')
			outputfile.write(',')
		# To get a return character on the end of the title
		outputfile.write('Locus'+str(loci-1+1)+'a')
		outputfile.write(',')
		outputfile.write('Locus'+str(loci-1+1)+'b')
		outputfile.write('\n')	
					
		# Loop through each grid spot and output
		for i in xrange(nogrids):
			
			outputfile.write(subpop_cdpop[i]+',')
			outputfile.write(str(float(x_cdpop[i]))+',')
			outputfile.write(str(float(y_cdpop[i]))+',')
			outputfile.write(str(id_cdpop[i])+',')
			outputfile.write(str(sex_cdpop[i])+',')
			outputfile.write(str(age_cdpop[i])+',')
			outputfile.write(str(infection_cdpop[i])+',')
			
			# Loop through each locus
			for ithloci in xrange(loci-1):
			
				# Loop through each allele spot at that locus
				for ithallele in xrange(2):
				
					outputfile.write(str(GenFormgenes[i][ithloci][ithallele])+',')
				
			# Return charater on end
			outputfile.write(str(GenFormgenes[i][loci-1][0])+',')
			outputfile.write(str(GenFormgenes[i][loci-1][1])+'\n')
											
		# Logging message
		stringout = 'The file grid'+outputfilename[0]+'/general'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	
	print('General grid format file conversion complete.')
	# End::DoGridOut_general()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_genalex(loci,alleles,ithmcrundir,logfHndl,subgridtotal):
	'''
	DoGridOut_genalex()
	Output grid.csv in genalex genotype format	
	'''	
	subpopno = len(subgridtotal)	
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in xrange(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each grid
	for igrid in xrange(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		infection_cdpop = []
		subpop_cdpop = []
		for ispot in xrange(nogrids):
			subpop_cdpop.append(x[ispot+1][0])			
			x_cdpop.append(float(x[ispot+1][1]))
			y_cdpop.append(float(x[ispot+1][2]))
			id_cdpop.append(x[ispot+1][3])
			sex_cdpop.append(x[ispot+1][4])
			age_cdpop.append(x[ispot+1][5])
			infection_cdpop.append(x[ispot+1][6])
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in xrange(nogrids):
			genes_cdpop.append([])
			for jspot in xrange(loci):
				genes_cdpop[ispot].append(x[ispot+1][int(8+sum(alleles[0:jspot])):int(8+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(nogrids):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in xrange(loci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in xrange(alleles[ithloci]):
				
					# Check if allele spot is 1
					if genes_cdpop[ithind][ithloci][ithallele] == '1':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif genes_cdpop[ithind][ithloci][ithallele] == '2':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if NA in that spot
					elif genes_cdpop[ithind][ithloci][ithallele] == 'NA':
						
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append('NA')
						GenFormgenes[ithind][ithloci].append('NA')
						
					# Error check
					elif genes_cdpop[ithind][ithloci][ithallele] != '0':
						print('Something wrong in gene genalex format. Email Erin.')
						sys.exit(-1)					
				
		# Create file to write matrix to
		outputfilename = filename.split('grid')
		outputfile = open(outputfilename[0]+'/genalexgrid'+outputfilename[1],'w')
			
		# Write out the first and second line of GENALEX format
		outputfile.write(str(loci)+',')
		outputfile.write(str(len(genes_cdpop))+',')
		outputfile.write(str(subpopno)+',')
		outputfile.write(str(len(genes_cdpop))+'\n')
		outputfile.write(filename+'\n')		
			
		# Write out the third line of GENALEX format
		outputfile.write('Individual ID,Population,')
		for i in range(loci):
			outputfile.write('locus'+str(i+1)+'a,')
			outputfile.write('locus'+str(i+1)+'b,')
		# The trailing white space and XY information
		outputfile.write(',X,Y\n')
							
		# Loop through each grid spot and output
		for i in xrange(nogrids):
			
			outputfile.write('indiv'+str(i)+',')
			outputfile.write(str(subpop_cdpop[i])+',')
			
			# Loop through each locus
			for ithloci in xrange(loci):
			
				# Loop through each allele spot at that locus
				for ithallele in xrange(2):
				
					outputfile.write(str(GenFormgenes[i][ithloci][ithallele])+',')
				
			# Write out trailing information - white space and x,y information
			outputfile.write(',')
			outputfile.write(str(x_cdpop[i]).strip('[').strip(']')+',')
			outputfile.write(str(y_cdpop[i]).strip('[').strip(']')+'\n')	
											
		# Logging message
		stringout = 'The file grid'+outputfilename[0]+'/genalex'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	
	print('GENALEX grid format file conversion complete.')
	# End::DoGridOut_genalex()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_structure(loci,alleles,ithmcrundir,logfHndl):
	'''
	DoGridOut_structure()
	Output grid.csv in structure genotype format	
	'''	
		
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in xrange(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each grid
	for igrid in xrange(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		infection_cdpop = []
		subpop_cdpop = []
		for ispot in xrange(nogrids):
			subpop_cdpop.append(x[ispot+1][0])			
			x_cdpop.append(float(x[ispot+1][1]))
			y_cdpop.append(float(x[ispot+1][2]))
			id_cdpop.append(x[ispot+1][3])
			sex_cdpop.append(x[ispot+1][4])
			age_cdpop.append(x[ispot+1][5])
			infection_cdpop.append(x[ispot+1][6])
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in xrange(nogrids):
			genes_cdpop.append([])
			for jspot in xrange(loci):
				genes_cdpop[ispot].append(x[ispot+1][int(8+sum(alleles[0:jspot])):int(8+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(nogrids):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in xrange(loci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in xrange(alleles[ithloci]):
				
					# Check if allele spot is 1
					if genes_cdpop[ithind][ithloci][ithallele] == '1':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif genes_cdpop[ithind][ithloci][ithallele] == '2':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if NA in that spot
					elif genes_cdpop[ithind][ithloci][ithallele] == 'NA':
						
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append('NA')
						GenFormgenes[ithind][ithloci].append('NA')
						
					# Error check
					elif genes_cdpop[ithind][ithloci][ithallele] != '0':
						print('Something wrong in gene genalex format. Email Erin.')
						sys.exit(-1)					
				
		# Create file to write matrix to
		outputfilename = filename.split('grid')
		outputfile = open(outputfilename[0]+'/structuregrid'+outputfilename[1].strip('.csv')+'.stru','w')	
		
		# Write out the first line of structure format
		for i in range(loci):
			outputfile.write('locus'+str(i+1)+' ')
		outputfile.write('\n\n')
							
		# Loop through each grid spot and output
		for i in xrange(nogrids):
		
			# Loop through each allele spot at that locus
			for ithallele in xrange(2):
				
				outputfile.write(str(id_cdpop[i])+' ')
				outputfile.write(str(subpop_cdpop[i])+' ')
								
				# Loop through each locus
				for ithloci in xrange(loci):
					
					outputfile.write(str(GenFormgenes[i][ithloci][ithallele])+' ')
					
				# Return character
				outputfile.write('\n\n')
															
		# Logging message
		stringout = 'The file grid'+outputfilename[0]+'/structure'+outputfilename[1]+'.stru has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	
	print('STRUCTURE grid format file conversion complete.')
	# End::DoGridOut_structure()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_genepop(loci,alleles,ithmcrundir,logfHndl,subgridtotal,subpop):
	'''
	DoGridOut_genalex()
	Output grid.csv in genalex genotype format	
	'''	
	subpopno = len(subgridtotal)	
	
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in xrange(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each grid
	for igrid in xrange(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		infection_cdpop = []
		subpop_cdpop = []
		for ispot in xrange(nogrids):
			subpop_cdpop.append(x[ispot+1][0])			
			x_cdpop.append(float(x[ispot+1][1]))
			y_cdpop.append(float(x[ispot+1][2]))
			id_cdpop.append(x[ispot+1][3])
			sex_cdpop.append(x[ispot+1][4])
			age_cdpop.append(x[ispot+1][5])
			infection_cdpop.append(x[ispot+1][6])
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in xrange(nogrids):
			genes_cdpop.append([])
			for jspot in xrange(loci):
				genes_cdpop[ispot].append(x[ispot+1][int(8+sum(alleles[0:jspot])):int(8+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(nogrids):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in xrange(loci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in xrange(alleles[ithloci]):
				
					# Check if allele spot is 1
					if genes_cdpop[ithind][ithloci][ithallele] == '1' or genes_cdpop[ithind][ithloci][ithallele] == '1\n':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif genes_cdpop[ithind][ithloci][ithallele] == '2' or genes_cdpop[ithind][ithloci][ithallele] == '2\n':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if NA in that spot
					elif genes_cdpop[ithind][ithloci][ithallele] == 'NA' or genes_cdpop[ithind][ithloci][ithallele] == 'NA\n':
						
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append('NA')
						GenFormgenes[ithind][ithloci].append('NA')
						
					# Error check
					elif genes_cdpop[ithind][ithloci][ithallele] != '0' and genes_cdpop[ithind][ithloci][ithallele] != '0\n':
						print('Something wrong in gene genepop format. Email Erin.')
						sys.exit(-1)					
				
		# Create file to write matrix to
		outputfilename = filename.split('grid')
		outputfile = open(outputfilename[0]+'/genepopgrid'+outputfilename[1].strip('.csv')+'.gen','w')
			
		# Write out the first and second line of GENEPOP format
		outputfile.write(outputfilename[0]+'grid'+outputfilename[1]+'\n')
		for i in xrange(loci):
			outputfile.write('LOCUS-'+str(i+1)+'\n')
					
		# Write out the genes of each individual by population
		for ipop in xrange(subpopno):
			outputfile.write('POP\n')
			
			# Locate the index for subpop for case that is unordered.
			popindex = np.where(np.asarray(subpop) == str(ipop+1))[0]

			# Loop through each grid spot and output
			for i in xrange(len(popindex)):
				# Individual ID
				outputfile.write(id_cdpop[popindex[i]]+', ')								
				# Loop through each locus
				for ithloci in xrange(loci-1):
					templociname = ''
					# Loop through each allele spot at that locus
					for ithallele in xrange(2):
						if GenFormgenes[popindex[i]][ithloci][ithallele] != 'NA':
							# Add 100
							templociname = templociname + str(GenFormgenes[popindex[i]][ithloci][ithallele]+100)
						else:
							templociname = templociname + str(GenFormgenes[popindex[i]][ithloci][ithallele])
					outputfile.write(templociname+' ')
				# Loop through last loci alleles for return character
				templociname = ''
				for ithallele in xrange(2):	
					if GenFormgenes[popindex[i]][loci-1][ithallele] != 'NA':
						# Add 100
						templociname = templociname + str(GenFormgenes[popindex[i]][loci-1][ithallele]+100)
					else:
						templociname = templociname + str(GenFormgenes[popindex[i]][loci-1][ithallele])
				outputfile.write(templociname+'\n')		
															
		# Logging message
		stringout = 'The file grid'+outputfilename[0]+'/genepop'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	print('GENEPOP grid format file conversion complete.')
	# End::DoGridOut_genepop()	
	
# ---------------------------------------------------------------------------------------------------	 
def DoOutput(nogrids,FID,OffDisperseIN,xgridcopy,ygridcopy,gen,\
id,sex,age,xgrid,ygrid,genes,nthfile,ithmcrundir,loci,alleles,subpop,\
logfHndl,gridformat,infection,Infected,cdinfect,opengrids,OffDispDistCD,geneswap):
	'''
	DoOutput()
	Generate .txt file of old+new+Immigration generations
	Input: ithmcrundir
	Output: ithmcrundir will have .csv files of x,y coord location values of
	cost distance dispersal of the old+new generation with gene info	
	'''	
	
	# ----------------------------------------------------------------------------------------
	# Order the grids back from the random processes - 0,1,2,...nogrids
	# ----------------------------------------------------------------------------------------
	
	# Storage lists for ordering id and no
	orderofgridid = []
	orderofgridno = []
	
	# checks
	if (len(FID) + len(OffDisperseIN) + len(opengrids)) != nogrids:
		pdb.set_trace()
	
	# Loop through all grid points
	for i in xrange(nogrids):
	
		# Loop through the FID values if any - existing adults
		for jFID in xrange(len(FID)):
			if int(FID[jFID])==i:
				orderofgridid.append('FID'+str(jFID))
				orderofgridno.append(jFID)
				
		# Loop through the dispersal values if any
		for jFG in xrange(len(OffDisperseIN)):
			if int(OffDisperseIN[jFG][1])==i:
				orderofgridid.append('FG'+str(jFG))
				orderofgridno.append(jFG)
				
		# Loop through the immigrant values if any
		for jopen in xrange(len(opengrids)):
			if int(opengrids[jopen])==i:
				orderofgridid.append('OPEN'+str(jopen))
				orderofgridno.append(jopen)
	
	# ------------------------------------------------------------------------------------------------------------------
	# Store the OffDisperseIN and Newimmigrants generation information: Update grid...
	# ------------------------------------------------------------------------------------------------------------------
	
	# Store new grid values
	FIDnew = []
	subpopnew = []
	xgridnew = []
	ygridnew = []
	idnew = []
	sexnew = []	
	agenew = []	
	genesnew = []
	infectionnew = []
	AllDispDistCD = []
	
	
	# Keep genes in within burn in period
	if gen < geneswap:
		genesnew = genes		
	
	# Extract information from OffDisperseIN
	if len(OffDisperseIN)!=0:
	
		# Temp variables storage
		xtempoff = []
		ytempoff = []
		offsex=[]
		offid=[]
		tempinf = []		
		
		# Extract grid x and y location, id,sex
		for dispot in xrange(len(OffDisperseIN)):
			xtempoff.append(xgridcopy[OffDisperseIN[dispot][1]])
			ytempoff.append(ygridcopy[OffDisperseIN[dispot][1]])
			offsex.append(OffDisperseIN[dispot][0][4])
			offid.append(OffDisperseIN[dispot][2])
			tempinf.append(OffDisperseIN[dispot][0][5])
		
# check orderofgridno	
	# Loop through each grid spot to find the right output index...
	for jgrid in xrange(nogrids):
			
		# Get the ordered grid location for the jth FID value
		if orderofgridid[jgrid][0:3] == 'FID':
				
			# Write out the old generation
			FIDnew.append(FID[orderofgridno[jgrid]])
			idnew.append(id[orderofgridno[jgrid]])
			sexnew.append(str(sex[orderofgridno[jgrid]]))
			# This adds a year to the age class
			if orderofgridno[jgrid] == 'NA':
				pdb.set_trace()
			if age[orderofgridno[jgrid]] == 'NA':
				pdb.set_trace()
			agenew.append(int(age[orderofgridno[jgrid]])+1)
			xgridnew.append(float(xgrid[orderofgridno[jgrid]]))
			ygridnew.append(float(ygrid[orderofgridno[jgrid]]))
			# Make sure this is the correct order writing to file!
			# Write out gene info
			if gen >= geneswap:
				genesnew.append(genes[orderofgridno[jgrid]])
			subpopnew.append(subpop[jgrid])
			infectionnew.append(infection[orderofgridno[jgrid]])
			AllDispDistCD.append('NoMove')
					
		# Get the ordered grid location for the jth FG value (the OffDisperseIN)
		elif orderofgridid[jgrid][0:2] == 'FG':
				
			# Write out the OffDisperseIN generation
			FIDnew.append(OffDisperseIN[orderofgridno[jgrid]][1])
			idnew.append(offid[orderofgridno[jgrid]])
			sexnew.append(str(int(offsex[orderofgridno[jgrid]])))
			# This makes the age of the offspring 0
			agenew.append(1)
			xgridnew.append(xtempoff[orderofgridno[jgrid]])
			ygridnew.append(ytempoff[orderofgridno[jgrid]])
			subpopnew.append(subpop[jgrid])
			infectionnew.append(int(tempinf[orderofgridno[jgrid]]))
			# Write out gene info
			if gen >= geneswap:
				genesnew.append([])
				for jloci in xrange(loci):
					genesnew[jgrid].append([])
					for jalleles in xrange(alleles[jloci]):
						genesnew[jgrid][jloci].append(OffDisperseIN[orderofgridno[jgrid]][0][6][jalleles+sum(alleles[0:jloci])])
			AllDispDistCD.append(OffDispDistCD[orderofgridno[jgrid]])			
		# Get the ordered grid location for the jth IMM value (the Newimmigrants)
		elif orderofgridid[jgrid][0:4] == 'OPEN':
					
			# Write out the Immigration generation
			FIDnew.append(opengrids[orderofgridno[jgrid]])
			idnew.append('OPEN')
			sexnew.append('NA')
			agenew.append('NA')
			xgridnew.append(float(xgridcopy[opengrids[orderofgridno[jgrid]]]))
			ygridnew.append(float(ygridcopy[opengrids[orderofgridno[jgrid]]]))
			subpopnew.append(subpop[jgrid])
			infectionnew.append('NA')
			# Write out gene info
			if gen >= geneswap:
				genesnew.append([])
				for jloci in xrange(loci):
					genesnew[jgrid].append([])
					for jalleles in xrange(alleles[jloci]):
						genesnew[jgrid][jloci].append('NA')
							
			AllDispDistCD.append('NA')	
	
	# ------------------------------------------------------------
	# Write out text file for generations specified by nthfile
	# ------------------------------------------------------------
	
	# Check if nthfile == generation.
	for inthfile in xrange(len(nthfile)):				
		if gen == nthfile[inthfile]:
			
			# Call DoGridOut_cdpop()
			DoGridOut_cdpop(ithmcrundir,gen,loci,alleles,nogrids,\
			subpopnew,xgridnew,ygridnew,idnew,sexnew,agenew,genesnew,\
			logfHndl,infectionnew,AllDispDistCD)
	
	# Sum infectednew to Infected
	temp = []
	for i in xrange(len(infection)):
		if infection[i] == 1:
			temp.append(1)
	Infected.append(len(temp))
	del(temp)
	
	# Return variables from this argument
	tupDoOut = FIDnew,idnew,sexnew,agenew,xgridnew,ygridnew,\
	genesnew,subpopnew,infectionnew,Infected
	return tupDoOut
	
	# End::DoOutput()
	
# ---------------------------------------------------------------------------------------------------	 
def DoPostProcess(ithmcrundir,nogrids,\
xgridcopy,ygridcopy,gendmatans,\
loci,alleles,looptime,Population,ToTFemales,ToTMales,\
BreedFemales,BreedMales,Migrants,Births,\
Deaths,Alleles,He,Ho,AllelesMutated,\
MateDistED,FDispDistED,MDispDistED,MateDistCD,FDispDistCD,\
MDispDistCD,nthfile,logfHndl,p1,p2,q1,q2,Infected,subpop,\
MateDistEDstd,FDispDistEDstd,MDispDistEDstd,MateDistCDstd,\
FDispDistCDstd,MDispDistCDstd,subpopmigration,FAvgMate,MAvgMate,\
FSDMate,MSDMate,DisperseDeaths,Open,CouldNotDisperse,\
Female_BreedEvents,gridformat,subpopemigration,females_nomate,\
subgridtotal,OffDeaths,Population_age,BreedFemales_age,subpopmatemort,SelectionDeaths,MateDistances):
	'''
	DoPostProcess()
	Create Distance Matrices - Geographic, Genetic, and Cost
	and output.csv file.
	'''	
	
	# -----------------------------
	# Matrix calculations
	# -----------------------------	
	# Create Genetic distance matrix with Bray-Curtis algorithm
	if gendmatans == 'braycurtis':	
		DoBrayCurtisGeneticDistance(loci,nogrids,alleles,ithmcrundir,nthfile,logfHndl)
	
	# Create Genetic distance matrix with proportion of shared alleles algorithm
	elif gendmatans == 'Dps':		
		DoDpsGeneticDistance(loci,nogrids,alleles,ithmcrundir,nthfile,logfHndl)
	
	# Create Genetic distance matrix with Nei's algorithm
	elif gendmatans == 'Da':	
		DoDaGeneticDistance(loci,nogrids,alleles,ithmcrundir,nthfile,logfHndl)
		
	# ------------------------
	# Grid format options
	# ------------------------
	# General format
	if gridformat == 'general':		
		DoGridOut_general(loci,alleles,ithmcrundir,logfHndl)
		
	# GENALEX format
	elif gridformat == 'genalex':
		DoGridOut_genalex(loci,alleles,ithmcrundir,logfHndl,subgridtotal)
		
	# STRUCTURE format
	elif gridformat == 'structure':
		DoGridOut_structure(loci,alleles,ithmcrundir,logfHndl)
	
	# GENEPOP format
	elif gridformat == 'genepop':
		DoGridOut_genepop(loci,alleles,ithmcrundir,logfHndl,subgridtotal,subpop)
			
	# -------------------------------------------------
	# output.csv
	# -------------------------------------------------
			
	# Calculate other metrics here
	# ----------------------------
	tempPop = np.asarray(Population,dtype='float')[:,0]
	growthPop = tempPop[1:]/tempPop[0:(len(tempPop)-1)]
	
	# Create time array
	time = np.arange(0,len(Population),1)
		
	# Get unique number of subpops
	nosubpops = len(np.unique(subpop))
	
	# Store Ouput in file
	# Create file to write info to
	outputfile = open(ithmcrundir+'output.csv','w')
	
	# Write out the titles
	# Add Titles from xypoints
	outputtitle = ['Year','Population','Population_Age','GrowthRate','ToTFemales','ToTMales',\
	'BreedFemales','BreedFemales_Age','BreedMales','BreedEvents_Females','Females_NoMate',\
	'Migrants','DisperseDeaths','Births','EggDeaths',\
	'AgeDeaths','Alleles','He','Ho','Mutations',\
	'MateDistED','MateDistEDstd','Female_DispDistED','Female_DispDistEDstd',\
	'Male_DispDistED','Male_DispDistEDstd','MateDistCD','MateDistCDstd','AllMateCDistances',\
	'Female_DispDistCD','Female_DispDistCDstd','Male_DispDistCD','Male_DispDistCDstd',\
	'p1','p2','q1','q2','Infected','SubpopImmigration','SubpopEmigration','SubpopNoMate','FemalesMeanMate',\
	'MalesMeanMate','FemalesSDMate','MalesSDMate','OpenLocations','CouldNotDisperse','MatureSelectionDeaths']
	
	# Write out the title
	for i in xrange(len(outputtitle)-1):
		outputfile.write(outputtitle[i])
		outputfile.write(',')
	# To get return character on the end
	outputfile.write(str(outputtitle[len(outputtitle)-1]))				
	outputfile.write('\n')
	
	# Write to file
	for i in xrange(len(time)):
		outputfile.write(str(time[i])+',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Population[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(Population_age[i])):
			outputfile.write(str(Population_age[i][j])+'|')
		outputfile.write(',')
		if i != len(time)-1:
			outputfile.write(str(growthPop[i])+',')
		else:
			outputfile.write('NA,')
		for j in xrange(nosubpops+1):
			outputfile.write(str(ToTFemales[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(ToTMales[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(BreedFemales[i][j])+'|')
		outputfile.write(',')
		for j in xrange(len(BreedFemales_age[i])):
			outputfile.write(str(BreedFemales_age[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(BreedMales[i][j])+'|')
		outputfile.write(',')		
		outputfile.write(str(Female_BreedEvents[i])+',')
		outputfile.write(str(len(females_nomate[i]))+',')
		outputfile.write(str(Migrants[i])+',')
		outputfile.write(str(DisperseDeaths[i])+',')
		outputfile.write(str(Births[i])+',')
		outputfile.write(str(OffDeaths[i])+',')
		for j in xrange(len(Deaths[i])):
			outputfile.write(str(Deaths[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Alleles[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(He[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Ho[i][j])+'|')
		outputfile.write(',')
		outputfile.write(str(AllelesMutated[i])+',')
		outputfile.write(str(MateDistED[i])+',')
		outputfile.write(str(MateDistEDstd[i])+',')
		outputfile.write(str(FDispDistED[i])+',')
		outputfile.write(str(FDispDistEDstd[i])+',')
		outputfile.write(str(MDispDistED[i])+',')
		outputfile.write(str(MDispDistEDstd[i])+',')
		outputfile.write(str(MateDistCD[i])+',')
		outputfile.write(str(MateDistCDstd[i])+',')
		for j in xrange(len(MateDistances[i])):
			outputfile.write(str(MateDistances[i][j])+'|')
		outputfile.write(',')
		outputfile.write(str(FDispDistCD[i])+',')
		outputfile.write(str(FDispDistCDstd[i])+',')
		outputfile.write(str(MDispDistCD[i])+',')
		outputfile.write(str(MDispDistCDstd[i])+',')
		outputfile.write(str(p1[i])+',')
		outputfile.write(str(p2[i])+',')
		outputfile.write(str(q1[i])+',')
		outputfile.write(str(q2[i])+',')
		outputfile.write(str(Infected[i])+',')
		for j in xrange(nosubpops):
			outputfile.write(str(subpopmigration[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops):
			outputfile.write(str(subpopemigration[i][j])+'|')
		outputfile.write(',')
		outputfile.write(str(subpopmatemort[i])+',')
		outputfile.write(str(FAvgMate[i])+',')
		outputfile.write(str(MAvgMate[i])+',')
		outputfile.write(str(FSDMate[i])+',')
		outputfile.write(str(MSDMate[i])+',')
		outputfile.write(str(Open[i])+',')
		outputfile.write(str(CouldNotDisperse[i])+',')
		outputfile.write(str(SelectionDeaths[i]))
		outputfile.write('\n')
				
	# Logging message
	stringout = 'The file outputfile.csv has been created'
	logMsg(logfHndl,stringout)	
	
	# Close file
	outputfile.close()