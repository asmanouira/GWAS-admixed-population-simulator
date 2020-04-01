import os
import pandas as pd
from textwrap import wrap




def convert_phased(population):

	"""
	Convert phased HapMap3 data to GWAsimulator format where the biallelic SNP are encoded by 0 or 1:
	For 23 chromosomes (22 autosomal chromosomes, X chromosome) 
	Chromosome X is denoted by 'chr23'

	Args:

		population: str: specify the HapMap3 population name, such as "ceu" or yri".

	Returns:

		Print output phased file for each chromosome.

	"""

	for ch in range(1,24):

		# Input files path
		file = population+'/hapmap3_r2_b36_fwd.consensus.qc.poly.chr'+str(ch)+'_'+population+'.phased'

		# Transform phased file into pandas Dataframe object
		chrom= pd.read_csv(file, delim_whitespace=True)

		# Delete extra columns
		chrom = chrom.drop(chrom.columns[0], axis = 1)
		chrom = chrom.drop(chrom.columns[0], axis = 1)

		# Merge all columns 
		chrom = chrom.astype(str).sum(axis=1)

		# Convert it to list
		chrom_list =[]
		for ind in chrom.index:
		        chrom_list.append( chrom[ind] )

		# Determine the major and the minor allele for each SNP
		alleles = [ list(set(i)) for i in chrom_list]
		
		# Remove useless charecters
		for i in range(len(alleles)):
			if "-" in alleles[i]: 
				alleles[i].remove("-")
			if 'n' in alleles[i]:
				alleles[i].remove('n')

		# Converting 
		phased = []
		for ic in range(len(chrom_list)):
		    
		    if len(alleles[ic]) == 2:
		        allele_1 = alleles[ic][0]
		        allele_2 = alleles[ic][1]

		        count_allele_1 = chrom_list[ic].count(allele_1)
		        count_allele_2 = chrom_list[ic].count(allele_2)

		        if count_allele_1 >= count_allele_2:
		            major_allele = allele_1
		            minor_allele = allele_2
		        elif count_allele_1 < count_allele_2:
		            major_allele = allele_2
		            minor_allele = allele_1
		         
		    if len(alleles[ic]) == 1:
		        major_allele = alleles[ic][0]
		        
		    SNP = ""  
		    for jc in range(len(chrom_list[ic])):
		        
		        if chrom_list[ic][jc] == major_allele:
		            SNP = SNP+'1'
		        elif chrom_list[ic][jc] == minor_allele:
		            SNP = SNP+'0'
		        elif chrom_list[ic][jc] == "-":
		        	SNP = SNP+'0'
		        elif chrom_list[ic][jc] == "a":
		        	SNP = SNP + '0'
		    phased.append(SNP)

		phased_wrapped = [wrap(i,2) for i in phased]
		phased_tr = map(list(map(None,*phased_wrapped))

	    
	    # final dataframe
		df = pd.DataFrame(phased_tr)

		# save data in text files
		# each file contains data for one chromosome 
		result_dir = population+'_encoded/'
		
		if not os.path.exists(result_dir):
			os.makedirs(result_dir)
   
		output= result_dir+'chr'+str(ch)+'_'+population+'.phased'
		df.to_csv(output, sep=' ', header = None, index = False)

def map(population):
	"""
	Genarates .map PLINK files
	Markers data: each row represents a SNP for 4 columns:
	- Chromosome number
	- Marker ID
	- Genetic distance
	- Physical position

	Args: 
		population: str: population name such as 'ceu'/ 'yri'
	"""
	for i in range(1,24):

		# Path to phased data used to generate genotype
		chrom_file =  population+'/hapmap3_r2_b36_fwd.consensus.qc.poly.chr'+str(i)+'_'+population+'.phased'
		chrom = pd.read_csv(chrom_file, delim_whitespace=True)
		# Obtain the number of chromosome
		chr_num = [i]*len(chrom)
		# The genetic distance will be set as 
		dist = [0]*len(chrom)
		# Marker IDs 
		col = chrom.columns
		map_file = pd.DataFrame(list(zip(chr_num, list(chrom[col[0]]), dist, list(chrom[col[1]]))), 
	               columns =['chr', 'ID', 'dist', 'pos']) 	
		 
		
		output= 'chr'+str(i)+'_'+population+'.map'
		# print map file
		map_file.to_csv(output, sep='\t', header = None, index= False)


def updateID(N_samples):

	"""
	Print text file wich is used as input for plink in order to rename individual IDs.
	Args:
		N_samples: the number of samples for the generated data.
	Returns:
		text file containes 4 fields: old_FID, IID, new_FID, IID).
	"""
	new_FID = range(N_samples+1,(N_samples*2)+1)
	old_FID = range(1,N_samples+1)
	IID = [1]*N_samples

	df = pd.DataFrame(list(zip(old_FID, IID, new_FID, IID)), columns =['1', '2', '3', '4']) 
	df.to_csv('ids.txt', sep=' ', header = None, index = False)
