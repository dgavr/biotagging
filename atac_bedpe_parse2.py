#!/usr/bin/env python
####### ATAC bedpe parsing 
<<<<<<< HEAD
######  Daria Gavriouchkina
=======
###### (c) Daria Gavriouchkina
>>>>>>> 186ae74e6324f3d8c86db0510f18373d8d731ea3
###### Sauka-Spengler lab
import sys
from collections import defaultdict
import math


def make_dict_with_chrom_lengths(chromFilename):
	"""
	"""
	#print 'getting chromosome lengths ...'
	chrom_dict = defaultdict(int)
	chromFile = open(chromFilename,'r')
	for chrom_line in chromFile:
		l=chrom_line.rstrip().split('\t')
		chr_chr = str(l[0])
		chr_max_pos = int(l[1])
		chrom_dict[chr_chr]= int(chr_max_pos)
	##print chrom_dict['Zv8_NA26']
	#print 'finished doing that now ...'
	return chrom_dict

def go_thru_bedpe(bedfile, chrom_dict):
	ct_tot=0
	ct_disc=0
	ct_kept=0
	header='.'.join(bedfile.split('.')[:-1])

	pe_bed_filename=header+'.pebed'
	pe_bed=open(pe_bed_filename, 'w')

	for line in open(bedfile):
		ct_tot+=1
		l=line.rstrip().split('\t')
		if not (l[0:3] == ['.','-1','-1']) and not (l[3:6] == ['.','-1','-1']) and not (([l[0]]+[l[3]]) == ['chrM', 'chrM'])  and not ('Zv9_' in line) :
			ct_kept+=1
			(chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2)= (str(l[0]),int(l[1]),int(l[2]), str(l[3]), int(l[4]),int(l[5]),str(l[6]),str(l[7]),str(l[8]),str(l[9]))

			#print chrom1, start1,end1, chrom2, start2, end2, name, score, strand1, stand2
			new_start1= start1+4 -50
			new_end1 = new_start1 + 50
			new_start2 = start2+4 -50
			new_end2 = new_start2 + 50
			if new_start1 < 0 :
				new_start1= 0
			if new_start2 < 0:
				new_start2 = 0
			if new_end1 > int(chrom_dict[chrom1]):
				new_end1 = int(chrom_dict[chrom1])
			if new_end2 > int(chrom_dict[chrom1]):
				new_end2 = int(chrom_dict[chrom1])

			pe_bed.write('\t'.join([chrom1, str(new_start1), str(new_end1), name, score, strand1])+'\n'+'\t'.join([chrom2, str(new_start2), str(new_end2), name, score, strand2])+'\n')

		else:
			ct_disc +=1
	print 'total number of read pairs =', ct_tot
	print 'number of reads discarded =', ct_disc
	print 'number of reads kept =', ct_kept

if __name__ == "__main__":
	chrom_dict=make_dict_with_chrom_lengths(sys.argv[1])
	go_thru_bedpe(sys.argv[2], chrom_dict)
