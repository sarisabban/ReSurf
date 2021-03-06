#!/usr/bin/python3

# Author: Sari Sabban
# Email:  sari.sabban@gmail.com
# URL:    https://github.com/sarisabban
#
# Created By:   	Sari Sabban
# Created Date: 	20 February 2017
# Modified Date: 	4 May 2017

import re
import itertools
import numpy
import sys
from Bio.PDB import *

Design=sys.argv[1]					#Designed PDB file
start=int(sys.argv[2])					#Motif start residue
end=int(sys.argv[3])					#Motif end residue
RCSB=sys.argv[4]					#RCSB PDB file

							#_
p=PDBParser()						# |
structure_Design=p.get_structure('X',Design)		# |
structure_RCSB=p.get_structure('X',RCSB)		# |
model_Design=structure_Design[0]			# |Standard structure to setup biopython's DSSP to calculate SASA using Wilke constants
model_RCSB=structure_RCSB[0]				# |
dssp_Design=DSSP(model_Design,Design,acc_array='Wilke')	# |
dssp_RCSB=DSSP(model_RCSB,RCSB,acc_array='Wilke')	#_|



list_Design=list()
for x in dssp_Design:					#Loop to isolate SASA for each amino acid for Design structure
	if x[1]=='A':sasa=129*(x[3])
	elif x[1]=='V':sasa=174*(x[3])
	elif x[1]=='I':sasa=197*(x[3])
	elif x[1]=='L':sasa=201*(x[3])
	elif x[1]=='M':sasa=224*(x[3])
	elif x[1]=='P':sasa=159*(x[3])
	elif x[1]=='Y':sasa=263*(x[3])
	elif x[1]=='F':sasa=240*(x[3])
	elif x[1]=='W':sasa=285*(x[3])
	elif x[1]=='R':sasa=274*(x[3])
	elif x[1]=='N':sasa=195*(x[3])
	elif x[1]=='C':sasa=167*(x[3])
	elif x[1]=='Q':sasa=225*(x[3])
	elif x[1]=='E':sasa=223*(x[3])
	elif x[1]=='G':sasa=104*(x[3])
	elif x[1]=='H':sasa=224*(x[3])
	elif x[1]=='K':sasa=236*(x[3])
	elif x[1]=='S':sasa=155*(x[3])
	elif x[1]=='T':sasa=172*(x[3])
	elif x[1]=='D':sasa=193*(x[3])

	if (x[2]=='G' or x[2]=='H' or x[2]=='I' or x[2]=='B' or x[2]=='E') and sasa>=60: #Helix (GHI) and Sheet (BE)
		list_Design.append((x[0],x[1]))
	elif (x[2]=='-' or x[2]=='T' or x[2]=='S') and sasa>=40: #Loop (-TS)
		list_Design.append((x[0],x[1]))



list_RCSB=list()
for y in dssp_RCSB:
	list_RCSB.append((y[0],y[1]))



print('NATAA\nstart')

for x in list_Design:				#Double loop
	if start <= x[0] <= end:		#Do not include the morif sequence in the resurfacing code
		continue
	else:
		for y in list_RCSB:
			if x[0]==y[0]:			#If the two lists have the same positions
				if y[1] == 'X':
					print(x[0],'A','PIKAA','M')	#Some bacteria has MES instead of M, which is represented as X in DSSP, therefore replace the X with M when found
				else:
					print(x[0],'A','PIKAA',y[1])	#for everything else, print the RCSB amino acid
