import argparse
import numpy as np
import pandas as pd
import random as rd
import sys

parser = argparse.ArgumentParser(description='C.I. first step python script, last update: 10/2024')
parser.add_argument('--i', help='Input File', required=True, type=str)
parser.add_argument('--o', help='Output File Prefix', required=True, type=str)
parser.add_argument('--t', help='Test type (xa or aa)', required=True, type=str)
parser.add_argument('--s', help='Number of Sterile individuals', type=int, required=True)
parser.add_argument('--f', help='Number of Fertile individuals', type=int, required=True)
parser.add_argument('--p', help='Penatrance (25,50,75,100)', type=str, required=True)
parser.add_argument('--bs', help='Background Sterility (As integer value of percent)', type=int, required=True)
parser.add_argument('--sd', help='Random seed generator', type=int, required=False)
args = parser.parse_args()

number_sterile = args.s
number_fertile = args.f
np.seterr(all='raise') #Catch-all for numpy error handling

##### If seeding of random module is desired, uncomment line below
rd.seed(args.sd)

### Chromosome window boundaries (zero-indexed)
## Adjust as needed
X_end = 545
Chr2_end = 1523
Chr3_end = 2579

### Load Array, Transpose (Individuals read in row-wise individual and converted to column-wise)
### Needed for current sibsam output format from sibsam_flies_2579windows_300_ind_Dec2021.txt
num_array = np.loadtxt(args.i, dtype=int)
num_array = num_array.transpose()
df = pd.DataFrame(num_array)
del(num_array)

### Get individuals
s_id = []
f_id = []

#Make a random incompatibility
what_focal = rd.randint(0,1)
if what_focal == 0:
	focal_1 = 0
	focal_2 = 2
else:
	focal_1 = 2
	focal_2 = 0

### Get BD cells
if args.t == 'aa':
	len_1 = list(range(X_end, Chr2_end))
	len_2 = list(range(Chr2_end, Chr3_end))
elif args.t == 'xa':
	len_1 = list(range(0, X_end))
	len_2 = list(range(X_end, Chr3_end))
else:
	sys.exit('Incorrect test type entered. Exiting...')

inc_1 = rd.choice(len_1)
inc_2 = rd.choice(len_2)

pen_dict = {'25':'3', '50':'2', '75':'1', '100':'0'} 
pen_num = int(pen_dict.get(args.p))

### Simulate replicates
i=0
while (len(s_id) < number_sterile) or (len(f_id) < number_fertile):
	if i == int(df.shape[1]):
		sys.exit("Exiting program, number of individuals exhausted...") #Exit if groups are not filled after individuals are exhausted
	site_1 = df.iat[inc_1, i]
	site_2 = df.iat[inc_2, i]
	if (site_1 == focal_1) and (site_2 == focal_2): #If focal window1/window2
		if rd.randint(0,3) >= pen_num: #25% penetrance
			if len(s_id) < number_sterile:
				s_id.append(i)
			else:
				pass
		else:
			if rd.randrange(0,100) < args.bs: #Background sterility
				if len(s_id) < number_sterile:
					s_id.append(i)
				else:
					pass
			else:
				if len(f_id) < number_fertile:
					f_id.append(i)
				else:
					pass
	else:
		if rd.randrange(0,100) < args.bs:
			if len(s_id) < number_sterile:
				s_id.append(i)
			else:
				pass
		else:
			if len(f_id) < number_fertile:
				f_id.append(i)
			else:
				pass
	i = i + 1

bd_1 = [] # S, W1F, W2F
bd_2 = [] # S, W1NF, W2F
bd_3 = [] # F, W1F, W2F
bd_4 = [] # F, W1NF, W2F
bd_5 = [] # S, W1F, W2NF
bd_6 = [] # S, W1NF, W2NF
bd_7 = [] # F, W1F, W2NF
bd_8 = [] # F, W1NF, W2NF

# Forward count
for window_1 in len_1:
	for window_2 in len_2:
		b1 = 0
		b2 = 0
		b3 = 0
		b4 = 0
		b5 = 0
		b6 = 0
		b7 = 0
		b8 = 0
		for index in s_id: #Sterile
			if df.iat[window_1, index] == focal_1:
				if df.iat[window_2, index] == focal_2:
					b1 += 1 # S, W1F, W2F
				else:
					b5 += 1 # S, W1F, W2NF
			else: 
				if df.iat[window_2, index] == focal_2:
					b2 += 1 # S, W1NF, W2F
				else:
					b6 += 1 # S, W1NF, W2NF
		for index in f_id:
			if df.iat[window_1, index] == focal_1:
				if df.iat[window_2, index] == focal_2:
					b3 += 1 # F, W1F, W2F
				else:
					b7 += 1 # F, W1F, W2NF
			else:
				if df.iat[window_2, index] == focal_2:
					b4 += 1 # F, W1NF, W2F
				else:
					b8 += 1 # F, W1NF, W2NF
		bd_1.append(b1)
		bd_2.append(b2)
		bd_3.append(b3)
		bd_4.append(b4)
		bd_5.append(b5)
		bd_6.append(b6)
		bd_7.append(b7)
		bd_8.append(b8)
# Reverse count
for window_1 in len_2:
	for window_2 in len_1:
		b1 = 0
		b2 = 0
		b3 = 0
		b4 = 0
		b5 = 0
		b6 = 0
		b7 = 0
		b8 = 0
		for index in s_id: #Sterile
			if df.iat[window_1, index] == focal_1:
				if df.iat[window_2, index] == focal_2:
					b1 += 1 # S, W1F, W2F
				else:
					b5 += 1 # S, W1F, W2NF
			else: 
				if df.iat[window_2, index] == focal_2:
					b2 += 1 # S, W1NF, W2F
				else:
					b6 += 1 # S, W1NF, W2NF
		for index in f_id:
			if df.iat[window_1, index] == focal_1:
				if df.iat[window_2, index] == focal_2:
					b3 += 1 # F, W1F, W2F
				else:
					b7 += 1 # F, W1F, W2NF
			else:
				if df.iat[window_2, index] == focal_2:
					b4 += 1 # F, W1NF, W2NF
				else:
					b8 += 1 # F, W1NF, W2NF
		bd_1.append(b1)
		bd_2.append(b2)
		bd_3.append(b3)
		bd_4.append(b4)
		bd_5.append(b5)
		bd_6.append(b6)
		bd_7.append(b7)
		bd_8.append(b8)

assert len(bd_1) == len(bd_2) == len(bd_3) == len(bd_4) == len(bd_5) == len(bd_6) == len(bd_7) == len(bd_8) # Sanity

### Print to output
loci_df = pd.DataFrame([inc_1, focal_1, inc_2, focal_2])
cell_df = pd.DataFrame([bd_1, bd_2, bd_3, bd_4, bd_5, bd_6, bd_7, bd_8])
cell_df = cell_df.transpose()
cell_df.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
out_name = args.o + '_' + args.t + '_bd_cells_for_CI.csv'
out_name_2 = args.o + '_' + args.t + '_focal_loci.txt'
cell_df.to_csv(out_name, header=True, index=False)
loci_df.to_csv(out_name_2, header=False, index=False)
