import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
mpl.rcParams.update({'font.size': 14})
from scipy import optimize
import feather as ft

narray = [10,20,50,100,200]

#rd1 = np.load('radialdist_10.npy')
#rd2 = np.load('radialdist_200.npy')
#rd3 = np.load('radialdist_50.npy')
#rd4 = np.load('radialdist_100.npy')
#rd5 = np.load('radialdist_200.npy')


# arr = rd2
# c = 0


# cs = []
# cutoff = 45
# for i in narray:

# 	filename = 'radialdist_' + np.str(i) + '.npy'

# 	fnload = np.load(filename)

# 	arr = fnload

# 	c = 0

# 	while (True):

# 		idx = np.argmax(arr < 0.6)
# 		slc = arr[0:idx]
# 		if np.argmax(slc>cutoff) !=0:
# 			c+=1
# 		tmp = arr[idx:len(arr)]

# 		while (tmp[0]<1):
# 			print len(tmp)
# 			tmp = tmp[1:len(tmp)]
# 		arr = tmp
# 		if np.argmax(arr<1) == 0 and arr[0] > 1:
# 			break

# 	cs.append(c)
diffarr = np.arange(10,120,10)
carray = []

lenarray = []
for i in diffarr:
	print i
	print i+40
	cs = []
	start_point = i
	cutoff = i + 40
	lens = []	
	for i in narray:

		filename = 'r' + np.str(i) + '.feather'

		df = ft.read_dataframe(filename)



		arr = df.values
		c = 0
		l=0

		while (True):

			idx = np.argmax(arr < 0.6)
			idx2 = np.argmax(arr< start_point)
			
			slc = arr[idx2:idx]
			#idx2 = np.argmax(slc<start_point)
			#slc2 = slc[idx2:-1]
			#print len(slc)
			if len(slc) == 0:
				break
			#print idx2
			if slc[np.argmax(slc>cutoff)] > cutoff:
				c+= 1
				#l+= len(slc[0:(np.argmax(slc>cutoff))])
			tmp = arr[idx:(len(arr)-1)]

			l+=1
			#print idx2, idx, len(slc), len(tmp), cutoff

			#print len(tmp)
			if len(tmp)==0:
				break 

			while (tmp[0]<1.1):
				#print len(tmp)
				#print tmp
				tmp = tmp[1:len(tmp)]
			arr = tmp
			if np.argmax(arr<1) == 0 and arr[0] > 1:
				break
		print l		
		cs.append(c)
		lens.append(l)
	

	carray.append(cs)
	lenarray.append(lens)

carray = np.array(carray)
lenarray = np.array(lenarray)





