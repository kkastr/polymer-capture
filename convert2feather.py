import feather as ft
import pandas as pd 

import numpy as np 

n = [10,20,50,100,200]

for i in n:
	infile = 'tmpr' + str(i) +'.csv'
	df = pd.read_csv(infile)
	outfile = 'r' + str(i) + '.feather'	
	ft.write_dataframe(df,outfile)



