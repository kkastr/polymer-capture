import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
mpl.rcParams.update({'font.size': 14})
from scipy import optimize

#mpl.rcParams.update({'font.size': 14})



num = [10,20,50,100,200]
stdevb = []
mtaub = []
stdev = []
mtau = []
rgp = []
mrgt = []
mrg = []
#df = pd.read_csv('metric-standard-default.csv',sep=",",header=None,names = ["N,tau,rgxtrans,rgytrans,rgztrans,rgxycorrtrans,rgxequil,rgyequil,rgzequil,rgxycorrequil,tfirstthread,tthread,tlastthread"])
df = pd.read_csv('metric-standard-capture.csv')
#df = pd.read_csv('summary-small-pore.csv')
print df
#print df

for key, grp in df.groupby('N'):
	# print grp['tau'].std()
	# print grp['tau'].mean()
	print len(grp['tau']), key, grp['tau'].std(), grp['tau'].mean()
	stdev.append(grp['tau'].std())
	mtau.append(grp['tau'].mean())
	#mrgt.append( np.sqrt(grp['rgxtrans'].mean()**2 + grp['rgytrans'].mean()**2 + grp['rgztrans'].mean()**2)  )


# for key, grp in df_base.groupby('N'):
# 	# print grp['tau'].std()
# 	# print grp['tau'].mean()
# 	print len(grp['tau']), key, grp['tau'].std(), grp['tau'].mean(axis=1), grp['tau'].std()/grp['tau'].mean(axis=1)
# 	stdevb.append(grp['tau'].std())
# 	mtaub.append(grp['tau'].mean())




stdev = np.array(stdev)
mtau = np.array(mtau)
mrgt = np.array(mrgt)
# stdevb = np.array(stdevb)
# mtaub = np.array(mtaub)

cov = stdev/mtau
# covb = stdevb/mtaub

df['rgr2'] = df['rgxtrans']**2 + df['rgytrans']**2 + df['rgztrans']**2

df['rgr'] = np.sqrt(df['rgr2'])


for key, grp in df.groupby('N'):
	mrg.append(grp['rgr'].mean())



logN = np.log(num)
logtau = np.log(mtau)


fitfunc= lambda f, x: f[1] + f[2] * x
errfunc = lambda f, x, y: (y - fitfunc(f, x))

out,success = optimize.leastsq(errfunc,[0,0,0],args=(logN, logtau),maxfev=3000)

print out[2],out[1]


yerr = stdev/1000.0
print stdev,yerr

x = np.arange(10,200)
y = np.exp(out[1])*x**out[2]

#plt.errorbar(num,mtau,yerr=yerr,fmt='o')
# plt.loglog(num,mtau,'bo')
# plt.loglog(x,y,'r--')
#plt.title(r'$\tau \sim 1.186$')
# df['rgxtrans'].hist(by=df['N'])
# df['rgytrans'].hist(by=df['N'])
# df['rgztrans'].hist(by=df['N'])

#plt.show()


# df[df['N']==100].plot(x='rgxtrans',y='tau',s=100,kind='scatter')
# plt.title("x")
# plt.show()
# df[df['N']==100].plot(x='rgytrans',y='tau',s=100,kind='scatter')
# plt.title("y")
# plt.show()
# df[df['N']==100].plot(x='rgztrans',y='tau',s=100,kind='scatter')
# plt.title("z")
# plt.show()

df['tau'].hist(by=df['N'],bins=100)
# plt.show()
for key, grp in df.groupby('N'):
	
	grp.plot(x='contact',y='tau',kind='scatter')
	plt.show()

#plots for tau, stdev, cov, tau vs rgp (do this for like, N=150)

# plt.plot(num,cov,'ko',ms=10,label='Pre-Confined')
# plt.plot(num,covb,'bo',ms=10,label='Baseline')
# plt.xlabel("N")
# plt.ylabel("Coefficient of Variation")
# plt.xlim([np.min(num)-5,np.max(num)+5])
# plt.legend()
# plt.show()

# x = np.arange(10,80)
# x2 = np.arange(20,200)
# y = 0.3*x*0.5
# y2 = 0.004*x2**2
# plt.loglog(x,y,'r--',lw=3)
# plt.loglog(x2,y2,'r--',lw=3)
# plt.loglog(num,stdev,'ko',ms=10,label='Pre-Confined')
# plt.loglog(num,stdevb,'bo',ms=10,label='Baseline')
# plt.ylabel('Standard Deviation')
# plt.xlabel('N')
# plt.legend()
# plt.xlim([np.min(num)-5,np.max(num)+100])
# plt.show()

# logN = np.log(num)
# logtau = np.log(mtau)


# fitfunc= lambda f, x: f[1] + f[2] * x
# errfunc = lambda f, x, y: (y - fitfunc(f, x))

# out,success = optimize.leastsq(errfunc,[0,0,0],args=(logN, logtau),maxfev=3000)

# print out[1],out[2]

# x = np.arange(10,200)
# y = 0.15*x**out[2]
# plt.loglog(x,y,'r--',lw=3)
# plt.loglog(num,mtau,'ko',ms=10,label='Pre-Confined')
# plt.loglog(num,mtaub,'bo',ms=10,label='Baseline')
# plt.legend(loc='upper left')
# plt.xlabel('N')
# plt.ylabel('Mean Translocation Time')
# plt.xlim([np.min(num)-5,np.max(num)+100])
# plt.show()


#folding distributions, Ask Martin how to only get 1 N

# df_full['s_lt'].hist(by=df_full['N'])


# df[df['N']==100].plot(x='Rgr_lt',y='tau',s=100,kind='scatter')
# plt.xlabel('Radius of Gyration')
# plt.ylabel('Translocation Time')
# plt.show()

# df_full[df_full.N==200]['s_lt'].hist()
# plt.ylabel('Frequency')
# plt.xlabel('Monomer at last thread')
# plt.show()

# plt.show()