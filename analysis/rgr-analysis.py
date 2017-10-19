import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
mpl.rcParams.update({'font.size': 14})
from scipy import optimize


rg10_20 = np.loadtxt('tmprgr50-20.csv')
rg10_40 = np.loadtxt('tmprgr50-40.csv')
rg10_60 = np.loadtxt('tmprgr50-60.csv')
rg10_110 = np.loadtxt('tmprgr50-110.csv')



rgr10_20 = rg10_20[:,0]
rgr10_40 = rg10_40[:,0]
rgr10_60 = rg10_60[:,0]
rgr10_110 = rg10_110[:,0]

rgp10_20 = (rg10_20[:,1] + rg10_20[:,2]) * 0.5
rgp10_40 = (rg10_40[:,1] + rg10_40[:,2]) * 0.5
rgp10_60 = (rg10_60[:,1] + rg10_60[:,2]) * 0.5
rgp10_110 = (rg10_110[:,1] + rg10_110[:,2]) * 0.5


transportdist = [20,40,60,110]
meanrgr = [np.mean(rgr10_20),np.mean(rgr10_40),np.mean(rgr10_60),np.mean(rgr10_110)]
meanrgp = [np.mean(rgp10_20),np.mean(rgp10_40),np.mean(rgp10_60),np.mean(rgr10_110)]
stdrgr = [np.std(rgr10_20),np.std(rgr10_40),np.std(rgr10_60),np.std(rgr10_110)]
stdrgp = [np.std(rgp10_20),np.std(rgp10_40),np.std(rgp10_60),np.std(rgr10_110)]


plt.plot(transportdist,meanrgr,'bo')
plt.show()
plt.plot(transportdist,meanrgp,'ro')
plt.show()


# plt.plot(transportdist,stdrgr,'bo')
# plt.show()
# plt.plot(transportdist,stdrgp,'ro')
# plt.show()


# rg10_20 = np.loadtxt('tmprgr50-20.csv')
# rg10_40 = np.loadtxt('tmprgr50-40.csv')
# rg10_60 = np.loadtxt('tmprgr50-60.csv')


# rgr10_20 = rg10_20[:,0]
# rgr10_40 = rg10_40[:,0]
# rgr10_60 = rg10_60[:,0]

# rgp10_20 = (rg10_20[:,1] + rg10_20[:,2]) * 0.5
# rgp10_40 = (rg10_40[:,1] + rg10_40[:,2]) * 0.5
# rgp10_60 = (rg10_60[:,1] + rg10_60[:,2]) * 0.5


# transportdist = [20,40,60]
# meanrgr = [np.mean(rgr10_20),np.mean(rgr10_40),np.mean(rgr10_60)]
# meanrgp = [np.mean(rgp10_20),np.mean(rgp10_40),np.mean(rgp10_60)]
# stdrgr = [np.std(rgr10_20),np.std(rgr10_40),np.std(rgr10_60)]
# stdrgp = [np.std(rgp10_20),np.std(rgp10_40),np.std(rgp10_60)]


# plt.plot(transportdist,meanrgr,'bo')
# plt.show()
# plt.plot(transportdist,meanrgp,'ro')
# plt.show()


plt.show()