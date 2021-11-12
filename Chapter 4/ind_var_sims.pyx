# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 10:32:20 2021

@author: Edith Invernizzi
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import wall_funcs as w


filepath = "path/to/folder/"

S = 100
stones = 3000
#sd_list = [1] # for debugging
sd_list = [0.5,1,1.5,2,5]


for sd in sd_list:
    
    TList = list()
    rugList = list()
    for s in range(S):
        
        if s == S-1:
            outT = w.templateOnly(s, 0, 18, R1 = 18, ind_var = 'on', var_sd = sd)
            outrug = w.rugatulus(s, 0, 18, R1 = 18, ind_var = 'on', var_sd = sd)
            TList.append(outT)
            rugList.append(outrug)
            
            imnmT = filepath + "images/template-indvar-sd"+ str(round(sd*100)) + "-s99-stones" + str(stones) + ".png"
            imgT = outT[6]
            plt.imshow(imgT, cmap = 'gray')
            plt.colorbar()
            plt.axis('off')
            plt.savefig(imnmT, dpi = 350)
            plt.clf()
            
            imnmr = filepath + "images/rugatulus-indvar-sd" + str(round(sd*100)) + "-s99-stones" + str(stones) + ".png"
            imgrug = outrug[6]
            plt.imshow(imgrug, cmap = 'gray')
            plt.colorbar()
            plt.axis('off')
            plt.savefig(imnmr, dpi = 350)
            plt.clf()
            
        else:
            TList.append(w.templateOnly(s, 0, 18, R1 = 18, ind_var = 'on', var_sd = sd))
            rugList.append(w.rugatulus(s, 0, 18, R1 = 18, ind_var = 'on', var_sd = sd))
    
    filenameT = filepath + "out_T_sd" + str(round(sd*100))
    fileT = open(filenameT, 'wb')
    pickle.dump(TList, fileT)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
    fileT.close()
        
    filenamer = filepath + "out_rug_sd" + str(round(sd*100))
    filer = open(filenamer, 'wb')
    pickle.dump(rugList, filer)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
    filer.close()
        
    print("sd is at " + str(sd))


# analyse results
from pandas import DataFrame

statListT = list()
statListr = list()
fsteps = 10
for sd in sd_list:
    
    loadpathT = filepath + "out_T_sd" + str(round(sd*100))
    fileT = open(loadpathT, 'rb')
    allResT = pickle.load(fileT)
    fileT.close()
        
    CoVList_T = [sum(f[4][-fsteps:])/fsteps for f in allResT]
    aveCoV_T = np.mean(CoVList_T)
    stdCoV_T = np.std(CoVList_T)
    RbarList_T = [sum(f[5][-fsteps:])/fsteps for f in allResT]
    aveRbar_T = np.mean(RbarList_T)
    stdRbar_T = np.std(RbarList_T)
    statListT.append(["T",sd,aveCoV_T,stdCoV_T,aveRbar_T,stdRbar_T])
    
        
    loadpathr = filepath + "out_rug_sd" + str(round(sd*100))
    filer = open(loadpathr, 'rb')
    allResr = pickle.load(filer)
    filer.close()
        
    CoVList_r = [sum(f[4][-fsteps:])/fsteps for f in allResr]
    aveCoV_r = np.mean(CoVList_r)
    stdCoV_r = np.std(CoVList_r)
    RbarList_r = [sum(f[5][-fsteps:])/fsteps for f in allResr]
    aveRbar_r = np.mean(RbarList_r)
    stdRbar_r = np.std(RbarList_r)
    statListr.append(["rug",sd,aveCoV_r,stdCoV_r,aveRbar_r,stdRbar_r])
    
    data = {'model': ['T']*S + ['rug']*S, 'SD': [sd]*2*S, 'CoV': CoVList_T + CoVList_r, 'Rbar': RbarList_T + RbarList_r}
    listDF = DataFrame.from_dict(data)
    svlistnm = filepath + "ind_var_list-sd" + str(round(sd*100)) + "-mean18-" + str(stones) + "stones.csv"
    listDF.to_csv(svlistnm, index = False)
    
statListAll = statListT + statListr
df = DataFrame(statListAll, columns = ["model",'sd','CoV_ave','CoV_SD','Rbar_ave', 'Rbar_SD'])
svname = filepath + 'ind_var_mean18_3000stones.csv'
df.to_csv(svname)


#################
# change in tau #
#################
filepath = "path/to/folder/"

stones = 3000
sd = 5
tau = 0.075
outrug3000 = w.rugatulus(99, 0, 18, R1 = 18, tau = tau, stones = stones, ind_var = 'on', var_sd = sd)
imnmrug = filepath + "images/ind var/rugatulus-indvar-sd" + str(round(sd*100)) + "-tau" + str(round(tau*1000)) +  "-s99-stones" + str(stones) + ".png"
imgrug = outrug3000[6]
plt.imshow(imgrug, cmap = 'gray')
plt.colorbar()
plt.axis('off')
plt.savefig(imnmrug, dpi = 350)
plt.clf()

stones = 1000
outrug1000 = w.rugatulus(99, 0, 18, R1 = 18, tau = tau, stones = stones, ind_var = 'on', var_sd = sd)
imnmrug = filepath + "images/ind var/rugatulus-indvar-sd" + str(round(sd*100)) + "-tau" + str(round(tau*1000)) +  "-s99-stones" + str(stones) + ".png"
imgrug = outrug1000[6]
plt.imshow(imgrug, cmap = 'gray')
plt.colorbar()
plt.axis('off')
plt.savefig(imnmrug, dpi = 350)
plt.clf()


