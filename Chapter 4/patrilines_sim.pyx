# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 21:26:06 2021

@author: Edith Invernizzi
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import wall_funcs as w

####################
## four patrilines #
####################

filepath = "path/to/folder/"

S = 10
stones = 3000
#sd_list = [1] # for debugging
#flistList = [[0.25,0.25,0.25]]
sd_list = [0.5,2,5]
flistList = [[0.25,0.25,0.25], [0.5,0.2,0.2], [0.5,0.4,0.05]]

for sd in sd_list:
    
    for flist in flistList:
        
        f1 = flist[0]
        f2 = flist[1]
        f3 = flist[2]
        
        TList = list()
        rugList = list()
        for s in range(S):
            
            outT = w.templateOnly_patrilines(s, 18, sd, f1 = f1, f2 = f2, f3 = f3)
            outrug = w.rugatulus_patrilines(s, 18, sd, f1 = f1, f2 = f2, f3 = f3)
            #TList.append(outT)
            #rugList.append(outrug)
            
            imnmT = filepath + "images/template-patrilines-sd"+ str(round(sd*100)) + "-s" + str(s) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100)) + str(stones) + ".png"
            imgT = outT[6]
            plt.imshow(imgT, cmap = 'gray')
            plt.colorbar()
            plt.axis('off')
            plt.savefig(imnmT, dpi = 350)
            plt.clf()
            
            imnmr = filepath + "images/rugatulus-patrilines-sd" + str(round(sd*100)) + "-s" + str(s) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100)) + str(stones) + ".png"
            imgrug = outrug[6]
            plt.imshow(imgrug, cmap = 'gray')
            plt.colorbar()
            plt.axis('off')
            plt.savefig(imnmr, dpi = 350)
            plt.clf()

                
    
        #filenameT = filepath + "out_T_sd" + str(round(sd*100)) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100))
        #fileT = open(filenameT, 'wb')
        #pickle.dump(TList, fileT)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
        #fileT.close()
        
        #filenamer = filepath + "out_rug_sd" + str(round(sd*100)) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100))
        #filer = open(filenamer, 'wb')
        #pickle.dump(rugList, filer)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
        #filer.close()
        
        print("sd is at " + str(sd) + " and flist is at " + str(flist))
        


# analyse results
import pandas as pd
from pandas import DataFrame

statListT = list()
statListr = list()
fsteps = 10
for sd in sd_list:
    
    sddt = {'model':[],'SD':[],'f1':[],'f2':[],'f3':[],'R1':[],'R2':[],'R3':[],'R4':[],'CoV':[],'Rbar':[]}
    sdDF = DataFrame(sddt)
    for flist in flistList:
        
        f1 = flist[0]
        f2 = flist[1]
        f3 = flist[2]
    
        loadpathT = filepath + "out_T_sd" + str(round(sd*100)) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100))
        fileT = open(loadpathT, 'rb')
        allResT = pickle.load(fileT)
        fileT.close()
        
        CoVList_T = [sum(f[4][-fsteps:])/fsteps for f in allResT]
        aveCoV_T = np.mean(CoVList_T)
        stdCoV_T = np.std(CoVList_T)
        RbarList_T = [sum(f[5][-fsteps:])/fsteps for f in allResT]
        aveRbar_T = np.mean(RbarList_T)
        stdRbar_T = np.std(RbarList_T)
        Rlist_T = [i[7] for i in allResT]
        R1list_T = [R1[0] for R1 in Rlist_T]
        min_R1_T = np.min(R1list_T)
        max_R1_T = np.max(R1list_T)
        R2list_T = [R2[1] for R2 in Rlist_T]
        min_R2_T = np.min(R2list_T)
        max_R2_T = np.max(R2list_T)
        R3list_T = [R3[2] for R3 in Rlist_T]
        min_R3_T = np.min(R3list_T)
        max_R3_T = np.max(R3list_T)
        R4list_T = [R4[3] for R4 in Rlist_T]
        min_R4_T = np.min(R4list_T)
        max_R4_T = np.max(R4list_T)
        statListT.append(["T",sd,f1,f2,f3,min_R1_T,max_R1_T,min_R2_T,max_R2_T,min_R3_T,max_R3_T,min_R4_T,max_R4_T,aveCoV_T,stdCoV_T,aveRbar_T,stdRbar_T])
        
        
        loadpathr = filepath + "out_rug_sd" + str(round(sd*100)) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100))
        filer = open(loadpathr, 'rb')
        allResr = pickle.load(filer)
        filer.close()
        
        CoVList_r = [sum(f[4][-fsteps:])/fsteps for f in allResr]
        aveCoV_r = np.mean(CoVList_r)
        stdCoV_r = np.std(CoVList_r)
        RbarList_r = [sum(f[5][-fsteps:])/fsteps for f in allResr]
        aveRbar_r = np.mean(RbarList_r)
        stdRbar_r = np.std(RbarList_r)
        Rlist_r = [i[7] for i in allResr]
        R1list_r = [R1[0] for R1 in Rlist_r]
        min_R1_r = np.min(R1list_r)
        max_R1_r = np.max(R1list_r)
        R2list_r = [R2[1] for R2 in Rlist_r]
        min_R2_r = np.min(R2list_r)
        max_R2_r = np.max(R2list_r)
        R3list_r = [R3[2] for R3 in Rlist_r]
        min_R3_r = np.min(R3list_r)
        max_R3_r = np.max(R3list_r)
        R4list_r = [R4[3] for R4 in Rlist_r]
        min_R4_r = np.min(R4list_r)
        max_R4_r = np.max(R4list_r)
        statListr.append(["rug",sd,f1,f2,f3,min_R1_r,max_R1_r,min_R2_r,max_R2_r,min_R3_T,max_R3_T,min_R4_T,max_R4_T,aveCoV_r,stdCoV_r,aveRbar_r,stdRbar_r])
        
        
        data = {'model':["T"]*S + ["rug"]*S,'SD':[sd]*2*S,'f1':[f1]*2*S,'f2':[f2]*2*S,'f3':[f3]*2*S,'R1':R1list_T + R1list_r,'R2':R2list_T + R2list_r,'R3':R3list_T + R3list_r,'R4':R4list_T + R4list_r,'CoV':CoVList_T + CoVList_r,'Rbar':RbarList_T + RbarList_r}
        valueDF = DataFrame(data)
        sdDF = pd.concat([sdDF,valueDF])
        
    #svvalnm = filepath + "allsimvals-sd" + str(round(sd*100)) + "-f1-"+ str(round(f1*100)) + "-f2-" + str(round(f2*100)) + ".csv"
    svvalnm = filepath + "sdsimvals-sd" + str(round(sd*100)) + ".csv"
    sdDF.to_csv(svvalnm, index = False)
        
        
statListAll = statListT + statListr
df = DataFrame(statListAll, columns = ["model",'sd', 'f1', 'f2', 'f3', 'minR1', 'maxR1', 'minR2', 'maxR2', 'minR3', 'maxR3', 'minR4', 'maxR4', 'CoV_ave','CoV_SD','Rbar_ave', 'Rbar_SD'])
svname = filepath + 'four-patrilines_mean18_3000stones.csv'
df.to_csv(svname, index = False)


###############################################
# four patrilines with spatial specialisation #
###############################################

filepath = "path/to/folder/"

S = 10
stones = 3000
T = 10000
#sd_list = [1] # for debugging
#flistList = [[0.25,0.25,0.25]]
sd_list = [0.5,2,5]
flistList = [[0.25,0.25,0.25], [0.5,0.2,0.2], [0.5,0.4,0.05]]


for sd in sd_list:
    
    for flist in flistList:
        
        f1 = flist[0]
        f2 = flist[1]
        f3 = flist[2]
        
        TList = list()
        rugList = list()
        for s in range(S):
            
            outT = w.templateOnly_patrilines(s, 18, sd, f1 = f1, f2 = f2, f3 = f3, T = T, spatSpec = 'on')
            outrug = w.rugatulus_patrilines(s, 18, sd, f1 = f1, f2 = f2, f3 = f3, T = T, spatSpec = 'on')
                
            #TList.append(outT)
            #rugList.append(outrug)
            
            imnmT = filepath + "images/template-patrilines-spsp-sd"+ str(round(sd*100)) + "-s" + str(s) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100)) + str(stones) + "-T" + str(T) + ".png"
            imgT = outT[6]
            plt.imshow(imgT, cmap = 'gray')
            plt.colorbar()
            plt.axis('off')
            plt.savefig(imnmT, dpi = 350)
            plt.clf()
            
            imnmr = filepath + "images/rugatulus-patrilines-spsp-sd" + str(round(sd*100)) + "-s" + str(s) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100)) + str(stones) + "-T" + str(T) + ".png"
            imgrug = outrug[6]
            plt.imshow(imgrug, cmap = 'gray')
            plt.colorbar()
            plt.axis('off')
            plt.savefig(imnmr, dpi = 350)
            plt.clf()
            
            
    
        #filenameT = filepath + "out_T_sd" + str(round(sd*100)) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100))
        #fileT = open(filenameT, 'wb')
        #pickle.dump(TList, fileT)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
        #fileT.close()
        
        #filenamer = filepath + "out_rug_sd" + str(round(sd*100)) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100))
        #filer = open(filenamer, 'wb')
        #pickle.dump(rugList, filer)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
        #filer.close()
        
        print("sd is at " + str(sd) + " and flist is at " + str(flist))
        


# analyse results
from pandas import DataFrame
import pandas as pd

statListT = list()
statListr = list()
fsteps = 10
for sd in sd_list:
    
    sddt = {'model':[],'SD':[],'f1':[],'f2':[],'f3':[],'R1':[],'R2':[],'R3':[],'R4':[],'CoV':[],'Rbar':[]}
    sdDF = DataFrame(sddt)
    for flist in flistList:
        
        f1 = flist[0]
        f2 = flist[1]
        f3 = flist[2]
    
        loadpathT = filepath + "out_T_sd" + str(round(sd*100)) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100))
        fileT = open(loadpathT, 'rb')
        allResT = pickle.load(fileT)
        fileT.close()
        
        CoVList_T = [sum(f[4][-fsteps:])/fsteps for f in allResT]
        aveCoV_T = np.mean(CoVList_T)
        stdCoV_T = np.std(CoVList_T)
        RbarList_T = [sum(f[5][-fsteps:])/fsteps for f in allResT]
        aveRbar_T = np.mean(RbarList_T)
        stdRbar_T = np.std(RbarList_T)
        Rlist_T = [i[7] for i in allResT]
        R1list_T = [R1[0] for R1 in Rlist_T]
        min_R1_T = np.min(R1list_T)
        max_R1_T = np.max(R1list_T)
        R2list_T = [R2[1] for R2 in Rlist_T]
        min_R2_T = np.min(R2list_T)
        max_R2_T = np.max(R2list_T)
        R3list_T = [R3[2] for R3 in Rlist_T]
        min_R3_T = np.min(R3list_T)
        max_R3_T = np.max(R3list_T)
        R4list_T = [R4[3] for R4 in Rlist_T]
        min_R4_T = np.min(R4list_T)
        max_R4_T = np.max(R4list_T)
        statListT.append(["T",sd,f1,f2,f3,min_R1_T,max_R1_T,min_R2_T,max_R2_T,min_R3_T,max_R3_T,min_R4_T,max_R4_T,aveCoV_T,stdCoV_T,aveRbar_T,stdRbar_T])
        
        
        loadpathr = filepath + "out_rug_sd" + str(round(sd*100)) + "f1-" + str(round(f1*100)) + "-f2-" + str(round(f2*100))
        filer = open(loadpathr, 'rb')
        allResr = pickle.load(filer)
        filer.close()
        
        CoVList_r = [sum(f[4][-fsteps:])/fsteps for f in allResr]
        aveCoV_r = np.mean(CoVList_r)
        stdCoV_r = np.std(CoVList_r)
        RbarList_r = [sum(f[5][-fsteps:])/fsteps for f in allResr]
        aveRbar_r = np.mean(RbarList_r)
        stdRbar_r = np.std(RbarList_r)
        Rlist_r = [i[7] for i in allResr]
        R1list_r = [R1[0] for R1 in Rlist_r]
        min_R1_r = np.min(R1list_r)
        max_R1_r = np.max(R1list_r)
        R2list_r = [R2[1] for R2 in Rlist_r]
        min_R2_r = np.min(R2list_r)
        max_R2_r = np.max(R2list_r)
        R3list_r = [R3[2] for R3 in Rlist_r]
        min_R3_r = np.min(R3list_r)
        max_R3_r = np.max(R3list_r)
        R4list_r = [R4[3] for R4 in Rlist_r]
        min_R4_r = np.min(R4list_r)
        max_R4_r = np.max(R4list_r)
        statListr.append(["rug",sd,f1,f2,f3,min_R1_r,max_R1_r,min_R2_r,max_R2_r,min_R3_T,max_R3_T,min_R4_T,max_R4_T,aveCoV_r,stdCoV_r,aveRbar_r,stdRbar_r])
        
        data = {'model':["T"]*S + ["rug"]*S,'SD':[sd]*2*S,'f1':[f1]*2*S,'f2':[f2]*2*S,'f3':[f3]*2*S,'R1':R1list_T + R1list_r,'R2':R2list_T + R2list_r,'R3':R3list_T + R3list_r,'R4':R4list_T + R4list_r,'CoV':CoVList_T + CoVList_r,'Rbar':RbarList_T + RbarList_r}
        valueDF = DataFrame(data)
        sdDF = pd.concat([sdDF,valueDF])
        
    svvalnm = filepath + "spsp-sdsimvals-sd" + str(round(sd*100)) + ".csv"
    sdDF.to_csv(svvalnm, index = False)
        
        
statListAll = statListT + statListr
df = DataFrame(statListAll, columns = ["model",'sd', 'f1', 'f2', 'f3', 'minR1', 'maxR1', 'minR2', 'maxR2', 'minR3', 'maxR3', 'minR4', 'maxR4', 'CoV_ave','CoV_SD','Rbar_ave', 'Rbar_SD'])
svname = filepath + 'four-patrilines_spatial-specialisation_mean18_3000stonesT10000.csv'
df.to_csv(svname)

   