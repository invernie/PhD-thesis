import pickle
import numpy as np
import matplotlib.pyplot as plt
import wall_funcs as w


filepath = "path/to/folder/"

S = 100
stones = 3000


## mutation in germline - Ropt
qrange = np.linspace(0.05,0.5,10)
R2range = np.linspace(8,24,9)
#qrange = np.linspace(0,1,1) # for debugging
#R2range = np.linspace(8,10,2) # for debugging

for q in qrange:
    
    for R2 in R2range:
        
        outListFD = list()
        outListT = list()      
        outListr = list()
        
        for s in range(S):
            if s == (S-1):
                
                outFD = w.FD(s,q,R2, stones=stones)
                outListFD.append(outFD)
                imnmFD = filepath + "FD-R" + str(round(R2)) + "-q" + str(round(q*100)) + "-s99-stones" + str(stones) + ".png"
                imgFD = outFD[6]
                plt.imshow(imgFD, cmap = 'gray')
                plt.colorbar()
                plt.axis('off')
                plt.savefig(imnmFD, dpi = 350)
                plt.clf()
                
                outT = w.templateOnly(s,q,R2, stones=stones)
                outListT.append(outT)
                imnmT = filepath + "template-R" + str(round(R2)) + "-q" + str(round(q*100)) + "-s99-stones" + str(stones) + ".png"
                imgT = outT[6]
                plt.imshow(imgT, cmap = 'gray')
                plt.colorbar()
                plt.axis('off')
                plt.savefig(imnmT, dpi = 350)
                plt.clf()
                
                outr = w.rugatulus(s,q,R2, stones=stones)
                outListr.append(outr)
                imnmr = filepath + "rugatulus-R" + str(round(R2)) + "-q" + str(round(q*100)) + "-s99-stones" + str(stones) + ".png"
                imgr = outr[6]
                plt.imshow(imgr, cmap = 'gray')
                plt.colorbar()
                plt.axis('off')
                plt.savefig(imnmr, dpi = 350)
                plt.clf()
                
            else:
                outListFD.append(w.FD(s,q,R2, stones=stones))            
                outListT.append(w.templateOnly(s,q,R2, stones=stones))
                outListr.append(w.rugatulus(s,q,R2, stones=stones))
        
        filenameFD = filepath + "out_FD_q" + str(round(q*100)) + "_R" + str(round(R2))
        fileFD = open(filenameFD, 'wb')
        pickle.dump(outListFD, fileFD)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
        fileFD.close()
        
        filenameT = filepath + "out_T_q" + str(round(q*100)) + "_R" + str(round(R2))
        fileT = open(filenameT, 'wb')
        pickle.dump(outListT, fileT)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
        fileT.close()
        
        filenamer = filepath + "out_rug_q" + str(round(q*100)) + "_R" + str(round(R2))
        filer = open(filenamer, 'wb')
        pickle.dump(outListr, filer)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
        filer.close()
        
        print("q is at " + str(q) + " and R2 is at: " + str(R2))
        

# analyse results
from pandas import DataFrame

statListFD = list()
statListT = list()
statListr = list()
fsteps = 10
for q in qrange:
    
    for R2 in R2range:
        
        loadpathFD = filepath + "out_FD_q" + str(round(q*100)) + "_R" + str(round(R2))
        fileFD = open(loadpathFD, 'rb')
        allResFD = pickle.load(fileFD)
        fileFD.close()
        
        CoVList_FD = [sum(f[4][-fsteps:])/fsteps for f in allResFD]
        aveCoV_FD = np.mean(CoVList_FD)
        stdCoV_FD = np.std(CoVList_FD)
        RbarList_FD = [sum(f[5][-fsteps:])/fsteps for f in allResFD]
        aveRbar_FD = np.mean(RbarList_FD)
        stdRbar_FD = np.std(RbarList_FD)
        statListFD.append(["FD",q,R2,aveCoV_FD,stdCoV_FD,aveRbar_FD,stdRbar_FD])
        
        
        loadpathT = filepath + "out_T_q" + str(round(q*100)) + "_R" + str(round(R2))
        fileT = open(loadpathT, 'rb')
        allResT = pickle.load(fileT)
        fileT.close()
        
        CoVList_T = [sum(f[4][-fsteps:])/fsteps for f in allResT]
        aveCoV_T = np.mean(CoVList_T)
        stdCoV_T = np.std(CoVList_T)
        RbarList_T = [sum(f[5][-fsteps:])/fsteps for f in allResT]
        aveRbar_T = np.mean(RbarList_T)
        stdRbar_T = np.std(RbarList_T)
        statListT.append(["T",q,R2,aveCoV_T,stdCoV_T,aveRbar_T,stdRbar_T])
        
        
        loadpathr = filepath + "out_rug_q" + str(round(q*100)) + "_R" + str(round(R2))
        filer = open(loadpathr, 'rb')
        allResr = pickle.load(filer)
        filer.close()
        
        CoVList_r = [sum(f[4][-fsteps:])/fsteps for f in allResr]
        aveCoV_r = np.mean(CoVList_r)
        stdCoV_r = np.std(CoVList_r)
        RbarList_r = [sum(f[5][-fsteps:])/fsteps for f in allResr]
        aveRbar_r = np.mean(RbarList_r)
        stdRbar_r = np.std(RbarList_r)
        statListr.append(["rug",q,R2,aveCoV_r,stdCoV_r,aveRbar_r,stdRbar_r])
        
statListAll = statListFD + statListT + statListr
df = DataFrame(statListAll, columns = ["model",'q','R2','CoV_ave','CoV_SD','Rbar_ave', 'Rbar_SD'])
svname = filepath + 'Rq_sweep_3000stones.csv'
df.to_csv(svname)



#################
# change in tau #
#################

filepath = "path/to/folder/"

stones = 1000
tauList = [0.025,0.075]
qList = [0.1,0.3,0.5]
R2List = [8,12,16]

s = 99

for tau in tauList:
    
    for q in qList:
        for R2 in R2List:
            
            outT = w.templateOnly(s,q,R2, tau = tau, stones=stones)
            imnmT = filepath + "images/template-R" + str(round(R2)) + "-q" + str(round(q*100)) + "-tau" + str(round(tau*1000)) + "-s99-stones" + str(stones) + ".png"
            imgT = outT[6]
            plt.imshow(imgT, cmap = 'gray')
            plt.colorbar()
            plt.axis('off')
            plt.savefig(imnmT, dpi = 350)
            plt.clf()
            
            outr = w.rugatulus(s,q,R2, tau = tau, stones=stones)
            imnmr = filepath + "images/rugatulus-R" + str(round(R2)) + "-q" + str(round(q*100)) + "-tau" + str(round(tau*1000)) + "-s99-stones" + str(stones) + ".png"
            imgr = outr[6]
            plt.imshow(imgr, cmap = 'gray')
            plt.colorbar()
            plt.axis('off')
            plt.savefig(imnmr, dpi = 350)
            plt.clf()
       

stones = 3000    
tau = 0.075
for q in qList:
    for R2 in R2List:
        
        outT = w.templateOnly(s,q,R2, tau = tau, stones=stones)
        imnmT = filepath + "images/3000 stones/template-R" + str(round(R2)) + "-q" + str(round(q*100)) + "-tau" + str(round(tau*1000)) + "-s99-stones" + str(stones) + ".png"
        imgT = outT[6]
        plt.imshow(imgT, cmap = 'gray')
        plt.colorbar()
        plt.axis('off')
        plt.savefig(imnmT, dpi = 350)
        plt.clf()
            
        outr = w.rugatulus(s,q,R2, tau = tau, stones=stones)
        imnmr = filepath + "images/3000 stones/rugatulus-R" + str(round(R2)) + "-q" + str(round(q*100)) + "-tau" + str(round(tau*1000)) + "-s99-stones" + str(stones) + ".png"
        imgr = outr[6]
        plt.imshow(imgr, cmap = 'gray')
        plt.colorbar()
        plt.axis('off')
        plt.savefig(imnmr, dpi = 350)
        plt.clf()
 
        

