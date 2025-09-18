from optparse import OptionParser
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.integrate as integ

#pd.options.mode.chained_assignment = None

def entropy(x): # Use base-2 log (bits) for consistency with JS divergence.
      
  if ((len(x) <= 1) or (not np.count_nonzero(x))):
    return 0

  x = x[x != 0]
  s = -np.sum(x*np.log2(x))
  return s

def gini(x):

    mg = np.mean(x)
    if (mg == 0):
      g = 0

    else:
      mad = np.abs(np.subtract.outer(x, x)).mean()
      g = 0.5 * mad/np.mean(x)
      
    return g
  
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("--covin", dest="COVIN", help="Input coverage bed file name")
    parser.add_option("--mmin", dest="MMIN", help="Input mismatch bed file name")
    parser.add_option("--outfile", dest="OUTFILE", help="Output file/path")
    
    (options, args) = parser.parse_args()
    
    COVIN = pd.read_csv(options.COVIN, sep='\t', index_col='name')
    MMIN = pd.read_csv(options.MMIN, sep='\t', index_col='name')
    BB = pd.merge(COVIN,MMIN, on=['chr','start','end','name','strand'])

    testHead = ["chr","start","end","name","length","numSites","delta","meanCov1","meanCov2","meanCov3","meanCov4","meanCov5","meanCov6","log2FCcoverage","pvalCoverage","meanMMRate1","meanMMRate2","meanMMRate3","meanMMRate4","meanMMRate5","meanMMRate6","log2FCMMRate","pvalMMrate","entropy1","entropy2","entropy3","entropy4","entropy5","entropy6","r","pvalr","gini1","gini2","gini3","gini4","gini5","gini6","log2FCgini","pvalGini"]
          
    out = open(options.OUTFILE,"w")
    out.write(",".join(testHead) + "\n")

    p = 0.00001 # Set a pseudocount for logarithms.
    upper = 0.05 # Set top range of values to Winsorize.
    for row in BB.iterrows(): # Each row is a genomic region.
      sites = len(row[1].coverage1.split(","))
      length = abs(row[1].end - row[1].start)
      coverage = np.empty([sites,6],dtype=float)
      rates = np.empty([sites,6],dtype=float)
      wrates = np.empty([sites,6],dtype=float)
      swrates = np.empty([sites,6],dtype=float)
        
      coverage[:,0] = np.array(row[1].coverage1.split(",")).astype(float)
      coverage[:,1] = np.array(row[1].coverage2.split(",")).astype(float)
      coverage[:,2] = np.array(row[1].coverage3.split(",")).astype(float)
      coverage[:,3] = np.array(row[1].coverage4.split(",")).astype(float)
      coverage[:,4] = np.array(row[1].coverage5.split(",")).astype(float)
      coverage[:,5] = np.array(row[1].coverage6.split(",")).astype(float)
      rates[:,0] = np.array(row[1].mmRate1.split(",")).astype(float)
      rates[:,1] = np.array(row[1].mmRate2.split(",")).astype(float)
      rates[:,2] = np.array(row[1].mmRate3.split(",")).astype(float)
      rates[:,3] = np.array(row[1].mmRate4.split(",")).astype(float)
      rates[:,4] = np.array(row[1].mmRate5.split(",")).astype(float)
      rates[:,5] = np.array(row[1].mmRate6.split(",")).astype(float)        

      for i in range(6): # Winsorize within transcript.
        wrates[:,i] = stats.mstats.winsorize(rates[:,i],limits=[0,upper])

      for i in range(6): # Rescale rates within transcript.
        scale = np.max(wrates[:,i])
        if scale > 0: # Control for uniform zero rate cases.
          swrates[:,i] = wrates[:,i]/scale

        else:
          swrates[:,i] = wrates[:,i]

      deltai = (1.0/3.0)*(np.abs(rates[:,0]-rates[:,3]) # Compute per-base d.
               + np.abs(rates[:,1]-rates[:,4])
               + np.abs(rates[:,2]-rates[:,5]))
      delta = np.mean(deltai) # Compute per-region D.
      
      within = np.empty(6,dtype=float)
      between = np.empty(9,dtype=float)
      rwithin = np.empty(6,dtype=float)
      rbetween = np.empty(9,dtype=float)

      if len(rates[:,i]) > 1:
        k=0
        for i in range(2):
          for j in range(i+1,3):
            rwithin[k] = stats.pearsonr(swrates[:,i],swrates[:,j]).statistic
            k = k+1

        for i in range(3,5):
          for j in range(i+1,6):
            rwithin[k] = stats.pearsonr(swrates[:,i],swrates[:,j]).statistic
            k = k+1
          
        k=0
        for i in range(3):
          for j in range(3,6):
            rbetween[k] = stats.pearsonr(swrates[:,i],swrates[:,j]).statistic
            k = k+1

        meanRBetween = np.mean(rbetween) # Measure and test r between cond.
        pRBetween = stats.ranksums(rbetween,rwithin,alternative='less').pvalue

      else:
        meanRBetween = np.nan
        pRBetween = np.nan
        
      meanCov = np.empty(6,dtype=float)
      meanMMR = np.empty(6,dtype=float)
      entr = np.empty(6,dtype=float)
      gnt = np.empty(6,dtype=float)
      for i in range(6):
        meanCov[i] = np.mean(coverage[:,i])
        meanMMR[i] = np.mean(rates[:,i])
        entr[i] = entropy(swrates[:,i])
        gnt[i] = gini(swrates[:,i])
        
      covrp = stats.ttest_rel(meanCov[0:3], meanCov[3:6]).pvalue 
      mmrp = stats.ttest_rel(meanMMR[0:3], meanMMR[3:6]).pvalue
      ginip = stats.ttest_rel(gnt[0:3], gnt[3:6]).pvalue
      lrcov = np.log2(np.mean(meanCov[3:6])+p)-np.log2(np.mean(meanCov[0:3])+p)
      lrrate =np.log2(np.mean(meanMMR[3:6])+p)-np.log2(np.mean(meanMMR[0:3])+p)
      lrgini = np.log2(np.mean(gnt[3:6])+p)-np.log2(np.mean(gnt[0:3])+p)
      mci = ','.join(map(str, np.ndarray.tolist(meanCov)))
      mri = ','.join(map(str, np.ndarray.tolist(meanMMR)))
      entri = ','.join(map(str, np.ndarray.tolist(entr)))
      gnti = ','.join(map(str, np.ndarray.tolist(gnt)))
      testList = [row[1].chr,str(row[1].start),str(row[1].end),row[1].name,str(length),str(sites),str(delta),mci,str(lrcov),str(covrp),mri,str(lrrate),str(mmrp),entri,str(meanRBetween),str(pRBetween),gnti,str(lrgini),str(ginip)]

      out.write(",".join(testList) + "\n")

    out.close()
