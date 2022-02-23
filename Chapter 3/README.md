This sub-repository contains the files used for the statistical analysis of ant activity data, performed as part of my Ph.D.. A more tailored version of this analysis was later developed with the help of Dr Valentin Popov and Dr Théo Michelot (CREEM, University of St Andrews) and is quickly approaching the submission-for-publication stage. The analysis contained here was written and performed by me.

The analysis run with these files uses hidden Markov models. This approach enables to split observations from the same time series into multiple distribution of origins (corresponding to different underlying individual "states" that are believed to have generated the observation), for cases where the difference between these distributions is large enough that this choice of approach generates a better, albeit more complex, model. The underlying process generating the states is modelled using a Markov chain with uknown (hidden) states that need to be inferred from the data. Regression analysis can be incorporated into this method to fit separate models to the data distributions.

Below you will find a description of the files.

### Data
`rates.csv`: processed data, on which the analysis is run
  
 `for data-processing replication` (folder): raw data files  
   
  `analyseData.R`: script used for data processing  (combines data files, calculates derived measures, _etc._)
    
    
 ### HMM analysis and model selection   
 `best_starting_values_HMMs.R`: fits HMM models over 100 iterations with different starting values, to find the global minimum (_i.e._, the maximum of the likelihood function). Uses the `momentuHMM` package.  
  
 `Activity analysis - overall variables - best model selection.rmd` and homonymous pdf file:  describes the model selection process
