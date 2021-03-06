# E_coli_Evolution  (Still Updating)
**How to run:**    
Step 1: download all files;  
Step 2: run R;  
Step 3: open "Run_shell_C.R", reset the working directory by the setwd command in the first line;  
Step 4: run program with "Rscript Run_shell_C.R *path_to_input_file seed traits threshold iteration_steps Num_expriment Num_replicate"  

**Parameters:**   
path_to_input_file: the location where the database is  
seed: random seed  
traits: number of traits  
threshold: stop the simulation when the objective reaches the threshold  
iteration_steps: number of iterations in the simulation  
Num_experiment: which experiment is excluded from the training set and used for testing  
Num_replicate: an experiment may have several biological replicates, what are their numbers in the database  


**Files generated:**  
run_report.txt: report the performance of the model every 10K steps  
training.jpg: plot the the origial data of training set and the corresponding simulating results  
testing.jpg: plot the test results 

**Files (functions) used:**    
**Merged_Data_4.csv**: the database with bacteria gene expression profile and environmental parameters  
**Main_new_C.R**: the main workflow of simulation  
**Main_new_C_par.R**: the main workflow of simulation, implemented with Rcpp package parallel computing features  
**Train_SA_C.cpp**: the training function using simulated annealing to find out a optimal solution, implemented in C++ to accelerate computing  
**Func_Support.R**: all supportive functions, including:  
*Kfold, None_zero_idx_same, Shuffle*: divide the data according to the papers they are from. For each paper that contributes num_exp samples to the dataset, a floor(num_exp/K) will be kept for testing  
*Validate*: apply the G and E matrices to testing data  
*Draw*: draw the figures of "training" and "testing"  
**GetpBad.cpp**: main part of simulated annealing, using the current possibility of making "bad" predictions to decide whether to update an element in the gene expression and environmental parameter matrices. Implemented in C++ to accelerate the computation  
**Getweight.cpp**: assign weights to each data point in the training dataset, based on how frequent its value (i.e. growth rates) is seen in the database. For example, most bacterial growth rates are between 0.5 and 0.8, then data points out of this range will be given larger weights so that the model will not only focus on those more common values    

**Brief introdution:**  
The growth rates of E.coli are co-determined by the environmental conditions and genetics. Some environmental conditions (supplementation of nutrients and maintenance of physiological pH in the medium) are beneficial, some are not. E.coli also constantly evolves by mutating their genomes. In general, the mutations are random, leading to beneficial or detrimental outcomes, but those carrying beneficial mutations could outgrow those carrying bad ones and become dominant. Therefore, the growth rate is a good metric of microorganisms' adaption to the environments. Modeling the interaction of environments and genetics to predict the growth rates under such situations is our goal.
  
Each experiment measuring the growth rate of an E.coli strain can be divided into two parts: environment and genetics. Both parts are converted to vectors. The first part has the information such as the medium used and the concentrations of nutrients added. The second part records all the gene mutations occured in this strain compared with the reference genome. Four possible mutations could occur to a gene: insertion, SNP, replication and deletion. Any mutation found to exist in the genome will be given value "1" in the vector, and "0" for not having such a mutation. The genetics vector is named X and environment one Y.    
  
Each gene or environmental condition may influence some traits of E.coli. For example, adding more glucose may facilitate cell division, but a hostile pH may slow it down again. Hence we define the gene-trait matrix G and environment-trait matrix E. G and E have the same number of rows, which corresponds to the number of traits. On the other hand, their column numbers correspond to the number of possible gene mutations environmental conditions. So G[a, b] describe the influence of gene mutaion b on trait a. Similarly, E[c, d] describes environmental condition d's influence on trait c.  
  
The simulation of growth data has two steps. First, TG = G * X and TE = E * Y are calculated, generating two vectors TG and TE with the length of trait number. Here, TG describes the overall strength of influences on all traits of a E.coli strain under such a gene mutation profile, and TE describes the overall strength of influences on all traits under a given environment. Then, the simulated growth rate is calculated as the dot product of TG and TE.  
  
Clearly, the aim of the study is to determine the G and E matrices. Their initial values are assigned randomly from a normal distribution. Then the optimization is given by simulated annealing so that the simulated growth data are closer and closer to the lab data. Good simulation dominated by G and E not only matches the growth rates of the training set, but could also predict the growth rates correctly based on the environmental conditions and genetic profiles provided. However, the simulated annealing is too slow in R so conversion to C/C++ to accelerate the running is needed.  
