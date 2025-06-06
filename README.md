*This reposotory contains data and code for the following manuscript:*

**Borden, J.B., J.P. Gibbs, J.P. Vanek, and B.J. Cosentino. Beyond urbanization metrics: Using graphical causal models to investigate mechanisms in urban ecology and evolution**

*The following files are included:*

#### measuredVariables.R

R Code to estimate direct effects among measured variables based on the DAG in Fig. 3A. 

#### measuredVariables.RData

Measured environmental variables and predator activity data at each site. Variables include building cover (`bld`, proportion), distance to city center (`dis`, km), forest fragmentation (`frg`, disjunct core area density), impervious cover (`imp`, proportion), human population density (`pdn`, persons per square km), predator activity (`prd`, percentage of days when predator detected), road cover (`rdl.c`, proportion), and tree cover (`tre`, proportion). The data frame `rai.sum` includes predator detection histories (detections by day by site, total number of days with detections, and total number of days the camera at each site was active).

#### squirrleModels.R

R code used to estimate direct effects, total effects, and univariate models of proportion melanic squirrels at each site.

#### squirrelObservations.RData

Squirrel detections for gray (`y.occ.g`) and melanic (`y.occ.m`) color morphs, squirrel counts for gray (`y.ct.g` and melanic (`y.ct.m`) color morphs, and temperature data used to fit models proportion melanic squirrels at each site (`temp.avg.ct` for count surveys and `temp.avg.occ` for camera surveys).

