Python script for dynamically load balancing parallel simulations.

Copy dynamicLoadBalancing.py to your simulation folder.

Then run it as: python3 dynamicLoadBalancing.py -i=maximumAllowedImbalance -w=weightField -l=logFile -s=solver

The default user defined variables are:
maximumAllowedImbalance = 10%
logFileName = log.uspFoam
solver = uspFoam
The weightField variable must always be given as an argument.

Eaxmple:
./dynamicLoadBalancing -i=5 -w=uspRhoN_Ar -l=log.uspFoam -s=uspFoam

For this example:
The maximum allowed imbalance is 5%.
The weightField used for rebalancing is uspRhoN_Ar.
The log file name where the imbalance is read from is log.uspFoam.
The solver that is rebalanced is uspFoam.


