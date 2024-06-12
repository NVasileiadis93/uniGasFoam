Python script for dynamically load balancing parallel simulations.

Copy dynamicLoadBalancing.py to your simulation folder.

Then run it as: python3 dynamicLoadBalancing.py -i=maximumAllowedImbalance -w=weightField -l=logFile -s=solver

The default user defined variables are:
maximumAllowedImbalance = 10%
logFileName = log.uniGasFoam
solver = uniGasFoam
The weightField variable must always be given as an argument.

Eaxmple:
./dynamicLoadBalancing -i=5 -w=uniGasRhoN_Ar -l=log.uniGasFoam -s=uniGasFoam

For this example:
The maximum allowed imbalance is 5%.
The weightField used for rebalancing is uniGasRhoN_Ar.
The log file name where the imbalance is read from is log.uniGasFoam.
The solver that is rebalanced is uniGasFoam.


