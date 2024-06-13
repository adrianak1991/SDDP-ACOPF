In order to run the SDDP algoritm for the case study, run the file: SDDP_parallel.m . 

At the begining of this file, you can set the parameters of the case, like number of stages, number of iterations etc. However, changing number of nodes per stage or number of stages, substages requires creating new lattice which can be accomplished only if you possess QUASAR licence. Then you have to set the value of the parameter 'param.new_fit' to 'true'.

In order to run out of sample comparison of results and see how the solution is improved with recovering of feasible voltages, you can run the file: comparison.m .