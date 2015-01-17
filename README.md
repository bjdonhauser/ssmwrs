# README.md #

## Summary ##
* Matlab scripts to estimate state-space models with markov-switching.
* Intended as companion scripts to Kim, C. J., & Nelson, C. R. (1999). State-space models with regime switching. *Cambridge, Mass*.
* Tested on Matlab 2007. 

## Main Scripts ##
* `ssm_optimizer.m`: The main engine used to estimate markov-switching state-space models. This is intended to be run from a file 'model_main' where "model" is the identifier for a particular model (e.g., "lam").
* `ssm_objective.m`: The objective function that the function `ssm_optimizer` optimizes. 

## Example Model Scripts ##
The following scripts represent Lam's Generalized Hamilton Model so that its parameters can be estimated via the main scripts. To estimate alternative state-space models with markov-switching, one must write files `<myaltmodel>_main.m`, `<myaltmodel>_ssf.m`, `<myaltmodel>_constraints.m` and follow the conventions found in the `lam` scripts.
* `/models/lam/lam_main.m`: The main method for the Lam's Generalized Hamilton Model as presented in Kim and Nelson (1999, pg. 111). Returns parameter estimates, standard errors, state inferences, state variable inferences, likelihood value, and other output.
* `/models/lam/lam_ssf.m`:Takes the 7 parameters from the Lam model of output and places them into state-space form as indicated in Kim and Nelson (1999, pgs. 111-112).
* `/models/lam/lam_constraints.m`: Transforms the unconstrained parameters of the Lam model to constrained parameters. This function is intended to be called by `ssm_objective`. The purpose of this function is so that unconstrained optimization techniques may be used.
