Install using devtools::install_github("alasorous/futilityOC")

This program evaluates the operating characteristics (OCs) of a futility analysis in clinical trials
in the setting of model-based drug development (MBDD) following the rules proposed by Lalonde (2007).
With hypothesized treatment effects, the OCs are evaluated as the chance of a trial being terminated earlier
at the interim analysis if no interested treatment effect exists and the probabilities of Go/NoGo decisions
with a trial designed without a futility analysis versus the probabilities of those decisions based on a trial
with a futility analysis. The futility criterion is specified in the R function as the probability that the trial
would be erroneously terminated at the futility analysis when the true effect equals the value hypothesized under H1.
Commonly used futility scales (i.e. conditional power and predictive power) are provided with their thresholds equivalent
to the probability under H1 in the outcome of the function. Many other function arguments can be specified to closely
reflect the trial design, hypotheses and statistical information fraction available at the futility analysis.
For convenience of use, default values are set so that users only need to input or modify a few necessary
parameter values in order to use the function.
