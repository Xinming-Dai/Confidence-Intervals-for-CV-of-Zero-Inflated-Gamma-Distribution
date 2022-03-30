# Confidence-Intervals-for-CV-of-Zero-Inflated-Gamma-Distribution

Code for "Confidence Interval for the difference of Coefficients of Variation of Two Zero-Inflated Gamma Distributions"

1. Run "1_Bootstrap Method Code.R", 
	"2_Fiducial Method Code.R", 
	"3_MOVER Method Code_ver1.R", and 
	"4_MOVER Method Code_ver2.R" first， 
	in order to run ci.pb, ci.f, ci.mover1, and ci.mover function. 


2. These four functions can construct confidence intervals for the difference of coefficient of variance of two zero-inflated gamma distribution using bootstrap, Fiducial, MOVER1, and MOVER2 methods respectively.

3. Then, simu function in Simulation.R can be used to calculate coverage probabilities, average lengths, lower errors, and upper errors of the four methods in the mean time.

4. Simulation time.R is for calculate time spent of the four methods.
