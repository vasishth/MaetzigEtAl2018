# Next steps ICCM paper
#iccm17 #phd  #msc #masterarbeit

- [x] check for stability of ANS with multiple 1000 iterations
		* Measures to calculate for this check: failure rate, average reading times (TFT FFD e.g.), average accuracy
		* write script for this, copy parts from Felix’ script
- [x] adapt python script (acc calculation) so that it  works for SR + OR sentences simultaneously
- [x] adapt python script (acc calculation) with random guess in case of parsing failure
- [x] use analysis scripts on the guess data and look at the discrepancies between strict failure and guessing
- [ ] make plan of implementation of reflexive rule in ACT-R
- [ ] plan simulation on comparing ANS vs GA
		* hypothesis [Felix]: ANS and GA have very similar effects with more misretrievals and lower differences in retrieval latency >> effectively, lower GA makes all chunks less different to the retrieval target; higher ANS can result in the same, but with a more random element
		* The crucial difference in GA vs ANS, conceptually, is the random element in ANS (which is very transparent there), and which is much less transparent in GA
		* What is the random element in GA?
		* a lot of data would be beneficial to look at _correlations_
- [ ] plan simulation on comparing ANS and UN (utility noise)

## Adapted python script (random guess when parsing failure)
* original script always interpreted parsing failure (no message of created dependency in “n-relations.txt”) as “wrong answer / 0”
* this is a conservative measure, and it’s more a model of actual dependency creation in ACT-R, than a model of the response
* my idea was: an actual human would, if he had no idea of the dependency, randomly guess between two options 
* the result is an overly positive estimate of accuracy: 

## Stability of simulations (esp. ANS manipulation)
* taking the averaged accuracy (across SR and OR, which is a stupid measure) the accuracy stayed stable across 3 1000-iteration-simulations
* except RRT, all measures are very stable across the simulations with high ANS value

## ANS vs. UN simulation
* We expect that modelling ANS and GA at the same time does not yield a better model than substituting ANS for UN. Why?
	* different manipulations to the model, but very similar effects:
		* more misretrievals, lower differences in retrieval latency
* however, GA is much less random than ANS
* if we look at the results, GA seems to be the most variable parameter in the ICCM paper simulations
	* Would the same happen if we varied GA-DAT-UN simultaneously?
		* Expectation: **no**
	* What we need is correlation data, but for that we need more data

## Reflexive rule in ACT-R (for SSP/SOP and SSRP/SORP sentence types)
* SR / OR: _The girl who hugged the boy washed the woman._
* SRPRO / ORPRO: _The woman who hugged the girl_i washed her_i._
* SRREF / ORREF: _The woman_j who hugged the girl washed herself_j._
* attach pronoun as object of matrix V
* 

## Readings: mapping hypotheses to ACT-R parameters
* in my opinion, in the big paper this deserves more attention; we should explain clearly why we decide for one specific parameter to represent the implementation of a specific hypothesis on the aphasic SP deficit
* a lot of good stuff is in the Patil et AL (2016) paper in the literature review
* it would make sense to include discussion 