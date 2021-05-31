
ABOUT

	This code was written by Michael S. Crawford, for the journal article:

	Crawford, M., K. Barry, A. Clark, C. Farrior, J. Hines, E. Ladouceur, J. Lichstein, I. Marechaux, 
	F. May, B. Reineking, L. Turnbull, C. Wirth, and N. Rüger. 2021. The function-dominance correlation drives 
	the direction and strength of biodiversity-ecosystem functioning relationships. Ecology Letters.


INSTRUCTIONS

	To reproduce the analysis associated with the Crawford, Berry et al. (2021), source the
	runscript (R/SimNet_Runscript.R), within RStudio or hard code the
	'base_dir' to the outermost directory of this analysis. You will need
	to ensure that you have all the required packages installed on your machine
	as well (see the separate analysis scripts for the exact packages.)

	Because the analysis can take a while, if you're running it repeatedly and 
	would like to save your intermediate steps, you may want to use the rudimentary
	caching functionality (see Runscript.R).

	All output is placed within the tmp folder. The final datasets - including
	the derived `brms` models used in our analysis - are contained in the data folder. 


LICENSE

	We invite others to use this analysis and data, allowing for its proper citation. 
	All code herein is subject to the Creative Commons Attribution 4.0 International license.
