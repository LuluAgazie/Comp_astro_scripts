# Comp_astro_scripts

#Scripts written for my Computational Astrophysics course during the Fall 2019 semester at West Virginia University. 

#Some of the scripts require a txt file to be in the same directory as the script and the ones that I have used are included.

#Acronyms: HI- Neutral Hydrogen; HIMF-Neutral Hydrogen Mass Function; FRB - Fast Radio Burst; MCMC- Monte Carlo Markov Chain; SNR-Signal to noise ratio

#####Short descriptions of each code#####

##FRB_dist.py : A code that runs a 2 sample KS test to determine whether the current population (as of Nov 29, 2019) of FRBs has a number density vs distance relation that follows the relation predicted by the star formation rate equation. Requires the file frbcat_data.txt to run

##FRB_pop.py : A code that runs a linear regression analysis to determine whether the flux density and SNR of known FRBs has any relation. Calculates the Pearson correlation coefficient, population correlation coefficient, Spearman rank coefficient, and Fischer transform. Requires the file frbcat.txt to run

##HIMF.py : A code that simulates a random population of galaxies and uses the rejection method to determine which galaxies are detectable, and calculates the HIMF for each galaxy. 

##HIMF_bayes.py : This code takes the simulated galaxies from HIMF.py and uses bayes theorem to recover one of the parameters (alpha) of the Schechter function used to define the HIMF.

##HIMF_mcmc.py : This code takes the simulated galaxies from HIMF.py and uses an MCMC to recover three parameters (alpha, local galactic density, and charactistic HI mass) of the Schechter function.

##UGC_328_HI.py : This code uses HI spectral line data on galaxy UGC328 taken on the Green Bank Telescope to calculate HI mass and the dynamical mass of UGC 328. Requires UGC328.dat to run

##asteroid_sim.py: This code simulates asteroids (with various distances from the sun) within the asteroid belt that have the gravitational force of the Sun and Jupiter acting on them. It simulates the orbit of these asteroids over 10000 years in order to demonstrate the asteroid orbit instability at the Kirkwood gaps.  
