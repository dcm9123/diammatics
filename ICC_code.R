#Daniel Castaneda Mogollon, PhD
#March 22nd, 2025
#This script was made in order to determine the ICC of the genome consortia in the
#T1D project. I used the TP, FP matrix from the 16S ASV vs the genomes I identified
#to determine the sensitivity and PPV of the matrix and then see the level of
#agreement between the four consortia.

install.packages("psych")
library(psych)

sensitivity_matrix_species = matrix(c(
  0.64, 0.64, 0.64,
  0.65, 0.58, 0.62,
  0.64, 0.58, 0.63,
  0.58, 0.58, 0.54), ncol=3, byrow=TRUE)

sensitivity_matrix_genus = matrix(c(
  0.95, 0.95, 0.95,
  0.94, 0.94, 0.94,
  0.95, 0.95, 0.9,
  0.95, 1, 1), ncol=3, byrow=TRUE)

sensitivity_matrix_species_t1d = matrix(c(
  0.72, 0.76, 0.72,
  0.81, 0.77, 0.81,
  0.79, 0.79, 0.83,
  0.88, 0.75, 0.75), ncol=3, byrow=TRUE)

icc_sensitivity_species = ICC(sensitivity_matrix_species)
icc_sensitivity_species
icc_sensitivity_genus = ICC(sensitivity_matrix_genus)
icc_sensitivity_genus
icc_sensitivity_species_t1d = ICC(sensitivity_matrix_species_t1d)
icc_sensitivity_species_t1d
