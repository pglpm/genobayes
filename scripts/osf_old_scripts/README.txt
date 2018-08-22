### README.txt 
## Author: Porta Mana, Bachmann
## Created: 2018-03-06T17:18:11+0100
## Last-Updated: 2018-03-21T20:02:19+0100

-- NB: Please note that these codes will be updated from time to time. --

The sampling described in sect 2.5.1 of the paper is made in R with three
scripts:

- definitions.R: contains definitions of the functions in the scripts. Its
  main functions are 'logevidence', which does the sampling, and 'genplot'
  which plots the results (fig. 3)


- sample1.R: here the various parameters for the sampling are given, for a
  particular model, and then the sampling function is called, and the final
  output saved in an RDS file, eg 'samples_logit.rds'. It needs files with
  patient data for each health condition, see below.

**NOTE: with 520000 samples the output file is around 1.6 GB. Size scales
  more or less linearly with sample size **


- plotresults.R: plots the averages from the samples saved in the RDS
  files, as in fig 3 in the paper


- Data files must have one column for each graph property (connectivity
  weight in the paper), one row for each patient

-weights_con_40cons
Data for healthy (control) patients:

- weights_Schizo_40cons
Data for schizophrenic patients:


- weights_keys_40cons
(not used) IDs of the connectivity weights (eg, "f[36, 91]" is the
connectivity between regions 36 and 91):
