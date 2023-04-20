from hapi import *

##Look at https://hitran.org/docs/iso-meta/ for the molecular isotopologue

#fetch('data/H2O',1,1, 1000, 10000)#Download from 1 microns to 10 micron range. 1000 per cm = 10 microns and 10000 per cm is 1 microns
fetch('data/CO',5,1,1000, 10000)