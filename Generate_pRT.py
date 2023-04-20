#Code written by Prajwal Niraula
#email: prajwalniraula@gmail.com

import hapi
import os
import numpy as np
import itertools
import matplotlib.pyplot as plt

AvogadroNum = 6.0221409e23


#Specify the molecule
MoleculeName = "CO"

Wavelength = np.fromfile("pRT_Data/wlen.dat", dtype=np.float64)
WaveNumber = 1./Wavelength[::-1]


Pressure_Array = np.array([-6.00, -5.00, -4.00, -3.00, -2.00, -1.0,  0.0, 1.00,  2.00, 3.00])
Temp_Array = np.array([81, 110, 148, 200, 365, 493, 900, 1641, 2217, 2995])

hapi.db_begin('data')




#Now read the data from the molecule
AllData = hapi.LOCAL_TABLE_CACHE[MoleculeName]['data']
MoleculeNumberDB = hapi.LOCAL_TABLE_CACHE[MoleculeName]['data']['molec_id'][0]
IsoNumberDB = hapi.LOCAL_TABLE_CACHE[MoleculeName]['data']['local_iso_id'][0]
Molecularmass = hapi.molecularMass(MoleculeNumberDB,IsoNumberDB)
Abundance = hapi.abundance(MoleculeNumberDB, IsoNumberDB)
MolecularID = str(MoleculeNumberDB).zfill(2)
print("The molecular ID is given by", MolecularID)

if not(os.path.exists(MoleculeName)):
    os.system("mkdir %s" %MoleculeName)

with open(MoleculeName+"/molparam_id.txt", 'w') as f:
    Text2Write = "#### Species ID (A2) format\n%s\n#### molparam value\n%s" %(MolecularID, str(Abundance))
    f.write(Text2Write)

os.system("cp prT_Data/wlen.dat %s/" %MoleculeName)
 

for TValue, PValue in list(itertools.product(Temp_Array,Pressure_Array)):
    print("Chaing the molecular ID")
    PValue = "%7.6f" %(10**PValue)
    SaveName = MoleculeName+"/sigma_"+MolecularID+"_"+str(TValue)+".K"+"_"+PValue+"bar.dat"
    print("The Save Name is given by:", SaveName)
    

    #Now generate the cross-section
    #Get molecular ID
    #Converting bars to atmosphere
    nu_Voigt, coef_Voigt = hapi.absorptionCoefficient_Voigt(SourceTables=MoleculeName, Diluent={'air':1}, \
				Environment={'p':float(PValue)/0.986923, 'T':TValue}, OmegaGrid=WaveNumber, HITRAN_units=True)
    coef_Voigt = coef_Voigt[::-1]
   
    #Now properly normalize the cross-section
    NormalizationFactor = AvogadroNum/Molecularmass
    coef_Voigt*=NormalizationFactor
   

    pRT_CS = np.fromfile("/media/prajwal/b3feb060-a565-44ab-a81b-7dd59881cba0/petitRADTRANS/petitRADTRANS/input_data/opacities/lines/line_by_line/CO/sigma_05_81.K_0.000001bar.dat", dtype=np.float64)

    coef_Voigt.tofile(SaveName)
    
   