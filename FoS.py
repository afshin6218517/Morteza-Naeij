# -*- coding: mbcs -*-
#### PLEASE DO NOT CHANGE THE FIRST LINE OF THIS CODE ####
# ***********************************************************************************************************   
# ** This code was prepared for calculating FoS of geotechnical slopes based on strength reduction procedure*
# ** This code could calculate FoS>1 only and for the case of FoS<=1 some modifications should be done      *
# ** This code is applicable with Mohr-Coulomb VUMAT and for soil geometry with one or two layers           *
# ** This code calculate FoS in 0.01 intervals and for having an efficient simulation runtime applying      *
# ** sophesticaed algorithems is useful                                                                     *
# ** The final result (minimum FoS) will be written in FoS.csv file                                         * 
# ** In order to use PYTHON file, FORTRAN VUMAT, and ABAQUS .cae file, please follow the instruction file   *
# ***********************************************************************************************************
# ***********************************************************************************************************
# ** This code was written by:                                                                              *
# ** Morteza Naeij (morteza.naeij@ut.ac.ir)                                                                 *
# ** Hussein Ghasemi                                                                                        *
# ** Danial Ghaffarian                                                                                      *
# ** Yousef Javanmardi (Y.javanmardi@ucl.ac.uk)                                                             *
# ** and is a part of supplementary data for manuscript with title:                                         *
# ** A novel procedure for finite element method to calculate the minimum factor of safety                  *
# ** 
# ***********************************************************************************************************  
# ** Please make sure that the ABAQUS file with name "slope.cae" is available in ABAQUS repository file,    *
# ** usually c\temp\                                                                                        *
# ***********************************************************************************************************
# ** Apply initial material parameters in this file (rad_fi and ci for first layer and rad_fi1 and ci for   *
# **  second layer)                                                                                         *
# *********************************************************************************************************** 
from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import *
import sys
import math
from array import *
openMdb('slope.cae')
for rfs in range (0,250) :
    xx=rfs/100.0
    rad_fi=25*pi/180.0 # down
    rad_fi1=30*pi/180.0 # mid
    fi=(math.atan(math.tan(rad_fi)/(1.0+xx)))*180/pi
    fi1=(math.atan(math.tan(rad_fi1)/(1.0+xx)))*180/pi
    ci=10000.0/(1.0+xx) # down
    ci1=25000.0/(1.0+xx) # mid
    mdb.models['Model-1'].materials['soil'].UserMaterial(mechanicalConstants=(
        140000000.0, 0.3, ci, fi, 1.0, 0.0, 0.2))
#    mdb.models['Model-1'].materials['soil1'].UserMaterial(mechanicalConstants=(
#        140000000.0, 0.3, ci1, fi1, 1.0, 0.0, 0.2))
    mdb.jobs['slope'].submit(consistencyChecking=OFF)
    mdb.jobs['slope'].waitForCompletion()
# -*- coding: mbcs -*-
    import sys
    import csv
    from odbAccess import *
    odb = openOdb ('slope.odb')
#: Model: C:/Temp/slope.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       4
#: Number of Node Sets:          3
#: Number of Steps:              2
    step = odb.steps['Step-1']
    #region=step.historyRegions['ElementSet SOIL-1.SET-1']
    #print step.historyRegions
    region=step.historyRegions['ElementSet SOIL']
    # SOIL-1.SET-1
    results=region.historyOutputs['ALLKE'].data
    initial_energy = 0.0
    for nn in range (0,180):
        if results[nn][1] >= initial_energy:
            initial_energy=results[nn][1]
#   initial_energy = results[180][1]
    end_energy = results[199][1]
    time_ratio = end_energy/initial_energy
    print time_ratio
    FoS=1.0
    if  time_ratio >= 1.1 :
                        resultsFile = open('FoS.csv','w')
                        FoS=rfs+1.0
                        resultsFile.write('%10.4E\n'% (FoS))
                        resultsFile.close()
                        odb.save()
                        odb.close()
                        sys.exit()
    odb.save()
    odb.close()
    #os.remove('C:\Temp\slope.odb')