Pyrolysis Kinetic Preprocessor
================
**PKP** (Pyrolysis Kinetic Preprocessor) is software for preprocessing
kinetics of coal and biomass pyrolysis.

----

## Introduction

**PKP** was developed at TU Bergakademie Freiberg and Universitaet Stuttgart. It interfaces different *detailed pyrolysis software*:

* CPD
* FG-DVC
* Flashchain
* Polimi

These softwares can be used for calibrating empirical pyrolysis models, which can be easily included and used for CFD simulations.
The following empirical models are available:

* SFOM (Single First Order reaction Model)
* C2SM (Competing 2 Step Model) 
* DAEM (Distributed Activation Energy Model)

More information about PKP can be found in [Vascellari et al.]

## Installation

```bash
# python setup.py install --user
```

## Use

```bash
pkp-cli

Usage:
      pkp-cli -h | --help
      pkp-cli generate (--json-input=<string> | --file-input=<loc>) --results-folder=<loc>
      pkp-cli generate-only (--json-input=<string> | --file-input=<loc>) --results-folder=<loc>
      pkp-cli fit-only      (--json-input=<string> | --file-input=<loc>) --fit-target=<string> --results-folder=<loc>

Options:
    -h  --help              Shows this screen
    --json-input=<string>   Foo
    --file-input=<loc>      Bar
    --results-folder=<loc>  Baz
    --fit-target=<string>   Bla
```

### Input file

Input files are written using **YAML** syntax:

```yaml
#Proximate Analysis (in percent, as recieved):
Coal:
    Proximate Analysis:
        Fixed Carbon    : 36.12
        Volatile Matter : 45.26
        Moisture        : 12.92
        Ash             : 5.70
    Ultimate Analysis:
        Carbon   : 68.342
        Hydrogen : 4.928
        Nitrogen : 0.695
        Oxygen   : 25.579
        Sulphur  : 0.456
    #Higher Heating Value, as received, in J/kg:
    hhv :   0.0 
    #Tar Molecule weight, MTar:
    rhoDry : 1310 # kg/m3

FIT:
    #Model: ['constantRate', 'arrheniusRate'] #, 'ArrheniusNoB', 'Kobayashi', 'DAEM', 'None']
    Model: arrheniusRate
    #Model: constantRate
    # [ftot, fh20, fco, fco2, fch4, ftar, fgas, fother]
    Species: [ftot] # empty list for all availible species
    constantRate:
        k :  1.0
        kBounds :  [0.0, 1.0]
        tstart : 0.0
        tstartBounds : [0.0, 1.0]
        finalYield : 1.0
        finalYieldBounds: [0.0, 1.0]
    arrheniusRate:
        preExp: 0.0
        preExpBounds: [0, 10000.0]
        beta : 0.0
        activationEnergy : 0
        activationEnergyBounds : [100.0, 100000.0]
        lowerDevolTemp : 600.0

# TODO MV: names of empirical models should be modified for being more coherent with publication in Fuel 2013 
# Arrhenius -> SFOR (single first order reaction)
# Kobayashi -> C2SM (Competing 2 steps model)
# DAEM
# constantRate -> SCR (single constant rate)
CPD:
    active : true
    deltaT : 1e-4
    MW_TAR : 130

FGDVC:
    active      : false
    fit         : None #?
    dir         : C:\FGDVC_8-2-3\
    dirOut      : C:\FGDVC_8-2-3\FGDVC\
    #Choose Coal: 0 interpolate between library coals and generate own coal. Set 1 to 8 for a library coal.
    refCoal     : 0 
    #Model tar cracking? If no, set tar residence time equal 0. For a partial tar cracking enter the tar residence time in s. For full tar cracking write -1.
    tarModel    : 0 
    # time step
    timeStep    : 1e-4
    
# I rename PCCL to FLASHCHAIN
FLASHCHAIN:
    active      :   false
    fit         :   None #?
    # PC Coal Lab Main Path:
    dir         :   C:\Users\vascella\Documents\PCCL\
    # PC Coal Lab executable name:
    execName    :   PCCoV41M1to7Par.exe
    # Known PC Coal Lab Coal Calibration Factor? None or Value(float): ?
    calibrationFactor   : None
    # Particle Size in micrometer:
    diameter    :   100.

PMSKD:  # Polimi Multi-Step Kinetic Devolatilization
    active      :   false
    fit         :   None #?
    # number of step
    nStep       :   100
    # Mechanism file for Cantera (only use xml!!!)
    mechFile    :   COAL.xml

Calibration:
    weightYield     :   100.
    weightRate      :   0.
    # species to fit: total, mainSpecies, allSpecies
    speciesToFit    :   total 

OperatingConditions:
    pressure    : 1.0 #atmosphere
    runs: 1
    run0        : [ [ 0, 400], [ 0.5, 2000] ]
    run1        : [ [ 0, 400], [ 0.1, 1400], [ 0.5, 2000] ]
    run2        : [ [ 0, 400], [ 0.1, 1400], [ 0.5, 2000] ]
    run3        : [ [ 0, 400], [ 0.1, 1400], [ 0.5, 2000] ]
    run4        : [ [ 0, 400], [ 0.1, 1400], [ 0.5, 2000] ]
    run5        : [ [ 0, 400], [ 0.1, 1400], [ 0.5, 2000] ]
```


[Vascellari et al.]: 10.1016/j.fuel.2013.06.014 "Vascellari M, Arora R, Pollack M, Hasse C. Simulation of entrained flow gasification with advanced coal conversion submodels. Part 1: Pyrolysis. Fuel 2013;113:654–669."

