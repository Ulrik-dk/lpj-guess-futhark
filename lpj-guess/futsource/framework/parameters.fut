-------------------------------------parameters.h--------------------------------------
open import "../futhark-extras"

-- Enums needed by some of the global instruction file parameters defined below

--/ Vegetation 'mode', i.e. what each Individual object represents
-- Can be one of:
--  1. The average characteristics of all individuals comprising a PFT
--     population over the modelled area (standard LPJ mode)
--  2. A cohort of individuals of a PFT that are roughly the same age
--  3. An individual plant
--/
type vegmodetype = #NOVEGMODE | #INDIVIDUAL | #COHORT | #POPULATION

--/ Land cover type of a stand. NLANDCOVERTYPES keeps count of number of items.
--  NB. set_lc_change_array() must be modified when adding new land cover types
--/
type landcovertype = #URBAN | #CROPLAND | #PASTURE | #FOREST | #NATURAL | #PEATLAND | #BARREN
let NLANDCOVERTYPES : int = 7 -- the c++ code had this as part of the enumerator - hack!

--/ Water uptake parameterisations
-- \see water_uptake in canexch.cpp
type wateruptaketype = #WR_WCONT | #WR_ROOTDIST | #WR_SMART | #WR_SPECIESSPECIFIC

--/bvoc: define monoterpene species used
type monoterpenecompoundtype = #APIN | #BPIN | #LIMO | #MYRC | #SABI | #CAMP | #TRIC | #TBOC | #OTHR
let NMTCOMPOUNDTYPES : int = 9 -- the c++ code had this as part of the enumerator - hack!

--/ Fire model setting. Either use
--	One of
--	BLAZE 		Use the BLAZE model to generate fire fluxes
--                      (must be accompanied by ignitionmode DEFAULT)
--	GLOBFIRM	fire parameterization following Thonicke et al. 2001
--	NOFIRE		no fire model
--/
type firemodeltype = #BLAZE | #GLOBFIRM | #NOFIRE

--/ Type of weathergenerator used
--     One of:
--      GWGEN           Global Weather GENerator (needed by BLAZE, due to
--                      additional rel. humidity and wind DEFAULT)
--      INTERP          use standard interpolation scheme
--      NONE            Should be set if daily input is used (e.g. in cfinput)
--/
type weathergeneratortype = #GWGEN | #INTERP | #NONE

--/How to determine root distribution in soil layers
type rootdisttype = #ROOTDIST_FIXED | #ROOTDIST_JACKSON

--- The global Paramlist object
--- C---------------------------------------------------------------------------------------
-- Global instruction file parameters

--- Title for this run
let title : xtring = 1

--- Vegetation mode (population, cohort or individual)
let vegmode : vegmodetype = #POPULATION --TODO

--- Default number of patches in each stand
--- Should always be 1 in population mode,
--  cropland stands always have 1 patch.
--  Actual patch number for stand objects may differ and
--  should always be queried by stand.npatch()

let npatch : int = 1

--- Number of patches in each stand for secondary stands
let npatch_secondarystand : int = 1

--- Whether to reduce equal percentage of all stands of a stand type at land cover change
let reduce_all_stands : int = 1

--- Minimum age of stands to reduce at land cover change
let age_limit_reduce : int = 1

--- Patch area (m2) (individual and cohort mode only)
let patcharea : real = 1.0

--- Whether background establishment enabled (individual, cohort mode)
let ifbgestab : bool = true

--- Whether spatial mass effect enabled for establishment (individual, cohort mode)
let ifsme : bool = true

--- Whether establishment stochastic (individual, cohort mode)
let ifstochestab : bool = true

--- Whether mortality stochastic (individual, cohort mode)
let ifstochmort : bool = true

--- Fire-model switch
let firemodel : firemodeltype = #BLAZE

--- Weather Generator switch
let weathergenerator : weathergeneratortype = #INTERP

--- Whether "generic" patch-destroying disturbance enabled (individual, cohort mode)
let ifdisturb : bool = true

--- Generic patch-destroying disturbance interval (individual, cohort mode)
let  distinterval : real = 1.0

--- Whether SLA calculated from leaf longevity (alt: prescribed)
let ifcalcsla : bool = true

--- Whether leaf C:N ratio minimum calculated from leaf longevity (alt: prescribed)
let ifcalccton : bool = true

--- Establishment interval in cohort mode (years)
let estinterval : int = 1

--- Whether C debt (storage between years) permitted
let ifcdebt : bool = true

--- Water uptake parameterisation
let wateruptake : wateruptaketype = #WR_SMART

--- Parameterisation of root distribution
let rootdistribution : rootdisttype = #ROOTDIST_FIXED

--- whether CENTURY SOM dynamics (otherwise uses standard LPJ formalism)
let ifcentury : bool = true

--- whether plant growth limited by available N
let ifnlim : bool = true

--- number of years to allow spinup without nitrogen limitation
let freenyears : int = 1

--- fraction of nitrogen relocated by plants from roots and leaves
let  nrelocfrac : real = 1.0

--- first term in nitrogen fixation eqn (Cleveland et al 1999)
let  nfix_a : real = 1.0

--- second term in nitrogen fixation eqn (Cleveland et al 1999)
let  nfix_b : real = 1.0

--- whether to use nitrification/denitrification in CENTURY SOM dynamics
let ifntransform : bool = true
--- Fraction of microbial respiration assumed to produce DOC, 0.0,0.3
let  frac_labile_carbon : real = 1.0

--- Soil pH (used for calculating N-transformation), 3.5,8.5
let  pH_soil : real = 1.0
--- Maximum nitrification rate, 0.03,0.15
let  f_nitri_max : real = 1.0
--- Constant in denitrification, 0.001,0.1
let  k_N : real = 1.0
--- Constant in temperature function for denitrification, 0.005,0.05
let  k_C : real = 1.0
--- Maximum gaseus losses in nitrification
let  f_nitri_gas_max : real = 1.0
--- Maximum fraction of NO3 converted to NO2
let  f_denitri_max : real = 1.0
--- Maximum fraction of NO2 converted to gaseus N
let  f_denitri_gas_max : real = 1.0


---------------------------------------------------------------------------------------
-- Landuse and crop settings

--- Whether other landcovers than natural vegetation are simulated.
let run_landcover : bool = true

--- Whether a specific landcover type is simulated (URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND, BARREN).
let run : [NLANDCOVERTYPES]bool = replicate NLANDCOVERTYPES true

--- Whether landcover fractions are not read from input file.
let lcfrac_fixed : bool = true

--- Whether fractions of stand types of a specific land cover are not read from input file.
let frac_fixed : [NLANDCOVERTYPES]bool = replicate NLANDCOVERTYPES true

--- Set to false by initio( ) if fraction input files have yearly data.
let all_fracs_const : bool = true

--- If a slow harvested product pool is included in patchpft.
let ifslowharvestpool : bool = true

-- If grass is allowed to grow between crop growingseasons
let ifintercropgrass : bool = true

-- Whether to calculate dynamic potential heat units
let ifcalcdynamic_phu : bool = true

-- Whether to use gross land transfer: simulate gross lcc (1) read landcover transfer matrix input file (2) read stand type transfer matrix input file (3), or not (0)
let gross_land_transfer : int = 1

-- Whether gross land transfer input read for this gridcell
let gross_input_present : bool = true

-- Whether to use primary/secondary land transition info in landcover transfer input file (1). or not (0)
let ifprimary_lc_transfer : bool = true

-- Whether to use primary-to-secondary land transition info (within land cover type) in landcover transfer input file (1). or not (0)
let ifprimary_to_secondary_transfer : bool = true

-- Pooling level of land cover transitions 0: one big pool 1: land cover-level 2: stand type-level
let transfer_level : int = 1

-- Whether to create new stands in transfer_to_new_stand() according to the rules in copy_stand_type()
let iftransfer_to_new_stand : bool = true

-- Whether to limit dynamic phu calculation to a period specified by nyear_dyn_phu
let ifdyn_phu_limit : bool = true

-- Number of years to calculate dynamic phu if dynamic_phu_limit is true
let nyear_dyn_phu : int = 1

--- number of spinup years
let nyear_spinup : int = 1

--- Whether to use sowingdates from input file
let readsowingdates : bool = true

--- Whether to use harvestdates from input file
let readharvestdates : bool = true

--- Whether to read N fertilization from input file
let readNfert : bool = true

--- Whether to read manure N fertilization from input file
let readNman : bool = true

--- Whether to read N fertilization (stand tyoe level) from input file
let readNfert_st : bool = true

--- Whether to print multiple stands within a land cover type (except cropland) separately
let printseparatestands : bool = true

--- Whether to simulate tillage by increasing soil respiration
let iftillage : bool = true

--- Use silt/sand fractions per soiltype
let textured_soil : bool = true

--- Whether pastures are affected by disturbance and fire (affects pastures' npatch)
let disturb_pasture : bool = true

--- Whether to simulate cropland as pasture
let grassforcrop : bool = true

---------------------------------------------------------------------------------------
-- Settings controlling the saving and loading from state files

--- Location of state files
let state_path : xtring = 1

--- Whether to restart from state files
let restart : bool = true

--- Whether to save state files
let save_state : bool = true

--- Save/restart year
let state_year : int = 1

--- The level of verbosity
let verbosity : int = 1

--- whether to vary mort_greff smoothly with growth efficiency (1) or to use the standard step-function (0)
let ifsmoothgreffmort : bool = true

--- whether establishment is limited by growing season drought
let ifdroughtlimitedestab : bool = true

--- rain on wet days only (1, true), or a little every day (0, false)
let ifrainonwetdaysonly : bool = true

--- whether BVOC calculations are included
let ifbvoc : bool = true

---------------------------------------------------------------------------------------
-- Arctic and wetland inputs

--- Use the original LPJ-GUESS v4 soil scheme, or not. If true, override many of the switches below.
let iftwolayersoil : bool = true

--- Use multilayer snow scheme, or the original LPJ-GUESS v4 scheme
let ifmultilayersnow : bool = true

--- whether to reduce GPP if there's inundation (1), or not (0)
let ifinundationstress : bool = true

--- Whether to limit soilC decomposition below 0 degC in upland soils (1), or not (0)
let ifcarbonfreeze : bool = true

--- Extra daily water input or output, in mm, to wetlands. Positive values are run ON, negative run OFF.
let wetland_runon : real = 1.0

--- Whether methane calculations are included
let ifmethane : bool = true

--- Whether soil C pool input is used to update soil properties
let iforganicsoilproperties : bool = true

--- Whether to take water from runoff to saturate low latitide wetlands
let ifsaturatewetlands : bool = true
