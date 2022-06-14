-------------------------------------soil.h-------------------------------------------
--------------------------------------------------------------------------------------/
--/ \file soil.h
--/ \brief Constants and parameters used in Arctic and wetland code, with references.
--/ NB: The class Soil and its member functions and variables are declared in guess.h,
--/ while its member functions are implemented in soil.cpp and in soilmethane.cpp.
--/
--/ \author Paul Miller
--/ $Date: 2019-03-10 16:34:14 +0100 (Sun, 10 Mar 2019) $
--/
--------------------------------------------------------------------------------------/

--#ifndef LPJ_GUESS_SOIL_H
--#define LPJ_GUESS_SOIL_H

-- DEBUGGING BOOLEANS

--let DEBUG_SOIL_WATER : real =  false
--let DEBUG_SOIL_TEMPERATURE : real =  false
--let DEBUG_METHANE : real =  false
open import "../framework/parameters"
open import "../futhark-extras"

-- Solving Century SOM pools

-- fraction of nyear_spinup minus freenyears at which to begin documenting for calculation of Century equilibrium
let SOLVESOMCENT_SPINBEGIN  : real = 0.1
-- fraction of nyear_spinup minus freenyears at which to end documentation and start calculation of Century equilibrium
let SOLVESOMCENT_SPINEND : real = 0.3

-- Kelvin to deg C conversion
let K2degC : real = 273.15

-- Maximum number of crop rotation items
let NROTATIONPERIODS_MAX : int = 3

-- Conversion factor for CO2 from ppmv to mole fraction
let CO2_CONV : real = 1.0e-6

-- Initial carbon allocated to crop organs at sowing, kg m-2
let CMASS_SEED : real = 0.01

-- Precision in land cover fraction input
let INPUT_PRECISION : real = 1.0e-14
let INPUT_ERROR : real = 0.5e-6
let INPUT_RESOLUTION : real = INPUT_PRECISION - INPUT_PRECISION * INPUT_ERROR

-- Averaging interval for average maximum annual fapar (SIMFIRE)
let AVG_INTERVAL_FAPAR : int = 3

-- Averaging interval for biome averaging (SIMFIRE)
let N_YEAR_BIOMEAVG : int = 3


-- Year at which to calculate equilibrium soil carbon
let SOLVESOM_END : int = 400

-- Year at which to begin documenting means for calculation of equilibrium soil carbon
let SOLVESOM_BEGIN : int = 350
-- SOIL DEPTH VALUES
-- soil upper layer depth (mm)
let SOILDEPTH_UPPER : real = 500.0
-- soil lower layer depth (mm)
let SOILDEPTH_LOWER : real = 1000.0

-- Depth of sublayer at top of upper soil layer, from which evaporation is
-- possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
-- Must be a multiple of Dz_soil
let SOILDEPTH_EVAP : real = 200.0

-- CONSTANTS
type insoltype = enum_type
let NOINSOL : insoltype = 0 -- No insolation type chosen
let SUNSHINE : insoltype = 1 -- Percentage sunshine
let NETSWRAD : insoltype = 2 -- Net shortwave radiation flux during daylight hours (W/m2)
let SWRAD : insoltype = 3 -- Total shortwave radiation flux during daylight hours (W/m2)
let NETSWRAD_TS : insoltype = 4 -- Net shortwave radiation flux during whole time step (W/m2)
let SWRAD_TS : insoltype = 5 -- Total shortwave radiation flux during whole time step (W/m2)
-- CENTURY pool names, NSOMPOOL number of SOM pools
type pooltype = enum_type
let SURFSTRUCT : pooltype = 0
let SOILSTRUCT : pooltype = 1
let SOILMICRO : pooltype = 2
let SURFHUMUS : pooltype = 3
let SURFMICRO : pooltype = 4
let SURFMETA : pooltype = 5
let SURFFWD : pooltype = 6
let SURFCWD : pooltype = 7
let SOILMETA : pooltype = 8
let SLOWSOM : pooltype = 9
let PASSIVESOM : pooltype = 10
let NSOMPOOL : pooltype = 11


--/ number of total soil layers. Must be at least NSOILLAYER + NLAYERS_SNOW + 1
let NLAYERS : int =  22

--/ Number of soil layers for soil temperature/water calculations. Typically 15 10mm layers, making up the 150cm-deep soil column
let NSOILLAYER : int =  15	-- rootdist in .ins file must have NSOILLAYER components

--/ Number of soil layers used in LPJ-GUESS v4.0 for soil water calculations. Typically 2 layers, 500mm + 10000mm, making up the 150cm-deep soil column
let NSOILLAYER_SIMPLE : int =  2

--/ number of depths at which we want the soil T output in outannual
let SOILTEMPOUT : int =  NSOILLAYER

--/ index of the first soil layer
let IDX_STD : int =  NLAYERS-NSOILLAYER

--/ number of padding layers in the soil
let PAD_LAYERS : int =  5

--/ number of total soil layers in the acrotelm
let NACROTELM : int =  3

--/ number of total soil layers in the catotelm
let NCATOTELM : int =  12

--/ Number of sublayers in the acrotelm, i.e. 1cm layers with this : real =  30
let NSUBLAYERS_ACRO : int =  30

--/ Maximum number of layers - for soil temperature calculations
let active_layersmax : int =  NLAYERS + PAD_LAYERS

--/ Total depth [mm] of padding layers - set to 8000 when analyticalSolutionTest : real =  true
let PAD_DEPTH : real =  48000

--/ Depth of each mineral soil layer [mm]
let Dz_soil : real =  100.0

--/ Depth of each acrotelm soil layer [mm]
let Dz_acro : real =  100.0

--/ Depth of each catotelm soil layer [mm]
let Dz_cato : real =  100.0

--/ Max height of standing water [mm]
let maxh : real =  0.0

--/ Slope of soil profile
let soil_slope : real =  -0.37

--/ Maximum density of soil organic carbon [KgC/m3]
let maxSOCdensity : real =  130.0 -- From Lawrence and Slater, 2008

-- SNOW PARAMETERS

--/ snow density at start of the snow season [kg m-3] (Wania et al. (2009a) have 150, Best et al. (2012) have 50)
let snowdens_start : real =  275.0

--/ snow density at end of the snow season [kg m-3] (Wania et al. 2009a)
let snowdens_end : real =  500.0

--/ maximum number of snow layers allowed (<= 5)
let NLAYERS_SNOW : int =  5

--/ ice density [kg m-3] - CLM value
let ice_density : real =  917.0

--/ water density [kg m-3]
let water_density : real =  1000.0

-- POROSITIES

--/ Porosity of organic material
let organic_porosity : real =  0.8 -- From Lawrence and Slater (2008) have 0.9, but 0.8 is consistent with organic soil code 8

--/ catotelm porosity
let catotelm_por : real =  0.92

--/ acrotelm porosity
let acrotelm_por : real =  0.98

--/ Gas fraction in peat
let Fgas : real =  0.00 -- Possible to reintroduce - was 0.08 in Wania et al (2010)

--/ Wilting point in peat
let peat_wp : real =  0.066

--/ First year when phase change is allowed
let FIRST_FREEZE_YEAR : real =  90

--/ time step [day]
let Dt : real =  1

-- HEAT CAPACITIES

--/ heat capacity of air [J m-3 K-1] - Bonan (2002)
let Cp_air : real =  1200

--/ heat capacity of water [J m-3 K-1] - Bonan (2002)
let Cp_water : real =  4180000

--/ heat capacity of ice [J m-3 K-1] : real =  2117.27 [J kg-1 K-1] * 917 [kg m-3 (ice density)] - CLM
let Cp_ice : real =  1941537

--/ heat capacity of organic matter [J m-3 K-1]
let Cp_org : real =  2500000

--/ heat capacity of dry peat (J m-3 K-1) - Bonan (2002), 0% water
let Cp_peat : real =  580000

--/ heat capacity of mineral soil [J m-3 K-1] - Bonan (2002)
let Cp_min : real =  2380000

--/ heat capacity of moss [J/ m-3 K-1] Ekici et al. (2015)
let Cp_moss : real =  2500000

-- NOTE: using the Cp_org values for Cp_peat and Korg for Kpeat does not
-- seem to influence the upper soil layer Ts, but increases the range
-- in the lower layers (2m)

-- THERMAL CONDUCTIVITIES
-- Values are from Hillel (1982) unless otherwise stated

--/ thermal conductivity of air [W m-1 K-1]
let Kair : real =  0.025

--/ thermal conductivity of water [W m-1 K-1]
let Kwater : real =  0.57

--/ thermal conductivity of ice [W m-1 K-1]
let Kice : real =  2.2

--/ thermal conductivity of organic matter [W m-1 K-1]
let Korg : real =  0.25

--/ thermal conductivity of dry peat [W m-1 K-1] - Bonan (2002)
let Kpeat : real =  0.06

--/ thermal conductivity of mineral soil [W m-1 K-1] - Wania et al. (2009a)
let Kmin : real =  2.0

--/ thermal conductivity of moss [W m-1 K-1]
let Kmoss : real =  0.25

--/ latent heat of fusion (Granberg et al. 1999) [J m-3]
let Lheat : real =  3.34E8

--/ A FILL-IN for layers above layer0.
let MISSING_VALUE : real =  -9999.0

--/ The number of subdaily timestep loops to perform.
-- Numerical instabilities can arise when there are thinner layers of snow and litter, so this should be > 1.
let TIMESTEPS : real =  2

--/ Thickness of the topmost air layer [m]
-- Note - keep air thickness low when using cnstep_full
let AIR_THICKNESS : real =  100.0


------------------------------------------------------------------------

-- SOIL CONSTANTS

--/ reduce infiltration and percolation rates sharply when there is ice present in the soil
-- See CLM4.5 - Swenson et al. (2012)
let ICE_IMPEDANCE: bool =  false

--/ min temp [deg C] for heterotrophic decomposition.
-- Clein & Schimel (1995) use -4 degC
let MIN_DECOMP_TEMP : real =  -8.0

--/ max CO2 [mimol/L] available to mosses in the acrotelm - see Wania et al (2009b)
-- value taken from Smolders et al. 2001 which give an average CO2 conc. of 70 sites as 934 mimol L-1
let PORE_WATER_CO2 : real =  934.0

-- METHANE CONSTANTS

--/ optimised moisture response under inundation (Wania et al. 2010, Table 5)
let RMOIST : real =  0.4

--/ Frolking et al (2001, 2010), Ise et al. (2008)
let RMOIST_ANAEROBIC : real=0.025

--/ time step for gas diffusion calculations [day]
let Dt_gas : real =  0.01

--/ CH4:CO2 ratio for peatland soils (> PEATLAND_WETLAND_LATITUDE_LIMIT N) - See Wania et al (2010)
-- Wania et al. (2010) optimal value: 0.1 (see Table 4).
-- McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.25, after optimisation
let CH4toCO2_peat : real =  0.085

--/ CH4:CO2 ratio for inundated soils (< PEATLAND_WETLAND_LATITUDE_LIMIT N) - See Spahni et al. (2011)
--let CH4toCO2_inundated : real =  0.024 -- SC1 value in Spahni et al. SC2 is 0.0415
let CH4toCO2_inundated : real =  0.027 -- Updated from Spahni et al. (2011) to match global emissions

--/ density of water [kg m-3]
let rho_H2O : real =  1000.0

--/ acceleration due to gravity [m s-2]
let gravity : real =  9.81

--/ Molecular mass of CH4 [g mol-1]
let mr_CH4 : real =  16.0

--/ Molecular mass of carbon [g mol-1]
let mr_C : real =  12.0

--/ molecular weight of water
let mr_h2o : real =  18.0

--/ universal gas constant [J mol-1 K-1]
let R_gas : real =  8.314472

--/ standard atmospheric pressure [Pa]
let atm_press : real =  101325.0

--/ coefficient for the calculation of the gas transport velocity, given in Riera et al. 1999
let n_coeff : real =  -0.5

--/ wind speed at 10m height [m s-1]
let U10 : real =  0.0

--/  Henry's Law constants [L atm mol-1] at 298.15K. Wania et al. (2010), Table 2
let henry_k_CO2 : real =  29.41
let henry_k_CH4 : real =  714.29
let henry_k_O2 : real =  769.23

--/ Constants [K] for CO2, CH4 and O2 for calculation of Henry's coefficient cited by Sander (1999). Wania et al. (2010), Table 2
let henry_C_CO2 : real =  2400.0
let henry_C_CH4 : real =  1600.0
let henry_C_O2 : real =  1500.0

--/ partial pressure of CH4 above water
let pp_CH4 : real =  1.7 -- micro atm

--/ partial pressure of O2 above water (value consistent with PO2 in canexch.h)
let pp_O2 : real =  209000 -- micro atm

--/ when ebullition occurs, the volumetric gas content (VGC) will drop to this level [unitless]
let vgc_low : real =  0.145

--/	ebullition occurs, when the volumetric gas content (VGC) exceeds this level [unitless, m3/m3]
let vgc_high : real =  0.15

--/ CH4 fraction of gas bubbles [unitless]
let bubble_CH4_frac : real =  0.57

--/ Fraction of oxygen used to oxidise CH4
-- Wania et al. (2010) optimal value: 0.5 (see Table 5).
-- McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.9, after optimisation
let oxid_frac : real =  0.5

--/ Fraction of ANPP used to calculate number of tillers
let ag_frac : real =  0.4

--/ Radius of an average tiller [m]
-- (tiller_radius : real =  0.004)  ! Schimel (1995) - Average over E. angustifolium (diam=7.9mm)
-- and C. aquatilis (diam=3.8mm)
-- Wania et al. (2010) optimal value: 0.003mm (see Table 5).
-- McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.0035, after optimisation
let tiller_radius : real =  0.0035


--/ Tiller porosity
-- (tiller_por : real =  0.6) ! Wetland plants book, eds. Cronk and Fennessy, p.90, values for 2 Erioph. spp.
let tiller_por : real =  0.7 -- Wania et al. (2010) optimal value: 0.7 (see Table 5)

--/ C content of biomass
let c_content : real =  0.45

--/ atomic mass of carbon [g/mol]
let atomiccmass : real =  12.0

--/ Individual tiller weight [g C]
let tiller_weight : real =  0.22

--/ a threshold factor for a minimum water content in the layer [unitless]
let water_min : real =  0.1

--/ Latitude (N). North of this and PEATLAND stands are treated as peatland as in Wania et al. (2009a, 2009b, 2010)
-- But south of this, then the PEATLAND stands are irrigated to avoid water stress, and can be a source of methane
let PEATLAND_WETLAND_LATITUDE_LIMIT : real =  40.0

--- Soiltype stores static parameters for soils and the snow pack.
--- One Soiltype object is defined for each Gridcell. State variables for soils
--  are held by objects of class Soil, of which there is one for each patch
--  (see below).



--- CENTURY SOIL POOL
type Sompool = {
  --- C mass in pool kgC/m2
  cmass : real,
  --- Nitrogen mass in pool kgN/m2
  nmass : real,
  --- (potential) decrease in C following decomposition today (kgC/m2)
  cdec : real,
  --- (potential) decrease in nitrogen following decomposition today (kgN/m2)
  ndec : real,
  --- daily change in carbon and nitrogen
  delta_cmass : real,
  delta_nmass : real,
  --- lignin fractions
  ligcfrac : real,
  --- fraction of pool remaining after decomposition
  fracremain : real,
  --- nitrogen to carbon ratio
  ntoc : real,

  -- Fire
  --- soil litter moisture flammability threshold (fraction of AWC)
  litterme : real,
  --- soil litter fire resistance (0-1)
  fireresist : real,

  -- Fast SOM spinup variables
  --- monthly mean fraction of carbon pool remaining after decomposition
  mfracremain_mean : [12]real
}

--- This struct contains litter for solving Century SOM pools.
-- We dont need getters, just read the field.
type LitterSolveSOM = {
  --- Carbon litter
  clitter : [NSOMPOOL]real,
  --- Nitrogen litter
  nlitter : [NSOMPOOL]real
}

type Soiltype = {
  soiltype_id : int,

  awc : [NSOILLAYER]real,
  --- fixed available water holding capacity of the standard Gerten soil layers [0=upper, 1 = lower] (mm)
  gawc : [2]real,

  --- coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996)
  perc_base : real,
  --- coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996) for the evaporation soil layer
  perc_base_evap : real,
  --- exponent in percolation calculation (=4 in Eqn 31, Haxeltine & Prentice 1996)
  perc_exp : real,

  --- thermal diffusivity at 0% WHC (mm2/s)
  thermdiff_0 : real,
  --- thermal diffusivity at 15% WHC (mm2/s)
  thermdiff_15 : real,
  --- thermal diffusivity at 100% WHC (mm2/s)
  thermdiff_100 : real,

  --- wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
  wp : [NSOILLAYER]real,
  --- saturation point. Cosby et al 1984
  wsats : [NSOILLAYER]real,

  --- organic soil fraction
  org_frac_gridcell : [NSOILLAYER]real,

  --- mineral soil fraction
  min_frac_gridcell : [NSOILLAYER]real,

  --- porosity of the soil
  porosity_gridcell : [NSOILLAYER]real,

  --- wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
  -- equivalents for the standard Gerten soil layers [0=upper, 1 = lower]
  gwp : [2]real,
  --- saturation point. Cosby et al 1984
  gwsats : [2]real,

  --- year at which to calculate equilibrium soil carbon
  solvesom_end : int,
  --- year at which to begin documenting means for calculation of equilibrium soil carbon
  solvesom_begin : int,

  --- water holding capacity plus wilting point for whole soil volume
  wtot : real,

  -- Sand, silt and clay fractions, should always add up to 1.
  --- fraction of soil that is sand
  sand_frac : real,
  --- fraction of soil that is clay
  clay_frac : real,
  --- fraction of soil that is silt
  silt_frac : real,
  --- pH
  pH : real,

  --- soilcode, 0 to 8
  soilcode : int,
  --- volumetric fraction of organic material (m3 m-3) (Hillel, 1998)
  organic_frac : real,
  -- water held below wilting point, important for heat conductance
  -- From the AGRMET Handbook, 2002
  water_below_wp : real,
  --- porosity from AGRMET Handbook, 2002
  porosity : real,
  -- volumetric fraction of mineral material (m3 m-3) (Hillel, 1998)
  -- = 1 - organic_frac - porosity
  mineral_frac : real,
  -- [cm], 3 depths at which monthly soil temperature is saved and output.
  -- Defaults: 25, 75 and 150 cm (if array values == 0)
  soiltempdepths : [10]real,
  -- Run on [mm/day]
  -- Currently only used for wetlands, but could be used for irrigation too
  runon : real,


  -- PEAT PROPERTIES (needed if there is a peat stand in this gridcell)

  --- fixed available water holding capacity of soil layers [0=upper layer] (mm)
  awc_peat : [NSOILLAYER]real,
  --- fixed available water holding capacity of the standard Gerten soil layers [0=upper, 1 = lower] (mm)
  gawc_peat : [2]real,
  --- wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
  wp_peat : [NSOILLAYER]real,
  --- saturation point. Cosby et al 1984
  wsats_peat : [NSOILLAYER]real,
  -- equivalents for the standard Gerten soil layers [0=upper, 1 = lower]
  --- wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
  gwp_peat : [2]real,
  --- saturation point. Cosby et al 1984
  gwsats_peat : [2]real,
  --- water holding capacity plus wilting point for whole soil volume
  wtot_peat : real,
  --- fraction of soil that is sand
  sand_frac_peat : real,
  --- fraction of soil that is clay
  clay_frac_peat : real,
  --- fraction of soil that is silt
  silt_frac_peat : real
}

--- Soil stores state variables for soils and the snow pack.
--- Initialised by a call to initdrivers. One Soil object is defined for each patch.
--  A reference to the parent Patch object (defined below) is included as a member
--  variable. Soil static parameters are stored as objects of class Soiltype, of which
--  there is one for each grid cell. A reference to the Soiltype object holding the
--  static parameters for this soil is included as a member variable.
--
--  NB: The class Soil and its member functions and variables are declared in guess.h,
--      while its member functions are implemented in soil.cpp and in soilmethane.cpp.
---
type Soil = {
  --- soil temperature today at 0.25 m depth (deg C)
  temp25 : real,
  --- water content of soil layers [0=upper layer] as fraction of available water holding capacity,
  wcont : [NSOILLAYER]real,
  --- water content of sublayer of upper soil layer for which evaporation from the bare soil surface is possible
  --- fraction of available water holding capacity--
  wcont_evap : real,

  --- reference to parent Patch object
  --patch : Patch, -- TODO circular definition
  --- reference to Soiltype object holding static parameters for this soil
  soiltype_id : int,
  --- the average wcont over the growing season, for each of the upper soil layers. Used in drought limited establishment.
  awcont_upper : real,
  --- daily water content in upper soil layer for each day of year
  dwcontupper : [Date_MAX_YEAR_LENGTH]real,
  --- mean water content in upper soil layer for last month
  --- (valid only on last day of month following call to daily_accounting_patch)--
  mwcontupper : real,
  --- stored snow as average over modelled area (mm rainfall equivalent)
  snowpack : real,
  --- total runoff today (mm/day)
  runoff : real,
  --- daily temperatures for the last month (deg C)
  --- (valid only on last day of month following call to daily_accounting_patch)--
  dtemp : [31]real,
  --- mean soil temperature for the last month (deg C)
  --- (valid only on last day of month following call to daily_accounting_patch)--
  mtemp : real,
  --- respiration response to today's soil temperature at 0.25 m depth
  -- incorporating damping of Q10 due to temperature acclimation (Lloyd & Taylor 1994)

  gtemp : real,
  --- soil organic matter (SOM) pool with c. 1000 yr turnover (kgC/m2)
  cpool_slow : real,
  --- soil organic matter (SOM) pool with c. 33 yr turnover (kgC/m2)
  cpool_fast : real,

  -- Running sums (converted to long term means) maintained by SOM dynamics module

  --- mean annual litter decomposition (kgC/m2/yr)
  decomp_litter_mean : real,
  --- mean value of decay constant for fast SOM fraction
  k_soilfast_mean : real,
  --- mean value of decay constant for slow SOM fraction
  k_soilslow_mean : real,


  -- Parameters used by function soiltemp and updated monthly

  alag : real,
  exp_alag : real,

  --- water content of soil layers [0=upper layer] as fraction of available water holding capacity
  mwcont : [12][NSOILLAYER]real,
  --- daily water content in lower soil layer for each day of year
  dwcontlower : [Date_MAX_YEAR_LENGTH]real,
  --- mean water content in lower soil layer for last month
  --- (valid only on last day of month following call to daily_accounting_patch)--
  mwcontlower : real,

  --- rainfall and snowmelt today (mm)
  rain_melt : real,
  --- upper limit for percolation (mm)
  max_rain_melt : real,
  --- whether to percolate today
  percolate : bool,

  ---------------------------------------------------------------------------------/
  -- Soil member variables

  --- sum of soil layers
  ngroundl : int,
  --- density of the snowpack, daily
  snowdens : real,


  -- Temperature variables:

  --- temperature in each layer today (after call to cnstep) [deg C]
  T : [NLAYERS]real,
  --- Recorded T each day of the year, at each level [deg C]
  T_soil : [NLAYERS]real,
  --- Record the monthly average soil temp at SOILTEMPOUT layers [deg C]
  T_soil_monthly : [12][SOILTEMPOUT]real,
  --- soil temperature from previous time step
  T_old : [NLAYERS]real,
  --- soil temperature in each layer YESTERDAY
  T_soil_yesterday : [NLAYERS]real,
  --- soil temperature at 25 cm depth, as calculated using previous versions of the model [deg C]
  temp_analyticsoln : real,


  -- Padding - The values initialised in first call to calctemp method

  --- temperature in padding layers [deg C]
  pad_temp : [PAD_LAYERS]real,
  --- thickness of padding layers [mm]
  pad_dz : [PAD_LAYERS]real,


  -- Ice and water variables for the soil layers, where Frac stands for Fraction

  --- (Frac)tion of ice in each layer: amount of ice / total volume of soil layer. Not associated with AWC.
  Frac_ice : [NLAYERS]real,

  --- ice fraction in each layer YESTERDAY
  Frac_ice_yesterday : [NLAYERS]real,

  --- fraction of water in each layer: water / total volume of soil layer
  Frac_water : [NLAYERS]real,


  -- Layer composition and properties:
  -- Soil layer information:

  --- Thickness of soil layers [mm]
  Dz : [NLAYERS]real,
  --- porosity of each soil layer
  por : [NLAYERS]real,
  --- organic soil fraction
  Frac_org : [NLAYERS]real,
  --- peat fraction
  Frac_peat : [NLAYERS]real,
  --- mineral soil fraction
  Frac_min : [NLAYERS]real,
  --- Fraction of water held below permanent wilting point when no freexing
  Fpwp_ref : [NLAYERS]real,
  --- fraction of water held below permanent wilting point.
  Frac_water_belowpwp : [NLAYERS]real,
  --- diffusivity of the soil layers [mm2 / day]
  Di : [NLAYERS]real,
  --- heat capacities of the soil layers [J mm-3 K-1]
  Ci : [NLAYERS]real,
  --- thermal conductivities of the soil layers [J day-1 mm-1 K-1]
  Ki : [NLAYERS]real,

  snow_active : bool, -- active snow layer(s)?
  snow_active_layers : int, -- <= NLAYERS_SNOW
  --- liquid water in the snow pack [kg m-2] == [mm]
  snow_water : [NLAYERS_SNOW]real,
  --- ice in the snowpack [kg m-2]
  snow_ice : [NLAYERS_SNOW]real,

  -- Log (l) of thermal conductivities. Used in the calculation of thermal conductivity and diffusivity.
  lKorg : real,
  lKpeat : real,
  lKmin : real,
  lKwater : real,
  lKice : real,
  lKair : real,

  --- records the first time T is calculated
  firstTempCalc : bool,


  -- Daily storage:

  --- Depth in mm of the maximum depth of thaw this year
  maxthawdepththisyear : real,
  --- daily thaw depth. The depth to the first soil layer with a temperature greater than 0 degrees C [mm]
  thaw : real,
  --- monthly thawing depth full, where ALL the ice has melted [mm]
  mthaw : [12]real,
  --- daily thawing depth full, where ALL the ice has melted [mm]
  dthaw : [Date_MAX_YEAR_LENGTH]real,

  --- depth of the acrotelm [mm]
  acro_depth : real,
  --- depth of the catotelm [mm]
  cato_depth : real,
  --- acrotelm CO2 level [mimil L-1]
  acro_co2 : real,


  -- Snow:

  --- days of continuous snow cover
  snow_days : int,
  --- previous days of continuous snow cover
  snow_days_prev : int,
  --- daily snow depth [mm]
  dsnowdepth : real,
  --- Monthly snow depth (average) [mm]
  msnowdepth : [12]real,
  --- Previous December's snowdepth [mm] - used in establishment - from Wolf et al. (2008)
  dec_snowdepth : real,


  -- Daily photosynthetic limits:

  --- Daily limit on moss photosynthetic activity due to dessication [0,1], but really [0.3,1]
  dmoss_wtp_limit : real,
  --- Daily limit on graminoid photosynthetic activity as WTP drops.
  dgraminoid_wtp_limit : real,
    --- acrotelm porosity MINUS Fgas (so, 0.98-0.08)
  acro_por : real,
  --- acrotelm porosity MINUS Fgas (so, 0.92-0.08)
  cato_por : real,


  -- Peatland hydrology variables:

  --- daily water table position [mm]
  wtp : [Date_MAX_YEAR_LENGTH]real,
  --- monthly average water table position [mm]
  mwtp : [12]real,
  --- annual average water table position [mm]
  awtp : real,
  --- water table depth from yesterday and updated today = -wtp [mm]
  wtd : real,
  --- Water in acrotelm plus standing water (up to a max of 100mm) [mm]
  Wtot : real,
  --- Standing water (up to a max of 100mm) [mm] - set to 0 in LPJG
  stand_water : real,
  --- Volumetric water content in the NSUBLAYERS_ACRO of the acrotelm
  sub_water : [NSUBLAYERS_ACRO]real,
  -- available water holding capacity of soil layers [0=upper layer] [mm], taking into
  -- account the unavailability of frozen water. Default value: soiltype.awc[]
  whc : [NSOILLAYER]real,
   -- available water holding capacity of evap soil layers [mm], taking into
  -- account the unavailability of frozen water. Default value: (2/5) * soiltype.awc[]
  whc_evap : real,
  --- Max water (mm) that can be held in each layer
  aw_max : [NSOILLAYER]real,
  -- Volumetric liquid water content. A fraction. Considers the entire (awc + Fpwp)
  -- volumetric water content MINUS the ice fraction. Updated daily.
  alwhc : [NSOILLAYER]real,
  -- Initial volumetric liquid water content. A fraction. Considers the entire
  -- (awc + Fpwp) volumetric water cont,ent MINUS the ice fraction.
  alwhc_init : [NSOILLAYER]real,


  -- Indices

  --- index for first active layer of soil
  IDX : int,
  --- index for snow layers
  SIDX : int,
  --- index for snow layer from previous day
  SIDX_old : int,
  --- index for mixed layer
  MIDX : int,

  -- Hydrology variables:

  --- records the first time hydrology routine is called
  firstHydrologyCalc : bool,

  -- These are the number of sublayers in the standard 0.5/1.0m
  -- hydrology laters. Set ONCE in hydrology routine
  nsublayer1 : int,
  nsublayer2 : int,
  num_evaplayers : int,

  --- root fractions per layer
  rootfrac : [NLAYERS]real,
  --- air fraction in each layer
  Frac_air : [NLAYERS]real,

  --- daily carbon flux to atmosphere from soil respiration
  --- Temporary storage of heterotrophic respiration on PEATLAND until it is reduced by allocation of a certain fraction to CH4 production.
  dcflux_soil : real,

  --- CH4 and CO2 stores in the soil layers - updated daily
  ch4_store : real,
  co2_store : real,

  --- dissolved CO2 concentration in each layer [g CO2-C layer-1 d-1]
  CO2_soil : [NLAYERS]real,
  --- dissolved CO2 concentration in each layer yesterday [g CO2-C layer-1 d-1]
  CO2_soil_yesterday : [NLAYERS]real,
  --- daily CO2 production in each layer [g CO2-C layer-1 d-1]
  CO2_soil_prod : [NLAYERS]real,
  --- dissolved CH4 concentration in each layer [g CH4-C layer-1 d-1]
  CH4 : [NLAYERS]real,
  --- dissolved CH4 concentration in each layer yesterday [g CH4-C layer-1 d-1]
  CH4_yesterday : [NLAYERS]real,
  --- daily CH4 production in each layer [g CH4-C layer-1 d-1]
  CH4_prod : [NLAYERS]real,
  --- daily CH4 oxidation in each layer [g CH4-C layer-1 d-1]
  CH4_oxid : [NLAYERS]real,
  --- CH4 which bubbles out [g CH4-C layer-1]
  CH4_ebull_ind : [NLAYERS]real,
  --- Volume of CH4 which bubbles out [m3]
  CH4_ebull_vol : [NLAYERS]real,
  --- gaseous CH4 concentration in each layer [g CH4-C layer-1 d-1]
  CH4_gas : [NLAYERS]real,
  --- dissolved CH4 concentration in each layer [g CH4-C layer-1 d-1]
  CH4_diss : [NLAYERS]real,
  --- gaseous CH4 concentration in each lapftidyer [g CH4-C layer-1] yesterday
  CH4_gas_yesterday : [NLAYERS]real,
  --- dissolved CH4 concentration in each layer [g CH4-C layer-1] yesterday
  CH4_diss_yesterday : [NLAYERS]real,
  --- gaseous CH4 volume in each layer [m3]
  CH4_gas_vol : [NLAYERS]real,
  -- volumetric CH4 content [unitless]
  CH4_vgc : [NLAYERS]real,
  --- dissolved O2 concentration in each layer [mol O2 layer-1 d-1]
  O2 : [NLAYERS]real,

  --- layer water volume [m3]
  volume_liquid_water : [NLAYERS]real,
  --- layer water + ice volume [m3]
  total_volume_water : [NLAYERS]real,
  --- tiller area
  tiller_area : [NLAYERS]real,


  -- Gas diffusion variables:

  --- gas transport velocity of O2 [m d-1]
  k_O2 : real,
  --- gas transport velocity of CO2 [m d-1]
  k_CO2 : real,
  --- gas transport velocity of CH4 [m d-1]
  k_CH4 : real,
  --- Equilibrium concentration of O2 [mol L-1]
  Ceq_O2 : real,
  --- Equilibrium concentration of CO2 [mol L-1]
  Ceq_CO2 : real,
  --- Equilibrium concentration of CH4 [mol L-1]
  Ceq_CH4 : real,


---------------------------------------------------------------------------------/
-- CENTURY SOM pools and other variables

  sompool : [NSOMPOOL]Sompool,

  --- daily percolation (mm)
  dperc : real,
  --- fraction of decayed organic nitrogen leached each day,
  orgleachfrac : real,
  --- soil NH4 mass in pool (kgN/m2)
  NH4_mass : real,
  --- soil NO3 mass in pool (kgN/m2)
  NO3_mass : real,
  --- soil NH4 mass input (kgN/m2)
  NH4_input : real,
  --- soil NO3 mass input (kgN/m2)
  NO3_input : real,
  --- annual sum of nitrogen mineralisation
  anmin : real,
  --- annual sum of nitrogen immobilisation
  animmob : real,
  --- annual leaching from available nitrogen pool
  aminleach : real,
  --- annual leaching of organics from active nitrogen pool
  aorgNleach : real,
  --- total annual nitrogen fixation
  anfix : real,
  --- calculated annual mean nitrogen fixation
  anfix_calc : real,
  --- annual leaching of organics nitrogen from carbon pool
  aorgCleach : real,

  -- Variables for fast spinup of SOM pools

  --- monthly fraction of available mineral nitrogen taken up
  fnuptake_mean : [12]real,
  --- monthly fraction of organic carbon/nitrogen leached
  morgleach_mean : [12]real,
  --- monthly fraction of available mineral nitrogen leached
  mminleach_mean : [12]real,
  --- annual nitrogen fixation
  anfix_mean : real,

  -- Solving Century SOM pools

  --- years at which to begin documenting for calculation of Century equilibrium
  solvesomcent_beginyr : int,
  --- years at which to end documentation and start calculation of Century equilibrium
  solvesomcent_endyr : int,

  --- Cumulative litter pools for one year.
  litterSolveSOM : LitterSolveSOM,

  --TODO: this one is a vector in c++, but it's seemingly never initialized or used
  --FIXME WAS UNUSED IN C++ CODE
  --solvesom : []LitterSolveSOM,

  --- stored NH4 deposition in snowpack
  snowpack_NH4_mass : real,
  --- stored NO3 deposition in snowpack
  snowpack_NO3_mass : real,

  --- pools of soil N species in transformation (nitrification & denitrifiacation)

  --- soil NH4 mass in pool (kgN/m2)
  -- NH4_mass : real,  -- total, definde above
  NH4_mass_w : real,  -- wet proportion
  NH4_mass_d : real,  -- dry...

  --- soil NO3 mass in pool (kgN/m2)
  -- NO3_mass : real, -- total, definde above
  NO3_mass_w : real,
  NO3_mass_d : real,
  --- soil NO2 mass in pool (kgN/m2)
  NO2_mass : real,
  NO2_mass_w : real,
  NO2_mass_d : real,
  --- soil NO mass in pool (kgN/m2)
  NO_mass : real,
  NO_mass_w : real,
  NO_mass_d : real,
  --- soil NO mass in pool (kgN/m2)
  N2O_mass : real,
  N2O_mass_w : real,
  N2O_mass_d : real,
  --- soil N2 mass in pool (kgN/m2)
  N2_mass : real,

  -- soil pH
  pH : real,  --TODO: pH - not used yet. Daily mean precip, based on annual average

  -- soil labile carbon availability daily (kgC/m2/day)
  labile_carbon : real,
  labile_carbon_w : real,
  labile_carbon_d : real
}


let Soiltype(soiltype_id : int) : Soiltype = {
  soiltype_id = soiltype_id,
  solvesom_end = SOLVESOM_END,
  solvesom_begin = SOLVESOM_BEGIN,
  organic_frac = 0.02,
  pH = -1.0,

  -- Assume no mineral content on peatlands
  sand_frac_peat = 0.0,
  clay_frac_peat = 0.0,
  silt_frac_peat = 0.0,

  runon = 0.0,
  soiltempdepths = replicate 10 0.0,

  -- NOTE: These were UNINITIALIZED in the c++ code:
  awc = replicate NSOILLAYER 0.0, --replicate NSOILLAYER nan, --- FIXME: awc is used EVERYWHERE but never initialized in the c++ code, how?
  awc_peat = replicate NSOILLAYER nan,
  clay_frac = nan,
  gawc = replicate 2 nan,
  gawc_peat = replicate 2 nan,
  gwp = replicate 2 nan,
  gwp_peat = replicate 2 nan,
  gwsats = replicate 2 nan,
  gwsats_peat = replicate 2 nan,
  min_frac_gridcell = replicate NSOILLAYER nan,
  mineral_frac = nan,
  org_frac_gridcell = replicate NSOILLAYER nan,
  perc_base = nan,
  perc_base_evap = nan,
  perc_exp = nan,
  porosity = nan,
  porosity_gridcell = replicate NSOILLAYER nan,
  sand_frac = nan,
  silt_frac = nan,
  soilcode = intnan,
  thermdiff_0 = nan,
  thermdiff_100 = nan,
  thermdiff_15 = nan,
  water_below_wp = nan,
  wp = replicate NSOILLAYER nan,
  wp_peat = replicate NSOILLAYER nan,
  wsats = replicate NSOILLAYER nan,
  wsats_peat = replicate NSOILLAYER nan,
  wtot = nan,
  wtot_peat = nan
}


--- Override the default SOM years with 70-80% of the spin-up period length
let updateSolveSOMvalues(this : Soiltype, nyrspinup : int) : Soiltype =
  let this = this with solvesom_end = intFromReal (0.8 * (realFromInt nyrspinup))
  let this = this with solvesom_begin = intFromReal (0.7 * (realFromInt nyrspinup))
  in this

--- Constructor
let Sompool() : Sompool = {
  -- Initialise pool
  cmass = 0.0,
  nmass = 0.0,
  ligcfrac = 0.0,
  delta_cmass = 0.0,
  delta_nmass = 0.0,
  fracremain = 0.0,
  litterme = 0.0,
  fireresist = 0.0,
  mfracremain_mean = replicate 12 0.0,

  -- FIXME: these were not initialized in the c++ code:
  cdec = nan,
  ndec = nan,
  ntoc = nan
}


let LitterSolveSOM() : LitterSolveSOM =
  { clitter = replicate NSOMPOOL 0.0
  , nlitter = replicate NSOMPOOL 0.0
  }

--- Add litter
let add_litter({clitter, nlitter} : *LitterSolveSOM, cvalue : real, nvalue : real, pool : int) : LitterSolveSOM =
  { clitter = clitter with [pool] = clitter[pool] + cvalue
  , nlitter = nlitter with [pool] = nlitter[pool] + nvalue
  }

let Soil(soiltype : Soiltype) : Soil = {
  soiltype_id = soiltype.soiltype_id,

  -- Initialises certain member variables
  alag = 0.0,
  exp_alag = 1.0,
  cpool_slow = 0.0,
  cpool_fast = 0.0,
  decomp_litter_mean = 0.0,
  k_soilfast_mean = 0.0,
  k_soilslow_mean = 0.0,
  wcont_evap = 0.0,
  snowpack = 0.0,
  orgleachfrac = 0.0,

  -- Extra initialisation
  aorgCleach = 0.0,
  aorgNleach = 0.0,
  anfix = 0.0,
  aminleach = 0.0,

  lKorg = 0.0,
  lKpeat = 0.0,
  lKmin = 0.0,
  lKwater = 0.0,
  lKice = 0.0,
  lKair = 0.0,

  mwcontupper = 0.0,
  mwcontlower = 0.0,

  --for (int mth=0, mth<12, mth++) {
  --  mwcont[mth][0] = 0.0,
  --  mwcont[mth][1] = 0.0,
  --  fnuptake_mean[mth] = 0.0,
  --  morgleach_mean[mth] = 0.0,
  --  mminleach_mean[mth] = 0.0,
  --}
  mwcont = replicate 12 (replicate NSOILLAYER 0.0),
  morgleach_mean = replicate 12 0.0,
  mminleach_mean = replicate 12 0.0,
  fnuptake_mean = replicate 12 0.0,

  dwcontupper = replicate Date_MAX_YEAR_LENGTH 0.0,
  dwcontlower = replicate Date_MAX_YEAR_LENGTH 0.0,
  -- for fire
  dthaw = replicate Date_MAX_YEAR_LENGTH 0.0,


  ----------------------------------------------------/
  -- Initialise CENTURY pools

  -- Set initial CENTURY pool N:C ratios
  -- Parton et al 1993, Fig 4

  -- inplace updates on records inside arrays is impossible?
  -- easier to do a scatter
  -- create the 'new' records and scatter them to their places
  sompool = (let somepool = Sompool()
             let dsts = replicate NSOMPOOL somepool
             let idx_vals = [(SOILMICRO, 1.0 / 15.0)
                            ,(SURFHUMUS, 1.0 / 15.0)
                            ,(SLOWSOM  , 1.0 / 20.0)
                            ,(SURFMICRO, 1.0 / 20.0)
                            -- passive has a fixed value
                            ,(PASSIVESOM,1.0 /  9.0)
                            ]
             let updated_Sompools = map (\(_, b) ->
                Sompool() with ntoc = b
               ) idx_vals
             in scatter dsts (map (\(a,_) -> a) idx_vals) updated_Sompools
             ),


  NO2_mass = 0.0,
  NO2_mass_w = 0.0,
  NO2_mass_d = 0.0,
  NO_mass = 0.0,
  NO_mass_w = 0.0,
  NO_mass_d = 0.0,
  N2O_mass = 0.0,
  N2O_mass_w = 0.0,
  N2O_mass_d = 0.0,
  N2_mass = 0.0,
  NH4_mass = 0.0,
  NO3_mass = 0.0,
  NH4_input = 0.0,
  NO3_input = 0.0,
  anmin = 0.0,
  animmob = 0.0,
  --aminleach = 0.0, -- these were initialied twice in c++
  --aorgNleach = 0.0,
  --aorgCleach = 0.0,
  --anfix = 0.0,
  anfix_calc = 0.0,
  anfix_mean = 0.0,
  snowpack_NH4_mass = 0.0,
  snowpack_NO3_mass = 0.0,
  labile_carbon = 0.0,
  labile_carbon_w = 0.0,
  labile_carbon_d = 0.0,

  pH = 6.5,

  dperc = 0.0,

  solvesomcent_beginyr = intFromReal (SOLVESOMCENT_SPINBEGIN * (realFromInt <| nyear_spinup - freenyears) + (realFromInt freenyears)),
  solvesomcent_endyr   = intFromReal (SOLVESOMCENT_SPINEND   * (realFromInt <|nyear_spinup - freenyears) + (realFromInt freenyears)),

  temp25 = 10,

  ----------------------------------------------------/
  -- Arctic and wetland initialisation

  -- Daily heterotrophic respiration. Only ever nonzero on peatland stands.
  dcflux_soil = 0.0,

  -- Indices
  IDX = 0,
  SIDX = 0,
  MIDX = 0,
  SIDX_old = 0,

  firstTempCalc = true,
  firstHydrologyCalc = true,
  --SIDX_old = 0, -- these were initialized multiple times
  --IDX = 0,
  --nsublayer1 = 0,
  --nsublayer2 = 0,
  --num_evaplayers = 0,

  -- initialise snow variables
  snow_active = false,
  snow_active_layers = 0,
  snow_days = 0,
  snow_days_prev = 365,
  dec_snowdepth = 0.0,

  -- Peatland hydrology variables
  awtp = 0.0,
--  acro_depth = 0.0,
--  cato_depth = 0.0,

  acro_co2 = 934.0, -- [mimol L-1]

  wtd = 0.0, -- [-100, +300] mm

  acro_por = acrotelm_por - Fgas,
  cato_por = catotelm_por - Fgas,

  Wtot = 0.0,
  stand_water = 0.0,

  dmoss_wtp_limit = 1.0,
  dgraminoid_wtp_limit = 1.0,

  k_O2 = 0.0,
  k_CO2 = 0.0,
  k_CH4 = 0.0,
  Ceq_O2 = 0.0,
  Ceq_CO2 = 0.0,
  Ceq_CH4 = 0.0,

  ch4_store = 0.0,
  co2_store = 0.0,

  maxthawdepththisyear = 0.0,
  thaw = 0.0,
  snowdens = snowdens_start, -- kg/m3
  dsnowdepth = 0.0,

  mthaw = replicate 12 0.0,
  msnowdepth = replicate 12 0.0,
  mwtp = replicate 12 0.0,
  T_soil_monthly = replicate 12 (replicate SOILTEMPOUT (-999.0)),

  sub_water = replicate NSUBLAYERS_ACRO 0.0,

  wtp = replicate Date_MAX_YEAR_LENGTH 0.0,

  -- FUNFACT: these were initialized to 0, 365 times in a row, in the c++ code
  Frac_ice = replicate NLAYERS 0.0,
  T_soil = replicate NLAYERS 0.0,

  -- FIX?: The original code was super wierd and inefficient here.
  Frac_water = replicate NLAYERS 0.0,
  T = replicate NLAYERS 0.0,
  T_old = replicate NLAYERS 0.0,
  Frac_org = replicate NLAYERS 0.0,
  Frac_peat = replicate NLAYERS 0.0,
  Frac_min = replicate NLAYERS 0.0,
  Fpwp_ref = replicate NLAYERS 0.0,
  Frac_water_belowpwp = replicate NLAYERS 0.0,
  Frac_ice_yesterday = replicate NLAYERS 0.0,
  Dz = replicate NLAYERS 0.0,
  por = replicate NLAYERS 0.0,
  rootfrac = replicate NLAYERS 0.0,
  CH4_yesterday = replicate NLAYERS 0.0,
  CO2_soil_yesterday = replicate NLAYERS 0.0,
  volume_liquid_water = replicate NLAYERS 0.0,
  total_volume_water = replicate NLAYERS 0.0,
  tiller_area = replicate NLAYERS 0.0001,  -- The max value
  CH4_ebull_ind = replicate NLAYERS 0.0,
  CH4_ebull_vol = replicate NLAYERS 0.0,
  CH4_gas_yesterday = replicate NLAYERS 0.0,
  CH4_diss_yesterday = replicate NLAYERS 0.0,
  CH4_gas_vol = replicate NLAYERS 0.0,
  CH4 = replicate NLAYERS 0.0,
  CO2_soil = replicate NLAYERS 0.0,
  O2 = replicate NLAYERS 0.0,
  CO2_soil_prod = replicate NLAYERS 0.0,
  CH4_prod = replicate NLAYERS 0.0,
  CH4_gas = replicate NLAYERS 0.0,
  CH4_diss = replicate NLAYERS 0.0,
  Frac_air = replicate NLAYERS 0.0,


  wcont = replicate NSOILLAYER 0.0,
  alwhc = replicate NSOILLAYER 0.0,
  alwhc_init = replicate NSOILLAYER 0.0,

  awcont_upper = 0.0,

  nsublayer1 = intFromReal(SOILDEPTH_UPPER / Dz_soil),    -- 5, typically
  nsublayer2 = intFromReal(SOILDEPTH_LOWER / Dz_soil),    -- 10, typically
  num_evaplayers = intFromReal(SOILDEPTH_EVAP / Dz_soil),  -- 2, typically

  whc_evap = soiltype.awc[0] + soiltype.awc[1],    -- mm

  -- Depths of the acrotelm & catotelm
  acro_depth = realFromInt <| NACROTELM * (intFromReal Dz_acro),
  cato_depth = SOILDEPTH_UPPER + SOILDEPTH_LOWER - (realFromInt <| NACROTELM * (intFromReal Dz_acro)),

  --- the following were not initialized in the original code:
  CH4_oxid = replicate NLAYERS nan,
  CH4_vgc = replicate NLAYERS nan,
  Ci = replicate NLAYERS nan,
  Di = replicate NLAYERS nan,
  Ki = replicate NLAYERS nan,
  NH4_mass_d = nan,
  NH4_mass_w = nan,
  NO3_mass_d = nan,
  NO3_mass_w = nan,
  T_soil_yesterday = replicate NLAYERS nan,
  aw_max = replicate NSOILLAYER nan,
  dtemp = replicate 31 nan,
  gtemp = nan,
  litterSolveSOM = LitterSolveSOM(),
  max_rain_melt = nan,
  mtemp = nan,
  ngroundl = intnan,
  pad_dz = replicate PAD_LAYERS nan,
  pad_temp = replicate PAD_LAYERS nan,
  percolate = false,
  rain_melt = nan,
  runoff = nan,
  snow_ice = replicate NLAYERS_SNOW nan,
  snow_water = replicate NLAYERS_SNOW nan,
  temp_analyticsoln = nan,
  whc = replicate NSOILLAYER nan
}



-- TODO: There are many SOIL member functions in soil.cpp

let get_soil_temp_25(this: Soil) : real =
	-- Return soil temperature at 25cm depth
	if (iftwolayersoil) then
		this.temp_analyticsoln
	else
		this.temp25




-- INLINE FUNCTIONS

--let tridiag(int n, long double a[], long double b[], long double c[], long double r[], long double u[]) {
--
--	-- Tridiagonal system solver from Numerical Recipes.
--
--	double gam[active_layersmax]
--	double bet
--
--	bet : real =  b[0]
--
--	u[0] : real =  r[0] / bet
--
--	for (int j : real =  1 j<n j++) {
--		gam[j] : real =  c[j - 1] / bet
--		bet : real =  b[j] - a[j] * gam[j]
--
--		u[j] : real =  (r[j] - a[j] * u[j - 1]) / bet
--	}
--
--	for (int j : real =  (n - 2) j >= 0 j--) {
--		u[j] -= gam[j + 1] * u[j + 1]
--	}
--}


------------------------------------------------------------------------------/

-- References

-- Aerts, R., Verhoeven, J.T.A., & Whigham, D.F. (1999) Plant-mediated controls on nutrient cycling
--   in temperate fenns and bogs, Ecology, 80 (7), 2170-2181.
-- Best et al. (2011) Geosci. Model Dev., 4, 677-699. www.geosci-model-dev.net/4/677/2011/
--   doi:10.5194/gmd-4-677-2011
-- Bonan, G.B. 2002 Ecological Climatology, Cambridge University Press
-- Clein, J. S., and J. P. Schimel, Microbial activity of tundra and taiga soils at sub-zero
--   temperatures, Soil Biol. Biochem., 27(9), 1231-1234, 1995.
-- Cronk, J. K. and Fennessy, M. S.: Wetland Plants: Biology and Ecology, CRC Press LLC, 2001.
-- Ekici, A., Chadburn, S., Chaudhary, N., Hajdu, L. H., Marmy, A., Peng, S., Boike, J., Burke, E.,
--   Friend, A. D., Hauck, C., Krinner, G., Langer, M., Miller, P. A., and Beer, C.: Site-level model
--   intercomparison of high latitude and high altitude soil thermal dynamics in tundra and barren
--   landscapes, The Cryosphere, 9, 1343-1361, https:--doi.org/10.5194/tc-9-1343-2015, 2015.
-- Frolking, S., Roulet, N. T., Moore, T. R., Richard, P. J. H., Lavoie, M., and Muller, S. D.: Modeling
--   northern peatland decomposition and peat accumulation, Ecosystems, 4, 479-498, 2001.
-- Frolking, S., Roulet, N. T., Tuittila, E., Bubier, J. L., Quillet, A., Talbot, J., and Richard,
--   P. J. H.: A new model of Holocene peatland net primary production, decomposition, water balance,
--   and peat accumulation, Earth Syst. Dynam., 1, 1-21, https:--doi.org/10.5194/esd-1-1-2010, 2010.
-- Fukusako, S. (1990) Thermophysical properties of ice, snow, and sea ice. Int. J. Thermophys., 11(2):
--	 353-372. doi:10.1007/BF01133567
-- Granberg, et al. 1999 A simple model for simulation of water content, soil frost,
--   and soil temperatures in boreal mixed mires. Water Resour. Res., 35(12), 3771-3782.
-- Hillel D., (1982) Introduction to soil physics. Academic Press, San Diego, CA, USA.
-- Ise, T., Dunn, A.L, Wofsy, S.C. & Moorcroft, P.R.: High sensitivity of peat decomposition to climate
--   change double Soil::nmass_avail(int pref) {through water-table feedback. Nature Geoscience volume 1, pages 763-766 (2008)
-- Jaehne, B., Heinz, G., and Dietrich, W.: Measurement of the diffusion coefficients of sparingly
--   soluble gases in water, J. Geophys. Res., 92, 10767-10776, 1987.
-- Koven, C.D, Ringeval, B., Friedlingstein, P., Ciais, P., Cadule, P., Khvorostyanov, D.,
--   Krinner, G., and Tarnocai C. (2011) Permafrost carbon-climate feedbacks accelerate global warming,
--   PNAS, vol. 108 no. 36, 14769-14774.
-- Lawrence, D. M., and A. G. Slater, 2008: Incorporating organic soil into a global climate model.
--   Climate Dynamics, 30, 145-160, doi:10.1007/s00382-007-0278-1.
-- Ling, F., and Zhang, T. (2006) Sensitivity of ground thermal regime and surface energy fluxes to
--   tundra snow density in northern Alaska. Cold Regions Science and Technology 44 (2006) 121-130
-- McGuire, A. D., Christensen, T. R., Hayes, D., Heroult, A., Euskirchen, E., Kimball, J. S., Koven, C.,
--   Lafleur, P., Miller, P. A., Oechel, W., Peylin, P., Williams, M., and Yi, Y.: An assessment of
--   the carbon balance of Arctic tundra: comparisons among observations, process models, and atmospheric
--   inversions, Biogeosciences, 9, 3185-3204, https:--doi.org/10.5194/bg-9-3185-2012, 2012.
-- Potter, C. S., Davidson, E. A., and Verchot, L. V.: Estimation of global biogeochemical controls and
--   seasonality in soil methane consumption, Chemosphere, 32, 2219-2246, 1996.
-- Riera, J. L., Schindler, J. E., and Kratz, T. K.: Seasonal dynamics of carbon dioxide and methane in
--   two clear-water lakes and two bog lakes in northern Wisconsin, USA,
--   Can. J. Fish. Aquat. Sci., 56, 265-274, 1999.
-- Sander, R.: Compilation of Henry's Law Constants for inorganic and organic species of potential
--   importance in environmental chemistry, Tech. Rep. Version 3, MPI Mainz, Air Chemistry Department,
--   Max-Planck Institute of Chemistry, 1999.
-- Schimel, J. P.: Plant transport and methane production as controls on methane flux from arctic wet
--   meadow tundra, Biogeochem., 28, 183-200, 1995.
-- Smolders, A. J. P., H. B. M. Tomassen, H. W. Pijnappel, L. P. M. Lamers, and J. G. M. Roelofs (2001),
--   Substrate-derived CO2 is important in the development of Sphagnum spp., New Phytol., 152(2), 325- 332.
-- Spahni, R., Wania, R., Neef, L., van Weele, M., Pison, I., Bousquet, P., Frankenberg, C., Foster,
--   P. N., Joos, F., Prentice, I. C., and van Velthoven, P.: Constraining global methane emissions and
--   uptake by ecosystems, Biogeosciences, 8, 1643-1665, https:--doi.org/10.5194/bg-8-1643-2011, 2011.
-- Sturm, M., Holmgren, J., Konig, M., and Morris, K. (1997) The thermal conductivity of seasonal snow.
--   Journal of Glaciology 43 (143), 26-41.
-- Swenson, S.C., Lawrence, D.M., and Lee, H. 2012. Improved Simulation of the Terrestrial Hydrological
--   Cycle in Permafrost Regions by the Community Land Model. JAMES, 4, M08002. DOI:10.1029/2012MS000165.
-- Tang, J., Miller, P. A., Persson, A., Olefeldt, D., Pilesjo, P., Heliasz, M., Jackowicz-Korczynski,
--   M., Yang, Z., Smith, B., Callaghan, T. V., and Christensen, T. R.: Carbon budget estimation of a
--   subarctic catchment using a dynamic ecosystem model at high spatial resolution,
--   Biogeosciences, 12, 2791-2808, doi:10.5194/bg-12-2791-2015, 2015.
-- Wania, R., Ross, I., & Prentice, I.C. (2009a) Integrating peatlands and permafrost
--   into a dynamic global vegetation model: I. Evaluation and sensitivity of physical
--   land surface processes. Global Biogeochemical Cycles, 23, GB3014, doi:10.1029/2008GB003412
-- Wania, R., Ross, I., & Prentice, I.C. (2009b) Integrating peatlands and permafrost
--   into a dynamic global vegetation model: II. Evaluation and sensitivity of vegetation
--   and carbon cycle processes. Global Biogeochemical Cycles, 23, GB015, doi:10.1029/2008GB003413
-- Wania, R., Ross, I., & Prentice, I.C. (2010) Implementation and evaluation of a new methane
--   model within a dynamic global vegetation model: LPJ-WHyMe v1.3.1, Geosci. Model Dev., 3, 565-584.
-- Wisser, D., Marchenko, S., Talbot, J., Treat, C., and Frolking, S. (2011) Soil temperature response
--   to 21st century global warming: the role of and some implications for peat carbon in thawing
--   permafrost soils in North America, Earth Syst. Dynam., 2, 121-138
-- Wolf, A., Callaghan T.V., & Larson K. (2008) Future changes in vegetation and ecosystem
--   function of the Barents Region. Climatic Change, 87:51-73 DOI 10.1007/s10584-007-9342-4
-- Yurova, A., Wolf, A., Sagerfors, J., & Nilsson, M. (2007) Variations in net ecosystem
--   exchange of carbon dioxide in a boreal mire: Modeling mechanisms linked to water table
--   position, Journal of Geophysical Research, 112, art. no. G02025, doi:10.1029/2006JG000342.
-- Zhang, W., Miller, P.A., Smith, B., Wania, R., Koenigk, T. & Doscher, R., 2013,
--   Tundra shrubification and tree-line advance amplify arctic climate warming: results from an
--   individual-based dynamic vegetation model. Environmental Research Letters 8: 034023.

--#endif --LPJ_GUESS_SOIL_H
