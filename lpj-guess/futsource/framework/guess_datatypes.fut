open import "../framework/guessmath"
--#include "archive.h"
open import "../framework/parameters"
--#include "guesscontainer.h"
open import "../modules/soil"
open import "../futhark-extras"



----------------------------------------------------------
-- GLOBAL CONSTANTS

-- number  of soil layers modelled
let NSOILLAYER_UPPER : int = 5
let NSOILLAYER_LOWER : int = NSOILLAYER - NSOILLAYER_UPPER

-- bvoc: number of monoterpene species used
let NMTCOMPOUNDS : int = NMTCOMPOUNDTYPES


-- Number of years to average growth efficiency over in function mortality
let NYEARGREFF : int = 5

-- Coldest day in N hemisphere (January 15)
-- Used to decide when to start counting GDD's and leaf-on days
--  for summergreen phenology.
--/
let COLDEST_DAY_NHEMISPHERE : int = 14

-- Coldest day in S hemisphere (July 15)
-- Used to decide when to start counting GDD's and leaf-on days
--  for summergreen phenology.
--/
let COLDEST_DAY_SHEMISPHERE : int = 195

-- Warmest day in N hemisphere (same as COLDEST_DAY_SHEMISPHERE)
let WARMEST_DAY_NHEMISPHERE : int = COLDEST_DAY_SHEMISPHERE

-- Warmest day in S hemisphere (same as COLDEST_DAY_NHEMISPHERE)
let WARMEST_DAY_SHEMISPHERE : int = COLDEST_DAY_NHEMISPHERE

-- number of years to average aaet over in function soilnadd
let NYEARAAET : int = 5

-- number of years to average max snow depth over in function soilnadd
let NYEARMAXSNOW : int = 20

-- Priestley-Taylor coefficient (conversion factor from equilibrium evapotranspiration to PET)
let PRIESTLEY_TAYLOR : real = 1.32


----------------------------------------------------------
-- GLOBAL ENUMERATED TYPE DEFINITIONS

-- Life form class for PFTs (trees, grasses)
type lifeformtype = enum_type
let NOLIFEFORM : lifeformtype = 0
let TREE : lifeformtype = 1
let GRASS : lifeformtype = 2
let MOSS : lifeformtype = 3

-- Phenology class for PFTs
type phenologytype = enum_type
let NOPHENOLOGY : phenologytype = 0
let EVERGREEN : phenologytype = 1
let RAINGREEN : phenologytype = 2
let SUMMERGREEN : phenologytype = 3
let CROPGREEN : phenologytype = 4
let ANY : phenologytype = 5
-- Biochemical pathway for photosynthesis (C3 or C4)
type pathwaytype = enum_type
let NOPATHWAY : pathwaytype = 0
let C3 : pathwaytype = 1
let C4 : pathwaytype = 2
-- Leaf physiognomy types for PFTs
type leafphysiognomytype = enum_type
let NOLEAFTYPE : leafphysiognomytype = 0
let NEEDLELEAF : leafphysiognomytype = 1
let BROADLEAF : leafphysiognomytype = 2
-- The level of verbosity of LPJ-GUESS. Decides the amount of information that is written to the log-file.
type verbositylevel = enum_type
let ERROR : verbositylevel = 0
let WARNING : verbositylevel = 1
let INFO : verbositylevel = 2
let DEBUG_WARNING : verbositylevel = 3
-- Units for insolation driving data
-- Insolation can be expressed as:
--
--  - Percentage sunshine
--  - Net instantaneous downward shortwave radiation flux (W/m2)
--  - Total (i.e. with no correction for surface albedo) instantaneous downward
--    shortwave radiation flux (W/m2)
--
--  Radiation flux can be interpreted as W/m2 during daylight hours, or averaged
--  over the whole time step which it represents (24 hours in daily mode). For
--  this reason there are two enumerators for these insolation types (e.g. SWRAD
--  and SWRAD_TS).
--/


-- Irrigation type for PFTs
type hydrologytype = enum_type
let RAINFED : hydrologytype = 0
let IRRIGATED : hydrologytype = 1

-- Intercrop type for PFTs
type intercroptype = enum_type
let NOINTERCROP : intercroptype = 0
let NATURALGRASS : intercroptype = 1

--  0:SEASONALITY_NO      No seasonality
--  1:SEASONALITY_PREC      Precipitation seasonality only
--  2:SEASONALITY_PRECTEMP    Both temperature and precipitation seasonality, but "weak" temperature seasonality (coldest month > 10degC)
--  3:SEASONALITY_TEMP      Temperature seasonality only
--  4:SEASONALITY_TEMPPREC    Both temperature and precipitation seasonality, but temperature most important (coldest month < 10degC)
--  5:    Temperature seasonality, always above 10 degrees (currently not used)
type seasonality_type = enum_type
let SEASONALITY_NO : seasonality_type = 0
let SEASONALITY_PREC : seasonality_type = 1
let SEASONALITY_PRECTEMP : seasonality_type = 2
let SEASONALITY_TEMP : seasonality_type = 3
let SEASONALITY_TEMPPREC : seasonality_type = 4
--let SEASONALITY_TEMPWARM : seasonality_type = 5 -- TODO: this one was in the original comment but not source

-- Precipitation seasonality type of gridcell
-- 0:DRY            (minprec_pet20<=0.5 && maxprec_pet20<=0.5)
--  1:DRY_INTERMEDIATE      (minprec_pet20<=0.5 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
--  2:DRY_WET          (minprec_pet20<=0.5 && maxprec_pet20>1.0)
--  3:INTERMEDIATE        (minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
--  4:INTERMEDIATE_WET      (minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>1.0)
--  5:WET            (minprec_pet20>1.0 && maxprec_pet20>1.0)
--/
type prec_seasonality_type = enum_type
let DRY : prec_seasonality_type = 0
let DRY_INTERMEDIATE : prec_seasonality_type = 1
let DRY_WET : prec_seasonality_type = 2
let INTERMEDIATE : prec_seasonality_type = 3
let INTERMEDIATE_WET : prec_seasonality_type = 4
let WET : prec_seasonality_type = 5

-- Temperature seasonality type of gridcell
-- 0:COLD            (mtemp_max20<=10)
-- 1:COLD_WARM          (mtemp_min20<=10 && mtemp_max20>10 && mtemp_max20<=30)
-- 2:COLD_HOT          (mtemp_min20<=10 && mtemp_max20>30)
-- 3:WARM            (mtemp_min20>10 && mtemp_max20<=30)
-- 4:WARM_HOT          (mtemp_min20>10 && mtemp_max20>30)
-- 5:HOT            (mtemp_min20>30)
--/
type temp_seasonality_type = enum_type
let COLD : temp_seasonality_type = 0
let COLD_WARM : temp_seasonality_type = 1
let COLD_HOT : temp_seasonality_type = 2
let WARM : temp_seasonality_type = 3
let WARM_HOT : temp_seasonality_type = 4
let HOT : temp_seasonality_type = 5
-- Gas type (used in methane code)
type gastype = enum_type
let O2gas : gastype = 0
let CO2gas : gastype = 1
let CH4gas : gastype = 2


--- General purpose object for handling simulation timing.
-- In general, frameworks should use a single Date object for all simulation
--  timing.
--
--  Member variables of the class (see below) provide various kinds of calender
--  and timing information, assuming init has been called to initialise the
--  object, and next() has been called at the end of each simulation day.
type Date = {
  --- Maximum number of days in an LPJ-GUESS simulation year
  --- The standard version doesn't yet support leap years. */
  MAX_YEAR_LENGTH : int,

  --- number of days in each month (0=January - 11=December)
  ndaymonth : [12]int,

  --- julian day of year (0-364 : int, 0=Jan 1)
  day : int,

  --- day of current month (0=first day)
  dayofmonth : int,

  --- month number (0=January - 11=December)
  month : int,

  --- year since start of simulation (0=first simulation year)
  year : int,

  --- number of subdaily periods in a day (to be set in IO module)
  subdaily : int,

  --- julian day for middle day of each month
  middaymonth : [12]int,

  --- true if last year of simulation, false otherwise
  islastyear : bool,

  --- true if last month of year, false otherwise
  islastmonth : bool,

  --- true if last day of month, false otherwise
  islastday : bool,

  --- true if middle day of month, false otherwise
  ismidday : bool,

  --- The calendar year corresponding to simulation year 0
  first_calendar_year : int,

  nyear : int
}
-- type Day


type MassBalance = {
  start_year : int,
  ccont : real,
  ccont_zero : real,
  ccont_zero_scaled : real,
  cflux : real,
  cflux_zero : real,

  ncont : real,
  ncont_zero : real,
  ncont_zero_scaled : real,
  nflux : real,
  nflux_zero : real
}


-- This struct contains the result of a photosynthesis calculation.
-- see photosynthesis
type PhotosynthesisResult = {
  -- Constructs an empty result

  -- RuBisCO capacity (gC/m2/day)
  vm : real,

  -- gross daily photosynthesis (gC/m2/day)
  agd_g : real,

  -- leaf-level net daytime photosynthesis
  -- expressed in CO2 diffusion units (mm/m2/day)--
  adtmm : real,

  -- leaf respiration (gC/m2/day)
  rd_g : real,

  -- PAR-limited photosynthesis rate (gC/m2/h)
  je : real,

  -- optimal leaf nitrogen associated with photosynthesis (kgN/m2)
  nactive_opt : real,

  -- nitrogen limitation on vm
  vmaxnlim : real

  -- net C-assimilation (gross photosynthesis minus leaf respiration) (kgC/m2/day)
}


--- Class containing serializable variables for Weathergenerator GWGen
---
--Variables for the build-in random-generator and to keep track of whether past days
--were rain days

let QSIZ : int = 10
type WeatherGenState = {
  --- Random state variable q
  q : [QSIZ]int,
  --- Random state variable carry
  carry : int,
  --- Random state variable xcng
  xcng : int,
  --- Random state variable xs
  xs : uint,
  --- Random state variable indx
  indx : int,
  --- Random state variable have
  have : bool,
  --- Random state gamma
  gamma_vals : [2]real,
  --- Indicator for whether the recent two days were rein-days
  pday : [2]bool,
  --- Random state's residuals
  resid : [4]real
}



-- This struct contains the environmental input to a photosynthesis calculation.
-- \see photosynthesis--
type PhotosynthesisEnvironment = {
  -- atmospheric ambient CO2 concentration (ppmv)
  co2 : real,
  -- mean air temperature today (deg C)
  temp : real,
  -- total daily photosynthetically-active radiation today (J / m2 / day) (ALPHAA not yet accounted for)
  par : real,
  -- fraction of PAR absorbed by foliage
  fpar : real,
  -- day length, must equal 24 in diurnal mode(h)
  daylength : real
}

--- This struct contains the stresses used in a photosynthesis calculation.
type PhotosynthesisStresses = {
  --- whether nitrogen should limit Vmax
  ifnlimvmax: bool,

  ---  limit to moss photosynthesis. [0,1], where 1 means no limit
  moss_ps_limit: real,

  --- limit to graminoid photosynthesis. [0,1], where 1 means no limit
  graminoid_ps_limit: real,

  --- limit to photosynthesis due to inundation, where 1 means no limit
 inund_stress: real
}




--- The Climate for a grid cell
--- Stores all static and variable data relating to climate parameters, as well as
--  latitude, atmospheric CO2 concentration and daylength for a grid cell. Includes
--  a reference to the parent Gridcell object (defined below). Initialised by a
--  call to initdrivers.
---
type Climate = {

  -- MEMBER VARIABLES
  --- reference to parent Gridcell object
  climate_id: int,

  gridcell_id: int,

  --- values for randomisation in Weathergenerator GWGEN
  weathergenstate : WeatherGenState,

  --- mean air temperature today (deg C)
  temp : real,

  --- total daily net downward shortwave solar radiation today (J/m2/day)
  rad : real,

  --- total daily photosynthetically-active radiation today (J/m2/day)
  par : real,

  --- precipitation today (mm)
  prec : real,

  --- 10 m wind (km/h)
  u10 : real,

  --- rel. humidity (fract.)
  relhum : real,

  --- min and max daily temperature (deg C)
  tmin : real,
  tmax : real,

  --- day length today (h)
  daylength : real,

  --- atmospheric ambient CO2 concentration today (ppmv)
  co2 : real,

  --- latitude (degrees : real, +=north, -=south)
  lat : real,

  --- Insolation today, see also instype
  insol : real,

  --- Type of insolation
  --- This decides how to interpret the variable insol,
  --  see also documentation for the insoltype enum.
  instype : insoltype,

  --- equilibrium evapotranspiration today (mm/day)
  eet : real,

  --- mean temperature for the last 31 days (deg C)
  mtemp : real,

  --- mean of lowest mean monthly temperature for the last 20 years (deg C)
  mtemp_min20 : real,

  --- mean of highest mean monthly temperature for the last 20 years (deg C)
  mtemp_max20 : real,

  --- highest mean monthly temperature for the last 12 months (deg C)
  mtemp_max : real,

  --- accumulated growing degree day sum on 5 degree base
  --- reset when temperatures fall below 5 deg C */
  gdd5 : real,

  --- total gdd5 (accumulated) for this year (reset 1 January)
  agdd5 : real,

  --- accumulated growing degree day sum on 0 degree base (Wolf et al. 2008)
  gdd0 : real,

  --- total gdd0 (accumulated) for this year (reset 1 January)
  agdd0 : real,

  --- total gdd0 (accumulated) over each of the last 20 years
  --- climate.agdd0_20.mean() gives the average total gdd0 (accumulated) over the last 20 years
  --Historic<double, 20> agdd0_20 : real,

  --- number of days with temperatures <5 deg C
  --- reset when temperatures fall below 5 deg C,
  --  maximum value is number of days in the year */
  chilldays : int,

  --- true if chill day count may be reset by temperature fall below 5 deg C
  ifsensechill : bool,

  --- Respiration response to today's air temperature incorporating damping of Q10
  --  due to temperature acclimation (Lloyd & Taylor 1994)
  --/
  gtemp : real,

  --- daily temperatures for the last 31 days (deg C)
  --Historic<double, 31> dtemp_31 : real,

  --- daily precipitation for the last 31 days (deg C)
  --Historic<double, 31> dprec_31 : real,

  --- daily eet for the last 31 days (deg C)
  --Historic<double, 31> deet_31 : real,

  --- minimum monthly temperatures for the last 20 years (deg C)
  mtemp_min_20 : [20]real,

  --- maximum monthly temperatures for the last 20 years (deg C)
  mtemp_max_20 : [20]real,

  --- minimum monthly temperature for the last 12 months (deg C)
  mtemp_min : real,

  --- mean of monthly temperatures for the last 12 months (deg C)
  atemp_mean : real,


  -- BLAZE
  --- average annual rainfall (mm/a)
  avg_annual_rainfall : real,
  --- current sum of annual Rainfall (mm)
  cur_rainfall : real,
  --- Accumulated last rainfall (mm)
  last_rainfall : real,
  --- Days since last rainfall
  days_since_last_rainfall : real,
  --- Keetch-Byram-Drought-Index
  kbdi : real,
  --- McArthur forest fire index (FFDI)
  mcarthur_forest_fire_index : real,
  --- To keep track of running months daily FFDI
  months_ffdi : [30]real,

  -- Saved parameters used by function daylengthinsoleet

  sinelat : real,
  cosinelat : real,
  qo : [Date_MAX_YEAR_LENGTH]real,
  u : [Date_MAX_YEAR_LENGTH]real,
  v : [Date_MAX_YEAR_LENGTH]real,
  hh : [Date_MAX_YEAR_LENGTH]real,
  sinehh : [Date_MAX_YEAR_LENGTH]real,
  daylength_save : [Date_MAX_YEAR_LENGTH]real,
  --- indicates whether saved values exist for this day
  doneday : [Date_MAX_YEAR_LENGTH]bool,

  --- diurnal temperature range, used in daily/monthly BVOC (deg C)
  dtr : real,

  -- containers for sub-daily values of temperature, short-wave downward
  -- radiation, par, rad and gtemp (equivalent to temp, insol, par, rad and gtemp)
  -- NB: units of these variable are the same as their daily counterparts,
  -- i.e. representing daily averages (e.g. pars [J/m2/day])

  -- FIXME: these are never initialized, and never used?
  --- Sub-daily temperature (deg C) (\see temp)
  temps : [Date_subdaily]real,

  --- Sub-daily insolation (\see insol)
  insols : [Date_subdaily]real,

  --- Sub-daily PAR (\see par)
  pars : [Date_subdaily]real,

  --- Sub-daily net downward shortwave solar radiation (\see rad)
  rads : [Date_subdaily]real,

  --- Sub-daily respiration response (\see gtemp)
  gtemps : [Date_subdaily]real,

  --- Variables used for crop sowing date or seasonality calculation

  --- daily precipitations for the last 10 days (mm)
  dprec_10 : [10]real,
  --- daily 10 day-sums of precipitations for today and yesterday (mm)
  sprec_2 : [2]real,
  --- max temperature during the last test period
  maxtemp : real,
  --- summer day when we test last year's crossing of sowing temperature limits : real, NH:June 30(day 180), SH:Dec.31(day 364), set in getgridcell()
  testday_temp : int,
  --- last day of dry month when we test last year's crossing of sowing precipitation limits : real, NH:Dec.31(day 364), SH:June 30(day 180), set in getgridcell()
  testday_prec : int,
  --- date used for sowing if no frost or spring occured during the year between the testmonths : real, NH:14, SH:195, set in getgridcell()
  coldestday : int,
  --- used to adapt equations to hemisphere, set in getgridcell()
  adjustlat : int,
  --- accumulated monthly pet values for this year
  mpet_year : [12]real,
  --- past 20 years monthly temperature values
  mtemp_20 : [20][12]real,
  --- past 20 years monthly precipitation values
  mprec_20 : [20][12]real,
  --- past 20 years monthly PET values
  mpet_20 : [20][12]real,
  --- past 20 years monthly precipitation to PET ratios
  mprec_pet_20 : [20][12]real,
  --- past 20 years minimum of monthly precipitation to PET ratios
  mprec_petmin_20 : [20]real,
  --- past 20 years maximum of monthly precipitation to PET ratios
  mprec_petmax_20 : [20]real,
  --- 20-year running average monthly temperature values
  mtemp20 : [12]real,
  --- 20-year running average monthly precipitation values
  mprec20 : [12]real,
  --- 20-year running average monthly PET values
  mpet20 : [12]real,
  --- 20-year running average monthly precipitation to PET ratios
  mprec_pet20 : [12]real,
  --- 20-year running average of minimum monthly precipitation to PET ratios
  mprec_petmin20 : real,
  --- 20-year running average of maximum monthly precipitation to PET ratios
  mprec_petmax20 : real,

  --Historic<double, 20> hmtemp_20[12] : real,
  --Historic<double, 20> hmprec_20[12] : real,
  --Historic<double, 20> hmeet_20[12] : real,

  --- seasonality type (SEASONALITY_NO, SEASONALITY_PREC, SEASONALITY_PRECTEMP, SEASONALITY_TEMP, SEASONALITY_TEMPPREC)
  seasonality : seasonality_type,
  seasonality_lastyear : seasonality_type,

  --- precipitation seasonality type (DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET)
  --- based on the extremes of the 20-year monthly means
  --/
  prec_seasonality : prec_seasonality_type,
  prec_seasonality_lastyear : prec_seasonality_type,

  --- precipitation range (DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET)
  --- based on the average of the 20-year monthly extremes
  --/
  prec_range : prec_seasonality_type,
  prec_range_lastyear : prec_seasonality_type,

  --- temperature seasonality (COLD, COLD_WARM, COLD_HOT, WARM, WARM_HOT, HOT)
  temp_seasonality : temp_seasonality_type,
  temp_seasonality_lastyear : temp_seasonality_type,

  --- whether several months with precipitation maxima exists (remains to be implemented)
  biseasonal : bool,

  --- variation coefficient of 20-year mean monthly temperatures
  var_prec : real,
  --- variation coefficient of 20-year mean monthly precipitation to PET ratios
  var_temp : real,

  --- annual precipitation sum
  aprec : real,
  --- annual average precipitation (last year) (mm)
  aprec_lastyear : real
}




--- Stores accumulated monthly and annual fluxes.
-- This class handles the storage and accounting of fluxes for a single patch.
--  Different fluxes can be stored in different ways, depending on what kind of
--  flux it is and what kind of output we want. The details of whether fluxes
--  are stored per PFT or just as a patch total, or per day, month or only a
--  yearly sum, is hidden from the 'scientific' code, which merely reports the
--  fluxes generated.
--/
type PerPatchFluxType = enum_type
type PerPFTFluxType = enum_type


--- Fluxes stored as totals for the whole patch
--- Carbon flux to atmosphere from burnt vegetation and litter (kgC/m2)
let FIREC : PerPatchFluxType = 0
--- Carbon flux to atmosphere from soil respiration (kgC/m2)
let SOILC : PerPatchFluxType = 1
--- Flux from atmosphere to vegetation associated with establishment (kgC/m2)
let ESTC : PerPatchFluxType = 2
--- Flux to atmosphere from consumed harvested products (kgC/m2)
let HARVESTC : PerPatchFluxType = 3
--- Flux from atmosphere to vegetation associated with sowing (kgC/m2)
let SEEDC : PerPatchFluxType = 4
--- Flux from atmosphere to vegetation associated with manure addition (kgC/m2)
let MANUREC : PerPatchFluxType = 5
--- Flux to vegetation associated with manure addition (kgN/m2)
let MANUREN : PerPatchFluxType = 6
--- Flux to vegetation associated with N addition (kgN/m2)
let NFERT : PerPatchFluxType = 7
--- Nitrogen flux to atmosphere from consumed harvested products (kgN/m2)
let HARVESTN : PerPatchFluxType = 8
--- Nitrogen flux from atmosphere to vegetation associated with sowing (kgC/m2)
let SEEDN : PerPatchFluxType = 9
--- NH3 flux to atmosphere from fire
let NH3_FIRE : PerPatchFluxType = 10
--- NOx flux to atmosphere from fire
let NOx_FIRE : PerPatchFluxType = 11
--- N2O flux to atmosphere from fire
let N2O_FIRE : PerPatchFluxType = 12
--- N2 flux to atmosphere from fire
let N2_FIRE : PerPatchFluxType = 13
------ Soil N transformation -----
--- NH3 flux from soil (ntransform)
let NH3_SOIL : PerPatchFluxType = 14
--- NO flux from soil (ntransform)
let NO_SOIL : PerPatchFluxType = 15
--- N2O flux in soil (ntransform)
let N2O_SOIL : PerPatchFluxType = 16
--- N2 flux from soil (ntransform)
let N2_SOIL : PerPatchFluxType = 17
--- DOC flux from soil (ntransform)
let DOC_FLUX : PerPatchFluxType = 18
--- Net nitrification (ntransform)
let NET_NITRIF : PerPatchFluxType = 19
--- Net denitrification (ntransform)
let NET_DENITRIF : PerPatchFluxType = 20
--- Gross nitrification (ntransform)
let GROSS_NITRIF : PerPatchFluxType = 21
--- Gross denitrification (ntransform)
let GROSS_DENITRIF : PerPatchFluxType = 22
--- Reproduction costs
let REPRC : PerPatchFluxType = 23
--- Total (i.e. CH4C_DIFF + CH4C_PLAN + CH4C_EBUL) CH4 flux to atmosphere from peatland soils (gC/m2).
let CH4C : PerPatchFluxType = 24
--- Diffused CH4 flux to atmosphere from peatland soils (gC/m2).
let CH4C_DIFF : PerPatchFluxType = 25
--- Plant-mediated CH4 flux to atmosphere from peatland soils (gC/m2).
let CH4C_PLAN : PerPatchFluxType = 26
--- CH4 flux to atmosphere from peatland soils due to ebullition (gC/m2).
let CH4C_EBUL : PerPatchFluxType = 27
--- Number of types must be last
let NPERPATCHFLUXTYPES : PerPatchFluxType = 28

--- Fluxes stored per pft
--- NPP (kgC/m2)
let NPP : PerPFTFluxType = 0
--- GPP (kgC/m2)
let GPP : PerPFTFluxType = 1
--- Autotrophic respiration (kgC/m2)
let RA : PerPFTFluxType = 2
--- Isoprene (mgC/m2)
let ISO : PerPFTFluxType = 3
--- Monoterpene (mgC/m2)
let MT_APIN : PerPFTFluxType = 4
let MT_BPIN : PerPFTFluxType = 5
let MT_LIMO : PerPFTFluxType = 6
let MT_MYRC : PerPFTFluxType = 7
let MT_SABI : PerPFTFluxType = 8
let MT_CAMP : PerPFTFluxType = 9
let MT_TRIC : PerPFTFluxType = 10
let MT_TBOC : PerPFTFluxType = 11
let MT_OTHR : PerPFTFluxType = 12
--- Number of types must be last
let NPERPFTFLUXTYPES : PerPFTFluxType = 13

-- emission ratios from fire (NH3, NOx, N2O, N2) Delmas et al. 1995

let NH3_FIRERATIO : real = 0.005
let NOx_FIRERATIO : real = 0.237
let N2O_FIRERATIO : real = 0.036
let N2_FIRERATIO  : real = 0.722

--- Reference to patch to which this Fluxes object belongs
--Patch& patch, --FUTHARK: ID instead of pointer? What is it used for?


type Fluxes = {
  ---- Private members, mutable
  -- Stores one flux value per PFT and flux type
  annual_fluxes_per_pft : [npft][NPERPFTFLUXTYPES]real,

  -- Stores one flux value per month and flux type
  -- For the fluxes only stored as totals for the whole patch--
  monthly_fluxes_patch : [12][NPERPATCHFLUXTYPES]real,

  -- Stores one flux value per month and flux type
  --- For the fluxes stored per pft for annual values--
  monthly_fluxes_pft : [12][NPERPFTFLUXTYPES]real,

  -- Stores one flux value per day and flux type
  daily_fluxes_patch : [365][NPERPATCHFLUXTYPES]real,

  -- Stores one flux value per day and flux type
  daily_fluxes_pft : [365][NPERPFTFLUXTYPES]real

}


type planting_system_type = enum_type
let planting_system_NONE : planting_system_type = 0
let MONOCULTURE : planting_system_type = 1
let SELECTION : planting_system_type = 2

type harvest_system_type = enum_type
let harvest_system_NONE : harvest_system_type = 0
let CLEARCUT : harvest_system_type = 1
let CONTINUOUS : harvest_system_type = 2


--- Storage class of crop management information for one rotation period for a stand type, read from the instruction file.
type ManagementType = {
  --- id code (should be zero based and sequential, 0...nst-1)
  managementtype_id : int,
  --- name of management type
  --xtring name,

  --- type of planting system ("", "MONOCULTURE", "SELECTION", etc.)
  planting_system : planting_system_type,
  --- type of harvest system ("", "CLEARCUT", "CONTINUOUS")
  harvest_system : harvest_system_type,
  --- name of crop pft
  --xtring pftname,
  --- identifier of pft selection
  --xtring selection,
  --- Rotation period in years
  nyears : real,
  --- hydrology (RAINFED,IRRIGATED)
  hydrology : hydrologytype,
  --- irrigation efficiency
--  double firr,
  --- forced sowing date, unless sdate_force read from file
  sdate : int,
  --- forced harvest date, unless hdate_force read from file
  hdate : int,
  --- Nitrogen fertilisation amount, unless Nfert_read read from file
  nfert : real,
  --- Whether grass is grown in fallow
  fallow : bool,
  --- Double cropping of one crop (e.g. rice)
  multicrop : bool
}






type CropRotation = {
  ncrops : int,
  firstrotyear : int
}

type naturalvegType = enum_type
let naturalvegNONE : naturalvegType = 0
let GRASSONLY : naturalvegType = 1
let naturalvegALL : naturalvegType = 2

type reestabType = enum_type
let reestabNONE : reestabType = 0
let RESTRICTED : reestabType = 1
let reestabALL : reestabType = 2

--- Stand type class for storing both static parameters, read from the instruction file,
--  and dynamic variables, updated in landcover_change()
--  Active stand types are stored in the stlist analogous to the pftlist.
type StandType = {
  --- id code (should be zero based and sequential, 0...nst-1)
  standtype_id : int,
  --- name of stand type
  --name : xtring,

  --- specifies type of landcover
  -- \see landcovertype */
  landcover : landcovertype,  -- specifies type of landcover (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND)
  --- Rotation information, read from the instruction file
  rotation : CropRotation,
  --- Management struct (static)
  management : ManagementType,
  --- Management types in a rotation cycle
  --xtring mtnames[NROTATIONPERIODS_MAX],
  --- First management year: sets time when common features for managed stands begin, e.g. relaxed establishment rules and absence of disturbance before harvest begins
  -- \this currently only applies for stands with wood havest */
  firstmanageyear : int,

  --- intercrop (NOINTERCROP,NATURALGRASS)
  intercrop : intercroptype,
  --- whether natural pft:s are allowed to grow in stand type
  naturalveg : naturalvegType, -- "", "GRASSONLY", "ALL"
  -- whether only pft:s defined in management are allowed (plus intercrop or naturalveg/grass)
  restrictpfts : bool,
  --- whether planted pft:s or all active pft:s are allowed to established after planting in a forest stand ("", "RESTRICTED", "ALL")
  reestab : reestabType
}


let nst : i64 = 5 -- num possible StandTypes

-- Holds static functional parameters for a plant functional type (PFT).
-- There should be one Pft object for each potentially occurring PFT. The same Pft object
-- may be referenced (via the pft member of the Individual object see below) by different
-- average individuals. Member functions are included for initialising SLA given leaf
-- longevity, and for initialising sapling/regen characteristics (required for
-- population mode).

type Pft  = {
  -- MEMBER VARIABLES
  -- id code (should be zero based and sequential, 0...npft-1)
  pft_id: int,
  -- name of PFT
  name: xtring,
  -- life form (tree or grass)
  lifeform: lifeformtype,
  -- leaf phenology (raingreen, summergreen, evergreen, rain+summergreen, cropgreen)
  phenology: phenologytype,
  -- leaf physiognomy (needleleaf, broadleaf)
  leafphysiognomy: leafphysiognomytype,
  -- growing degree sum on 5 degree base required for full leaf cover
  phengdd5ramp: real,
  -- water stress threshold for leaf abscission (range 0-1 raingreen PFTs)
  wscal_min: real,
  -- biochemical pathway for photosynthesis (C3 or C4)
  pathway: pathwaytype,
  -- approximate low temperature limit for photosynthesis (deg C)
  pstemp_min: real,
  -- approximate lower range of temperature optimum for photosynthesis (deg C)
  pstemp_low: real,
  -- approximate upper range of temperature optimum for photosynthesis (deg C)
  pstemp_high: real,
  -- maximum temperature limit for photosynthesis (deg C)
  pstemp_max: real,
  -- non-water-stressed ratio of intercellular to ambient CO2 partial pressure
  lambda_max: real,
  -- vegetation root profile in an array containing fraction of roots in each soil layer, [0=upper layer]
  rootdist: [NSOILLAYER]real,
  -- shape parameter for initialisation of root distribtion
  root_beta: real,
  -- canopy conductance component not associated with photosynthesis (mm/s)
  gmin: real,
  -- maximum evapotranspiration rate (mm/day)
  emax: real,
  -- maintenance respiration coefficient (0-1)
  respcoeff: real,

  -- minimum leaf C:N mass ratio allowed when nitrogen demand is determined
  cton_leaf_min: real,
  -- maximum leaf C:N mass ratio  allowed when nitrogen demand is determined
  cton_leaf_max: real,
  -- average leaf C:N mass ratio (between min and max)
  cton_leaf_avr: real,
  -- average fine root C:N mass ratio (connected cton_leaf_avr)
  cton_root_avr: real,
  -- maximum fine root C:N mass ratio (used when mass is negligible)
  cton_root_max: real,
  -- average sapwood C:N mass ratio (connected cton_leaf_avr)
  cton_sap_avr: real,
  -- maximum sapwood C:N mass ratio (used when mass is negligible)
  cton_sap_max: real,
  -- reference fine root C:N mass ratio
  cton_root: real,
  -- reference sapwood C:N mass ratio
  cton_sap: real,
  -- Maximum nitrogen (NH4+ and NO3- seperatly) uptake per fine root [kgN kgC-1 day-1]
  nuptoroot: real,
  -- coefficient to compensate for vertical distribution of fine root on nitrogen uptake
  nupscoeff: real,
  -- fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen longterm storage scalar
  fnstorage: real,

  -- Michaelis-Menten kinetic parameters
  -- Half saturation concentration for N uptake [kgN l-1] (Rothstein 2000)--
  km_volume: real,

  -- fraction of NPP allocated to reproduction
  reprfrac: real,
  -- annual leaf turnover as a proportion of leaf C biomass
  turnover_leaf: real,
  -- annual fine root turnover as a proportion of fine root C biomass
  turnover_root: real,
  -- annual sapwood turnover as a proportion of sapwood C biomass
  turnover_sap: real,
  -- sapwood and heartwood density (kgC/m3)
  wooddens: real,
  -- maximum tree crown area (m2)
  crownarea_max: real,
  -- constant in allometry equations
  k_allom1: real,
  -- constant in allometry equations
  k_allom2: real,
  -- constant in allometry equations
  k_allom3: real,
  -- constant in allometry equations
  k_rp: real,
  -- tree leaf to sapwood area ratio
  k_latosa: real,
  -- specific leaf area (m2/kgC)
  sla: real,
  -- leaf longevity (years)
  leaflong: real,
  -- leaf to root mass ratio under non-water-stressed conditions
  ltor_max: real,
  -- litter moisture flammability threshold (fraction of AWC)
  litterme: real,
  -- fire resistance (0-1)
  fireresist: real,
  -- minimum forest-floor PAR level for growth (grasses) or establishment (trees)
  -- J/m2/day, individual and cohort modes--
  parff_min: real,
  -- parameter capturing non-linearity in recruitment rate relative to
  -- understorey growing conditions for trees (Fulton 1991) (individual and
  -- cohort modes)

  alphar: real,
  -- maximum sapling establishment rate (saplings/m2/year) (individual and cohort modes)
  est_max: real,
  -- constant used in calculation of sapling establishment rate when spatial
  -- mass effect enabled (individual and cohort modes)

  kest_repr: real,
  -- constant affecting amount of background establishment
  -- \see ifbgestab--
  kest_bg: real,
  -- constant used in calculation of sapling establishment rate when spatial
  -- mass effect disabled (individual and cohort modes)

  kest_pres: real,
  -- expected longevity under non-stressed conditions (individual and cohort modes)
  longevity: real,
  -- threshold growth efficiency for imposition of growth suppression mortality
  -- kgC/m2 leaf/year, individual and cohort modes--
  greff_min: real,

  -- Bioclimatic limits (all temperatures deg C)

  -- minimum 20-year coldest month mean temperature for survival
  tcmin_surv: real,
  -- maximum 20-year coldest month mean temperature for establishment
  tcmax_est: real,
  -- minimum degree day sum on 5 deg C base for establishment
  gdd5min_est: real,
  -- minimum 20-year coldest month mean temperature for establishment
  tcmin_est: real,
  -- minimum warmest month mean temperature for establishment
  twmin_est: real,
  -- continentality parameter for boreal summergreen trees
  twminusc: real,
  -- constant in equation for budburst chilling time requirement (Sykes et al 1996)
  k_chilla: real,
  -- coefficient in equation for budburst chilling time requirement
  k_chillb: real,
  -- exponent in equation for budburst chilling time requirement
  k_chillk: real,
  -- array containing values for GDD0(c) given c=number of chill days
  -- Sykes et al 1996, Eqn 1
  -- gdd0 has one element for each possible value for number of chill days

  ---- FIXME TODO HACK! comes from Date_MAX_YEAR_LENGTH+1
  gdd0: [Date_MAX_YEAR_LENGTH_plusone]real,

  -- interception coefficient (unitless)
  intc: real,

  -- the amount of N that is applied (kg N m-2)
  N_appfert: real,
  -- 0 - 1 how much of the fertiliser is applied the first date, default 1.
  fertrate: (real, real),
  -- dates relative to sowing date
  fertdates: (int, int),
  fert_stages: (real, real),
  fertilised: (bool, bool),

  T_vn_min: real,
  T_vn_opt: real,
  T_vn_max: real,

  T_veg_min: real,
  T_veg_opt: real,
  T_veg_max: real,

  T_rep_min: real,
  T_rep_opt: real,
  T_rep_max: real,

  photo: (real, real, real),

  dev_rate_veg: real,
  dev_rate_rep: real,

  a1: real, b1: real, c1: real, d1: real, a2: real, b2: real, c2: real, d2: real, a3: real, b3: real, c3: real, d3: real,
  cton_stem_avr: real,
  cton_stem_max: real,

  -- Drought tolerance level (0 = very -> 1 = not at all) (unitless)
  -- Used to implement drought-limited establishment--
  drought_tolerance: real,

  -- bvoc

  -- aerodynamic conductance (m s-1)
  ga: real,
  -- isoprene emission capacity (ug C g-1 h-1)
  eps_iso: real,
  -- whether (1) or not (1) isoprene emissions show a seasonality
  seas_iso: bool,
  -- monoterpene emission capacity (ug C g-1 h-1) per monoterpene species
  eps_mon: [NMTCOMPOUNDS]real,
  -- fraction of monoterpene production that goes into storage pool (-) per monoterpene species
  storfrac_mon: [NMTCOMPOUNDS]real,

  -- Bioclimatic limits parameters from Wolf et al. 2008

  -- snow max [mm]
  max_snow: real,
  -- snow min [mm]
  min_snow: real,
  -- GDD0 min
  gdd0_min: real,
  -- GDD0 max
  gdd0_max: real,

  -- New parameters from parameters from Wania et al. (2009a, 2009b, 2010)

  -- Days per month for which inundation is tolerated
  inund_duration: int,
  -- Inundation stress is felt when the water table (mm) is above wtp_max
  wtp_max: real,
  -- Whether this PFT has aerenchyma through which O2 and CH4 can be transported (Wania et al. 2010 - Sec 2.6)
  has_aerenchyma: bool,

  -- Sapling/regeneration characteristics (used only in population mode)
  -- For trees, on sapling individual basis (kgC) for grasses, on stand area basis,
  -- kgC/m2--

  -- leaf C biomass
  regen_cmass_leaf: real,
  -- fine root C biomass
  regen_cmass_root: real,
  -- sapwood C biomass
  regen_cmass_sap: real,
  -- heartwood C biomass
  regen_cmass_heart: real,

  -- specifies type of landcover pft is allowed to grow in (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND)
  landcover: landcovertype,
  -- pft selection
  selection: xtring,
  -- fraction of residue outtake at harvest
  res_outtake: real,
  -- harvest efficiencytype leafphysiognomy = int

  harv_eff: real,
  -- harvest efficiency of intercrop grass
  harv_eff_ic: real,
  -- fraction of harvested products that goes into patchpft.harvested_products_slow
  harvest_slow_frac: real,
  -- yearly turnover fraction of patchpft.harvested_products_slow (goes to gridcell.acflux_harvest_slow)
  turnover_harv_prod: real,
  -- whether pft may grow as cover crop
  isintercropgrass: bool,
  -- whether autumn temperature dependent sowing date is calculated
  ifsdautumn: bool,
  -- upper temperature limit for autumn sowing
  tempautumn: real,
  -- lower temperature limit for spring sowing
  tempspring: real,
  -- default length of growing period
  lgp_def: int,
  -- upper minimum temperature limit for crop sowing
  maxtemp_sowing: real,
  -- default sowing date in the northern hemisphere (julian day)
  sdatenh: int,
  -- default sowing date in the southern hemisphere
  sdatesh: int,
  -- whether sowing date adjusting equation is used
  sd_adjust: bool,
  -- parameter 1 in sowing date adjusting equation
  sd_adjust_par1: real,
  -- parameter 2 in sowing date adjusting equation
  sd_adjust_par2: real,
  -- parameter 3 in sowing date adjusting equation
  sd_adjust_par3: real,
  -- latest date for harvesting in the northern hemisphere
  hlimitdatenh: int,
  -- latest date for harvesting in the southern hemisphere
  hlimitdatesh: int,
  -- default base temperature (°C) for heat unit (hu) calculation
  tb: real,
  -- temperature under which vernalisation is possible (°C)
  trg: real,
  -- default number of vernalising days required
  pvd: int,
    -- sensitivity to the photoperiod effect [0-1]
  psens: real,
  -- basal photoperiod (h) (pb<ps for longer days plants)
  pb: real,
  -- lag in days after sowing before vernalization starts
  vern_lag: int,
  -- saturating photoperiod (h) (ps<pb for shorter days plants)
  ps: real,
  -- default potential heat units required for crop maturity (degree-days)
  phu: real,
  -- whether quadratic equation used for calculating potential heat units (Bondeau method)
  phu_calc_quad: bool,
  -- whether linear equation used for calculating potential heat units (Bondeau method)
  phu_calc_lin: bool,
  -- minimum potential heat units required for crop maturity (Bondeau method) (degree-days)
  phu_min: real,
  -- maximum potential heat units required for crop maturity (Bondeau method) (degree-days)
  phu_max: real,
  -- reduction factor of potential heat units in spring crops (Bondeau method) (degree-days)
  phu_red_spring_sow: real,
  -- number of days of phu decrease in the linear phu equation (Bondeau method)
  ndays_ramp_phu: real,
  -- intercept for the linear phu equation (Bondeau method)
  phu_interc: real,
  -- fraction of growing season (phu) at which senescence starts [0-1]
  fphusen: real,
  -- type of senescence curve (see Bondeau et al. 2007)
  shapesenescencenorm: bool,
  -- fraction of maximal LAI still present at harvest [0-1]
  flaimaxharvest: real,
  -- default maximum LAI (only used for intercrop grass in the case where no pasture is present in any stand)
  laimax: real,
  -- whether harvestable organs are above ground
  aboveground_ho: bool,
  -- optimum harvest index
  hiopt: real,
  -- minimum harvest index
  himin: real,
  -- initial fraction of growing season's npp allocated to roots
  frootstart: real,
  -- final fraction of growing season's npp allocated to roots
  frootend: real,
  -- autumn/spring sowing of pft:s with tempautumn = 1
  forceautumnsowing: int,  --0 = NOFORCING,  1 = AUTUMNSOWING, 2 = SPRINGSOWING
  -- N limited version of pft
  nlim: bool
}

-- type cropindiv

--/ Container for crop-specific data at the individual level
type cropindiv_struct = {

  --Plant carbon biomass variables are all on patch area basis (kgC/m2)

  --/ year's harvestable organ C biomass (= ycmass_plant)
  cmass_ho : real,
  --/ above-ground pool C biomass (when calculating daily cmass_leaf from lai_crop) (= ycmass_agpool)
  cmass_agpool : real,
  cmass_stem : real,

  --/ year's maximum value of leaf C biomass
  cmass_leaf_max : real,
  --/ grs_cmass_leaf value saved at day before senescence (for LAI-calculation in allometry)
  cmass_leaf_sen : real,

  --/ today's increase of whole plant C biomass
  dcmass_plant : real,
  --/ today's increase of leaf C biomass
  dcmass_leaf : real,
  --/ today's increase of root C biomass
  dcmass_root : real,
  --/ today's increase of harvestable organ C biomass
  dcmass_ho : real,
  --/ today's increase of above-ground pool C biomass
  dcmass_agpool : real,
  dcmass_stem : real,

  --/ today's increase of leaf N biomass
  dnmass_leaf : real,
  --/ today's increase of root N biomass
  dnmass_root : real,
  --/ today's increase of harvestable organ N biomass
  dnmass_ho : real,
  --/ today's increase of above-ground pool N biomass
  dnmass_agpool : real,

  --/CARBON
  --/ daily updated whole plant C biomass, reset at harvest day
  grs_cmass_plant : real,
  --/ daily updated leaf C biomass, reset at harvest day
  grs_cmass_leaf : real,
  --/ daily updated root C biomass, reset at harvest day
  grs_cmass_root : real,
  --/ daily updated harvestable organ C biomass, reset at harvest day
  grs_cmass_ho : real,
  --/ daily updated above-ground pool C biomass, reset at harvest day
  grs_cmass_agpool : real,
  --/ daily updated dead leaf C biomass, reset at harvest day
  grs_cmass_dead_leaf : real,
  --/ daily updated stem pool C biomass, reset at harvest day
  grs_cmass_stem : real,

  --/ carbon content of harvestable organs saved on first day of land use change year
  grs_cmass_leaf_luc : real,
  --/ carbon content of harvestable organs saved on first day of land use change year
  grs_cmass_root_luc : real,
  --/ carbon content of harvestable organs saved on first day of land use change year
  grs_cmass_ho_luc : real,
  --/ carbon content of above-ground pool saved on first day of land use change year
  grs_cmass_agpool_luc : real,
  --/ carbon content of dead leaves saved on first day of land use change year
  grs_cmass_dead_leaf_luc : real,
  --/ carbon content of stem saved on first day of land use change year
  grs_cmass_stem_luc : real,

  --/ daily updated whole plant C biomass, reset at day 0
  ycmass_plant : real,
  --/ daily updated leaf C biomass, reset at day 0
  ycmass_leaf : real,
  --/ daily updated root C biomass, reset at day 0
  ycmass_root : real,
  --/ daily updated harvestable organ C biomass, reset at day 0
  ycmass_ho : real,
  --/ daily updated above-ground pool C biomass, reset at day 0
  ycmass_agpool : real,
  --/ daily updated dead leaf C biomass, reset at day 0
  ycmass_dead_leaf : real,
  --/ daily updated stem C biomass, reset at day 0
  ycmass_stem : real,

  --/ year's whole plant C biomass at time of harvest (cumulative if several harvest events)
  harv_cmass_plant : real,
  --/ year's leaf C biomass at time of harvest (cumulative if several harvest events)
  harv_cmass_leaf : real,
  --/ year's root C biomass at time of harvest (cumulative if several harvest events)
  harv_cmass_root : real,
  --/ year's harvestable organ C biomass at time of harvest (cumulative if several harvest events)
  harv_cmass_ho : real,
  --/ year's above-ground pool C biomass at time of harvest (cumulative if several harvest events)
  harv_cmass_agpool : real,
  --/ year's stem C biomass at time of harvest (cumulative if several harvest events)
  harv_cmass_stem : real,

  --/NITROGEN
  --/ nitrogen content of harvestable organs
  nmass_ho : real,
  --/ nitrogen content of above-ground pool
  nmass_agpool : real,
  --/ nitrogen content of dead leaves
  nmass_dead_leaf : real,

  --/ nitrogen content of harvestable organs saved on first day of land use change year
  nmass_ho_luc : real,
  --/ nitrogen content of above-ground pool saved on first day of land use change year
  nmass_agpool_luc : real,
  --/ nitrogen content of dead leaves saved on first day of land use change year
  nmass_dead_leaf_luc : real,

  --/ daily updated leaf N biomass, reset at day 0
  ynmass_leaf : real,
  --/ daily updated root N biomass, reset at day 0
  ynmass_root : real,
  --/ daily updated harvestable organ N biomass, reset at day 0
  ynmass_ho : real,
  --/ daily updated above-ground pool N biomass, reset at day 0
  ynmass_agpool : real,
  --/ daily updated dead leaf N biomass, reset at day 0
  ynmass_dead_leaf : real,

  --/ year's leaf N biomass at time of harvest (cumulative if several harvest events)
  harv_nmass_leaf : real,
  --/ year's root N biomass at time of harvest (cumulative if several harvest events)
  harv_nmass_root : real,
  --/ year's harvestable organ N biomass at time of harvest (cumulative if several harvest events)
  harv_nmass_ho : real,
  --/ year's above-ground pool N biomass at time of harvest (cumulative if several harvest events)
  harv_nmass_agpool : real,

  --/ dry weight crop yield harvested this year (cumulative if several harvest events), based on harv_cmass_xx
  harv_yield : real,

  --/ harvestable organ C biomass at the last two harvest events this year
  cmass_ho_harvest : [2]real,
  --/ harvestable organ N biomass at the last two harvest events this year
  nmass_ho_harvest : [2]real,
  --/ dry weight crop yield at the last two harvest events this year
  yield_harvest : [2]real,

  --/ dry weight crop yield grown this year (cumulative if several harvest events), based on ycmass_xx
  yield : real,

  --/ whether this pft is the main crop in the stand (pft.id==stand.pftid)
  isprimarycrop : bool,
  --/ whether this pft is allowed to compete with the main crop during the same growing period (for future use)
  isprimarycovegetation : bool,
  --/ whether this pft is grown during a second growing period, different from the primary (main) crop (for future use)
--  bool issecondarycrop : real,

  --/ set to true if pft.isintercropgrass is true and the stand's main crop pft.intercrop is "naturalgrass"
  isintercropgrass : bool
}

--- Object describing sub-daily periods
type Day = {
  --- Whether sub-daily period first/last within day (both true in daily mode)
  isstart : bool,
  isend : bool,
  --- Ordinal number of the sub-daily period [0, date.subdaily)
  period : int
}

--- A vegetation individual.
--- In population mode this is the average individual of a PFT population,
--  in cohort mode: the average individual of a cohort,
--  in individual mode: an individual plant. Each grass PFT is represented as a single
--  individual in all modes. Individual objects are collected within list arrays of
--  class Vegetation (defined below), of which there is one for each patch, and include
--  a reference to their 'parent' Vegetation object. Use the createobj member function
--  of class Vegetation to add new individuals.
--/
type Individual = {
  --- reference to Pft object containing static parameters for this individual
  pft_id : int,
  --- reference to Vegetation object to which this Individual belongs - unnecessary, since Vegetation and its parent - Patch - are one-to-one.
--Vegetation& vegetation,
  gridcell_id : int,
  stand_id :  int,
  patch_id :  int,
  --- id code (0-based, sequential)
  individual_id : int,
  --- leaf C biomass on modelled area basis (kgC/m2)
  cmass_leaf : real,
  --- fine root C biomass on modelled area basis (kgC/m2)
  cmass_root : real,
  --- sapwood C biomass on modelled area basis (kgC/m2)
  cmass_sap : real,
  --- heartwood C biomass on modelled area basis (kgC/m2)
  cmass_heart : real,
  --- C "debt" (retrospective storage) (kgC/m2)
  cmass_debt : real,
  --- Total C mass at land use change (kgC/m2)
  cmass_tot_luc : real,
  --- leaf C mass after tunrnover
  cmass_leaf_post_turnover : real,
  --- root C mass after turnover
  cmass_root_post_turnover : real,
  --- Latest tunover day for this individual
  last_turnover_day : int,

  --- leaf N biomass on modelled area basis (kgN/m2)
  nmass_leaf : real,
  --- root N biomass on modelled area basis (kgN/m2)
  nmass_root : real,
  --- sap N biomass on modelled area basis (kgN/m2)
  nmass_sap : real,
  --- heart N biomass on modelled area basis (kgN/m2)
  nmass_heart : real,

  --- leaf N biomass on modelled area basis saved on first day of land use change year
  nmass_leaf_luc : real,
  --- root N biomass on modelled area basis on first day of land use change year
  nmass_root_luc : real,
  --- sap N biomass on modelled area basis on first day of land use change year
  nmass_sap_luc : real,
  --- heart N biomass on modelled area basis on first day of land use change year
  nmass_heart_luc : real,
  --- total N biomass on modelled area basis on first day of land use change year
  nmass_tot_luc : real,

  --- foliar projective cover (FPC) under full leaf cover as fraction of modelled area
  fpc : real,
  --- foliar projective cover (FPC) this day as fraction of modelled area
  fpc_daily : real,
  --- fraction of PAR absorbed by foliage over projective area today, taking account of leaf phenological state
  fpar : real,
  --- average density of individuals over patch (indiv/m2)
  densindiv : real,
  --- vegetation phenological state (fraction of potential leaf cover)
  phen : real,
  --- annual sum of daily fractional leaf cover
  --- Equivalent number of days with full leaf cover
  --  (population mode only : real, reset on expected coldest day of year)
  --/
  aphen : real,
  --- annual number of days with full leaf cover) (raingreen PFTs only : real, reset on 1 January)
  aphen_raingreen : int,

  --- Photosynthesis values for this individual under non-water-stress conditions
  photosynthesis_result : PhotosynthesisResult,

  --- sub-daily version of the above variable (NB: daily units)
  phots : [Date_subdaily]PhotosynthesisResult, -- TODO FIXME confirm the size parameter

  --- accumulated NPP over modelled area (kgC/m2/year),
  --- annual NPP following call to growth module on last day of simulation year */
  anpp : real,
  --- actual evapotranspiration over projected area (mm/day)
  aet : real,
  --- annual actual evapotranspiration over projected area (mm/year)
  aaet : real,
  --- leaf to root mass ratio
  ltor : real,
  --- plant height (m)
  height : real,
  --- plant crown area (m2)
  crownarea : real,
  --- increment in fpc since last simulation year
  deltafpc : real,
  --- bole height, i.e. height above ground of bottom of crown cylinder (m)
  --- (individual and cohort modes only) */
  boleht : real,
  --- patch-level lai for this individual or cohort (function fpar)
  lai : real,
  --- patch-level lai for cohort in current vertical layer (function fpar)
  lai_layer : real,
  --- individual leaf area index (individual and cohort modes only)
  lai_indiv : real,
  --- patch-level individual leaf area index (individual and cohort modes only)
  lai_daily : real,
  --- daily individual leaf area index (individual and cohort modes only)
  lai_indiv_daily : real,
  --- growth efficiency (NPP/leaf area) for each of the last five simulation years (kgC/m2/yr)
  --Historic<double, NYEARGREFF> greff_5,
  --- increment of wood C for each of the last five simulation years (kgC/m2/yr)
  --Historic<double, 10> cmass_wood_inc_5,
  --- individual/cohort age (years)
  age : real,
  --- monthly LAI (including phenology component)
  mlai : [12]real,
  --- monthly maximum LAI (including phenology component)
  mlai_max : [12]real,
  --- FPAR assuming full leaf cover for all vegetation
  fpar_leafon : real,
  --- LAI for current layer in canopy (cohort/individual mode : real, see function fpar)
  lai_leafon_layer : real,
  --- non-water-stressed canopy conductance on FPC basis (mm/s)
  gpterm : real,
  --- sub-daily version of the above variable (mm/s)
  gpterms : [Date_subdaily]real,
  --- interception associated with this individual today (patch basis)
  intercep : real,

  --- accumulated mean fraction of potential leaf cover
  phen_mean : real,

  --- whether individual subject to water stress
  wstress : bool,

  --- leaf nitrogen that is photosyntetic active
  nactive : real,
  --- Nitrogen extinction scalar
  --- Scalar to account for leaf nitrogen not following the optimal light
  --- extinction, but is shallower.
  nextin : real,
  --- long-term storage of labile nitrogen
  nstore_longterm : real,
  --- storage of labile nitrogen
  nstore_labile : real,
  --- long-term storage of labile nitrogen saved on first day of land use change year
  nstore_longterm_luc : real,
  --- storage of labile nitrogen saved on first day of land use change year
  nstore_labile_luc : real,
  --- daily total nitrogen demand
  ndemand : real,
  --- fraction of individual nitrogen demand available for uptake
  fnuptake : real,
  --- annual nitrogen uptake
  anuptake : real,
  --- maximum size of nitrogen storage
  max_n_storage : real,
  --- scales annual npp to maximum nitrogen storage
  scale_n_storage : real,
  --- annual nitrogen limitation on vmax
  avmaxnlim : real,
  --- annual optimal leaf C:N ratio
  cton_leaf_aopt : real,
  --- annual average leaf C:N ratio
  cton_leaf_aavr : real,
  --- plant mobile nitrogen status
  cton_status : real,
  --- total carbon in compartments before growth
  cmass_veg : real,
  --- total nitrogen in compartments before growth
  nmass_veg : real,
  --- whether individual subject to nitrogen stress
  nstress : bool,
  --- daily leaf nitrogen demand calculated from Vmax (kgN/m2)
  leafndemand : real,
  --- daily root nitrogen demand based on leafndemand
  rootndemand : real,
  --- daily sap wood nitrogen demand based on leafndemand
  sapndemand : real,
  --- daily labile nitrogen demand based on npp
  storendemand : real,
  --- daily harvestable organ nitrogen demand
  hondemand : real,
  --- leaf fraction of total nitrogen demand
  leaffndemand : real,
  --- root fraction of total nitrogen demand
  rootfndemand : real,
  --- sap fraction of total nitrogen demand
  sapfndemand : real,
  --- store fraction of total nitrogen demand
  storefndemand : real,
  --- daily leaf nitrogen demand over possible uptake (storage demand)
  leafndemand_store : real,
  --- daily root nitrogen demand over possible uptake (storage demand)
  rootndemand_store : real,

  --- The daily C lossed from leaves due to senescense, only crops.
  daily_cmass_leafloss : real,
  --- The daily N lossed from leaves due to senescense, only crops.
  daily_nmass_leafloss : real,
  --- The daily C lossed from roots due to senescense, only crops.
  daily_cmass_rootloss : real,
  --- The daily N lossed from roots due to senescense, only crops.
  daily_nmass_rootloss : real,

  --- Number of days with non-negligible phenology this month
  nday_leafon : int,
  -- Whether this individual is truly alive.
  --- Set to false for first year after the Individual object is created, then true. */
  alive : bool,
  --- NPP this day
  dnpp : real,

  -- bvoc

  --- isoprene production (mg C m-2 d-1)
  iso : real,
  --- monoterpene production (mg C m-2 d-1) per monoteprene species
  mon : [NMTCOMPOUNDS]real,
  --- monoterpene storage pool (mg C m-2) per monoterpene species
  monstor : [NMTCOMPOUNDS]real,
  --- isoprene seasonality factor (-)
  fvocseas : real,

  --- Pointer to struct with crop-specific data
  cropindiv : cropindiv_struct

}

-- type Vegetation -- We dont use this as it is just an array





--- Container for crop-specific data at patchpft level
type cropphen_struct = {
  --- latest sowing date
  sdate : int,
  --- sowing date of growing period ending in latest harvest this year
  sdate_harv : int,
  --- sowing dates of growing periods ending in the two latest harvests this year
  sdate_harvest : [2]int,
  --- sowing dates of growing periods starting this year
  sdate_thisyear : [2]int,
  --- number of sowings this year
  nsow : int,
  --- latest harvest date
  hdate : int,
  --- two latest harvest dates this year
  hdate_harvest : [2]int,
  --- last date for harvest
  hlimitdate : int,
  --- last day of heat unit sampling period, set in Crop_sowing_date_new()
  hucountend : int,
  --- number of harvests this year
  nharv : int,
  --- whether sdate_harvest[0] happened last year
  sownlastyear : bool,
  --- latest senescence start date this year
  sendate : int,
  --- latest beginning of intercropseason (2 weeks after the harvest date)
  bicdate : int,
  --- latest end of intercropseason (2 weeks before the sowing date)
  eicdate : int,
  --- number of growing days this growing period
  growingdays : int,
  --- number of growing days this year (used for wscal_mean calculation)
  growingdays_y : int,
  --- length of growingseason ending in last harvest
  lgp : int,
  --- base temp for heat unit calculation (°C)
  tb : real,
  --- number of vernalising days required
  pvd : int,
  --- number of accumulated vernalizing days
  vdsum : int,
  --- heat unit reduction factor due to vernalization [0-1]
  vrf : real,
  --- heat unit reduction factor due to photoperiodism [0-1]
  prf : real,
  --- potential heat units required for crop maturity (°Cd)
  phu : real,
  --- potential heat units that would have been used without dynamic phu calculation
  phu_old : real,
  --- heat unit sum aquired during last growing period (°Cd)
  husum : real,
  --- heat unit sum aquired durin sampling period, starting with sdate
  husum_sampled : real,
  --- this year's heat unit sum aquired from sdate to hucountend
  husum_max : real,
  --- running mean of recent past's husum_max
  husum_max_10 : real,
  --- number of heat units sampling years
  nyears_hu_sample : int,
  --- fraction of growing season [0-1] (husum/phu)
  fphu : real,
  --- fraction of growing season at latest harvest
  fphu_harv : real,
  --- whether in period of heat unit sampling
  hu_samplingperiod : bool,
  --- number of heat unit sampling days
  hu_samplingdays : int,
  --- harvest index today [0-1, >1 if below-ground ho], harvestable organ/above-ground C for above-ground harvestable organs, dependent on fphu, reduced by water stress
  hi : real,
  --- fraction of harvest index today
  fhi : real,
  --- phenology (fphu) contribution of fraction of harvest index today
  fhi_phen : real,  --Phenology (fPHU) compoment of fhi
  --- water stress contribution of fraction of harvest index today
  fhi_water : real,
  --- fraction of harvest index at latest harvest
  fhi_harv : real,
  --- sum of crop patch demand (patch.wdemand) during crop growing period, reset on harvest day
  demandsum_crop : real,
  --- sum of crop supply (patchpft.wsupply) during crop growing period, reset on harvest day
  supplysum_crop : real,

  --- whether inside crop/intercrop grass growing period
  growingseason : bool,
  --- whether yesterday was inside crop/intercrop grass growing period
  growingseason_ystd : bool,
  --- whether inside crop senescence
  senescence : bool,
  --- whether yesterday was inside crop senescence
  senescence_ystd : bool,
  --- whether inside intercrop crass growing period (main crop pft variable)
  intercropseason : bool,

  vdsum_alloc : real,
  vd : real,

  --- The fraction of the daily assimilates allocated to roots.
  f_alloc_root : real,
  --- The fraction of the daily assimilates allocated to leaves.
  f_alloc_leaf : real,
  --- The fraction of the daily assimilates allocated to harvestable organs, seeds.
  f_alloc_horg : real,
  --- The fraction of the daily assimilates allocated to stem.
  f_alloc_stem : real,
  --- Development stage from Wang & Engel 1998
  dev_stage : real,
  -- A variable holding the memory of whether this field was fertilised or not.
  fertilised  : [3]bool
}



--/ State variables common to all individuals of a particular PFT in a particular patch
-- Used in individual and cohort modes only.--
type Patchpft = {
  -- MEMBER VARIABLES:
  --/ id code (equal to value of member variable id in corresponding Pft object)
  gridcell_id : int,
  stand_id : int,
  patch_id : int,
  patchpft_id : int,
  --/ reference to corresponding Pft object in PFT list
  pft_id : int,
  --/ potential annual net assimilation (leaf-level net photosynthesis) at forest floor (kgC/m2/year)
  anetps_ff : real,
  --/ water stress parameter (0-1 range, 1=minimum stress)
  wscal : real,
  --/ running sum (converted to annual mean) for wscal
  wscal_mean : real,
  --/ potential annual net assimilation at forest floor averaged over establishment interval (kgC/m2/year)
  anetps_ff_est : real,
  --/ first-year value of anetps_ff_est
  anetps_ff_est_initial : real,
  --/ annual mean wscal averaged over establishment interval
  wscal_mean_est : real,
  --/ vegetation phenological state (fraction of potential leaf cover), updated daily
  phen : real,
  --/ annual sum of daily fractional leaf cover
  -- equivalent number of days with full leaf cover
  --  (reset on expected coldest day of year)
  aphen : real,
  --/ whether PFT can establish in this patch under current conditions
  establish : bool,
  --/ running total for number of saplings of this PFT to establish (cohort mode)
  nsapling : real,
  --/ leaf-derived litter for PFT on modelled area basis (kgC/m2)
  litter_leaf : real,
  --/ fine root-derived litter for PFT on modelled area basis (kgC/m2)
  litter_root : real,
  --/ remaining sapwood-derived litter for PFT on modelled area basis (kgC/m2)
  litter_sap : real,
  --/ remaining heartwood-derived litter for PFT on modelled area basis (kgC/m2)
  litter_heart : real,
  --/ litter derived from allocation to reproduction for PFT on modelled area basis (kgC/m2)
  litter_repr : real,

  --/ leaf-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
  nmass_litter_leaf : real,
  --/ root-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
  nmass_litter_root : real,
  --/ remaining sapwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
  nmass_litter_sap : real,
  --/ remaining heartwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
  nmass_litter_heart : real,

  --/ non-FPC-weighted canopy conductance value for PFT under water-stress conditions (mm/s)
  gcbase : real,
  --/ daily value of the above variable (mm/s)
  gcbase_day : real,

  --/ evapotranspirational "supply" function for this PFT today (mm/day)
  wsupply : real,
  wsupply_leafon : real,
  --/ fractional uptake of water from each soil layer today
  fwuptake : [NSOILLAYER]real,

  --/ whether water-stress conditions for this PFT
  wstress : bool,
  --/ daily version of the above variable
  wstress_day : bool,

  --/ carbon depository for long-lived products like wood
  harvested_products_slow : real,
  --/ nitrogen depository for long-lived products like wood
  harvested_products_slow_nmass : real,
  --/ first and last day of crop sowing window, calculated in crop_sowing_patch() or Crop_sowing_date_new()
  swindow : [2]int,
  --/ daily value of water deficit, calculated in irrigated_water_uptake()
  water_deficit_d : real,
  --/ yearly sum of water deficit
  water_deficit_y : real,

  --/ Struct for crop-specific variables
  cropphen : cropphen_struct,

  -- INUNDATION STRESS TERMS:

  --/ Number of days a month that the water table is above this PFT's wtp_max (updated on the first day of the month)
  inund_count : int,
  --/ [0,1] - a measure of the inundation stress. Daily photosynthesis is reduced by this factor.
  inund_stress : real
}


-- Stores data for a patch.
-- In cohort and individual modes, replicate patches are
--  required in each stand to accomodate stochastic variation, in population mode there
--  should be just one Patch object, representing average conditions for the entire
--  stand. A reference to the parent Stand object (defined below) is included as a
--  member variable.

type Patch = {
  -- id code in range 0-npatch for patch
  gridcell_id : int,
  stand_id : int,
  patch_id : int,

  -- use patch_id to access the array of Patchpft in the global state record
  --pfts: [npft]Patchpft,
  -- use patch_id to access the array of Individual in the global state record
  --vegetation: [npft]Individual,
  -- soil for this patch
  --soil_id: Soil,
  -- fluxes for this patch
  fluxes: Fluxes,
  -- FPAR at top of grass canopy today
  fpar_grass: real,
  -- FPAR at soil surface today
  fpar_ff: real,
  -- mean growing season PAR at top of grass canopy (J/m2/day)
  par_grass_mean: real,
  -- number of days in growing season, estimated from mean vegetation leaf-on fraction
  -- \see function fpar in canopy exchange module--
  nday_growingseason: int,
  -- total patch FPC
  fpc_total: real,
  -- whether patch was disturbed last year
  disturbed: bool,
  -- patch age (years since last disturbance)
  age: int,
  -- probability of fire this year (GlobFIRM)
  fireprob: real,

  -- BLAZE Fire line intensity (kW/m)
  fire_line_intensity: real,

  -- BLAZE fire related carbon fluxes
  -- BLAZE-fire carbon flux: live wood to atmosphere (kgC/m2)
  wood_to_atm: real,
  -- BLAZE-fire carbon flux: leaves to atmosphere (kgC/m2)
  leaf_to_atm: real,
  -- BLAZE-fire carbon flux: leaves to litter (kgC/m2)
  leaf_to_lit: real,
  -- BLAZE-fire carbon flux: live wood to structural litter (kgC/m2)
  wood_to_str: real,
  -- BLAZE-fire carbon flux: live wood to fine woody debris (kgC/m2)
  wood_to_fwd: real,
  -- BLAZE-fire carbon flux: live wood to coarse woody debris (kgC/m2)
  wood_to_cwd: real,
  -- BLAZE-fire carbon flux: fine litter to atmosphere (kgC/m2)
  litf_to_atm: real,
  -- BLAZE-fire carbon flux: fine woody debris to atmosphere (kgC/m2)
  lfwd_to_atm: real,
  -- BLAZE-fire carbon flux: coarse woody debris to atmosphere (kgC/m2)
  lcwd_to_atm: real,

  -- Storage for averaging of different Fapars for biome mapping in SIMFIRE
  -- SIMFIRE fapar: Grasses
  avg_fgrass: [N_YEAR_BIOMEAVG]real,
  -- SIMFIRE fapar: Needle-leaf tree
  avg_fndlt: [N_YEAR_BIOMEAVG]real,
  -- SIMFIRE fapar: Broad-leaf tree
  avg_fbrlt: [N_YEAR_BIOMEAVG]real,
  -- SIMFIRE fapar: Shrubs
  avg_fshrb: [N_YEAR_BIOMEAVG]real,
  -- SIMFIRE fapar: Total fapar
  avg_ftot: [N_YEAR_BIOMEAVG]real,

  -- whether management has started on this patch
  managed: bool,
  -- cutting intensity (initial percent of trees cut, further selection at individual level has to be done in a separate function)
  man_strength: real,

  managed_this_year: bool,
  plant_this_year: bool,

  -- DLE - the number of days over which wcont is averaged for this patch
  -- i.e. those days for which daily temp > 5.0 degC--
  growingseasondays: int,


  -- Variables used by new hydrology (Dieter Gerten 2002-07)

  -- interception by vegetation today on patch basis (mm)
  intercep: real,
  -- annual sum of AET (mm/year)
  aaet: real,
  -- annual sum of AET (mm/year) for each of the last five simulation years
  --Historic<double, NYEARAAET> aaet_5, TODO FIXME
  -- annual sum of soil evaporation (mm/year)
  aevap: real,
  -- annual sum of interception (mm/year)
  aintercep: real,
  -- annual sum of runoff (mm/year)
  asurfrunoff: real,
  -- annual sum of runoff (mm/year)
  adrainrunoff: real,
  -- annual sum of runoff (mm/year)
  abaserunoff: real,
  -- annual sum of runoff (mm/year)
  arunoff: real,
  -- water added to wetlands today (mm)
  wetland_water_added_today: real,
  -- annual sum of water added to wetlands (mm/year)
  awetland_water_added: real,
  -- annual sum of potential evapotranspiration (mm/year)
  apet: real,

  -- equilibrium evapotranspiration today, deducting interception (mm)
  eet_net_veg: real,

  -- transpirative demand for patch, patch vegetative area basis (mm/day)
  wdemand: real,
  -- daily average of the above variable (mm/day)
  wdemand_day: real,
  -- transpirative demand for patch assuming full leaf cover today
  -- mm/day, patch vegetative area basis --
  wdemand_leafon: real,
  -- rescaling factor to account for spatial overlap between individuals/cohorts populations
  fpc_rescale: real,

  -- monthly AET (mm/month)
  maet: [12]real,
  -- monthly soil evaporation (mm/month)
  mevap: [12]real,
  -- monthly interception (mm/month)
  mintercep: [12]real,
  -- monthly runoff (mm/month)
  mrunoff: [12]real,
  -- monthly PET (mm/month)
  mpet: [12]real,

  -- daily nitrogen demand
  ndemand: real,

  -- annual nitrogen fertilization (kgN/m2/year)
  anfert: real,
  -- daily nitrogen fertilization (kgN/m2/day)
  dnfert: real,

  -- daily value of irrigation water (mm), set in irrigation(), derived from water_deficit_d
  irrigation_d: real,
  -- yearly sum of irrigation water (mm)
  irrigation_y: real,

  -- whether litter is to be sent to the soil today
  is_litter_day: bool,
  -- number of harvests and/or cover-crop killing or turnover events
  nharv: int,
  -- whether today is a harvest day and/or cover-crop killing or turnover day
  isharvestday: bool
}





--/ Container for variables common to individuals of a particular PFT in a stand.
-- Used in individual and cohort modes only
--
type Standpft = {

  -- MEMBER VARIABLES
  gridcell_id : int,
  standpft_id : int,
  stand_id : int,
  pft_id : int,
-- net C allocated to reproduction for this PFT in all patches of this stand this year (kgC/m2)
  cmass_repr : real,
-- maximum value of Patchpft::anetps_ff for this PFT in this stand so far in the simulation (kgC/m2/year)
  anetps_ff_max : real,

-- FPC sum for this PFT as average for stand
  fpc_total : real,

-- Photosynthesis values for this PFT under non-water-stress conditions
  photosynthesis_result : PhotosynthesisResult,

-- Whether this PFT is allowed to grow in this stand
  active : bool,
-- Whether this PFT is planted in this stand
  plant : bool,
-- Whether this PFT is allowed to establish (after planting) in this stand
  reestab : bool,

-- Whether this PFT is irrigated in this stand
  irrigated : bool,
-- sowing date specified in stand type or read from input file
  sdate_force : int,
-- harvest date specified in stand type or read from input file
  hdate_force : int
  -- MEMBER FUNCTIONS

-- Constructor: initialises various data members
}

--/ The stand class corresponds to a modelled area of a specific landcover type in a grid cell.
-- There may be several stands of the same landcover type (but with different settings).
--
type Stand = {
  --data : [num_patches]Patch,
  -- MEMBER VARIABLES
-- list array [0...npft-1] of Standpft (initialised in constructor)
  --standpft : [npft]Standpft, -- this was called 'pft' in the c++ code, which is confusing

-- A number identifying this Stand within the grid cell
  original : landcovertype,
  stand_id : int,
  gridcell_id : int,

  num_patches : int, -- since Patches is aggressively _______, in order to ensure regularity, use this field to avoid using those extras

-- stand type id
  standtype_id : int,

-- pft id of main crop, updated during rotation
  pft_id : int,

-- current crop rotation item
  current_rot : int,
-- number of days passed in current rotation item
  ndays_inrotation : int,
-- Returns true if stand is in fallow (with cover crop grass)
  infallow : bool,
-- Returns true if crop rotation item is to be updated today
  isrotationday : bool,
-- Returns true if current crop management hydrology == irrigated, updated during rotation
  isirrigated : bool,
-- Returns true if the stand's main crop pft intercrop==naturalgrass and a pft with isintercrop==true is in the pftlist.
  hasgrassintercrop : bool,
-- gdd5-value at first intercrop grass growth
  gdd5_intercrop : real,

-- old fraction of this stand relative to the gridcell before update
  frac_old : real,

-- used during land cover change involving several calls to reveiving_stand_change()
-- Set to frac_old in reduce_stands(), then modified in donor_stand_change() and receiving_stand_change().
  --
  frac_temp : real,
-- fraction unavailable for transfer to other stand types
  protected_frac : real,
-- net stand fraction change
  frac_change : real,
-- gross fraction increase
  gross_frac_increase : real,
-- gross fraction decrease
  gross_frac_decrease : real,
-- fraction that has been cloned from another stand
  cloned_fraction : real,
-- Returns true if this stand is cloned from another stand
  cloned : bool,
-- pointer to array of fractions transferred from this stand to other stand types
  transfer_area_st : [nst]real,
-- land cover origin of this stand
  origin : landcovertype,
-- used for output from separate stands
  anpp : real,
-- used for output from separate stands
  cmass : real,

-- Seed for generating random numbers within this Stand
-- The reason why Stand has its own seed, rather than using for instance
  -- a single global seed is to make it easier to compare results when using
  -- different land cover types.
  --
  --  Randomness not associated with a specific stand, but rather a whole
  --  grid cell should instead use the seed in the Gridcell class.
  --
  --  \see randfrac()
  --
  seed : int,

-- type of landcover
-- \see landcovertype
-- initialised in constructor
  --
  landcover : landcovertype,

-- The year when this stand was created.
-- Will typically be year zero unless running with dynamic
  --  land cover.
  --
  --  Needed to set patchpft.anetps_ff_est_initial
  --
  first_year : int,
  -- The year this stand was cloned from another stand
  clone_year : int,
-- scaling factor for stands that have grown in area this year (old fraction/new fraction)
  scale_LC_change : real,

  -- Pointer to parent object, could be a null pointer
  -- Prefer to access the gridcell through get_gridcell(), even internally
  --  within the Stand class.
  --gridcell : Gridcell,
  --gridcellIdx : int,

  -- Soil type to be used in this Stand
  soiltype_id : int,

  -- Fraction of this stand relative to the gridcell
  -- used by crop stands, initialized in constructor to 1,
  --  set in landcover_init()
  frac : real

}



--/ Storage of land cover fraction data and some land cover change-related pools and fluxes
type Landcover = {
  --/ The fractions of the different land cover types.
  --  landcoverfrac is read in from land cover input file or from
  --  instruction file in getlandcover().
  --/
  frac : [NLANDCOVERTYPES]real,

  --/ The land cover fractions from the previous year
  --  Used to keep track of the changes when running with dynamic
  --  land cover.
  --/
  frac_old : [NLANDCOVERTYPES]real,

  frac_change : [NLANDCOVERTYPES]real,

  --/ Transfer matrices
  frac_transfer : [NLANDCOVERTYPES][NLANDCOVERTYPES]real,
  primary_frac_transfer : [NLANDCOVERTYPES][NLANDCOVERTYPES]real,

  --/ Whether the land cover fractions changed for this grid cell this year
  --  \see landcover_dynamics
  --/
  updated : bool,

  --/ Gridcell-level C flux from slow harvested products
  acflux_harvest_slow : real,

  --/ Gridcell-level C flux from harvest associated with landcover change
  acflux_landuse_change : real,

  --/ Gridcell-level N flux from slow harvested products
  anflux_harvest_slow : real,

  --/ Gridcell-level N flux from harvest associated with landcover change
  anflux_landuse_change : real,

  --/ Landcover-level C flux from slow harvested products (donating landcover)
  acflux_harvest_slow_lc : [NLANDCOVERTYPES]real,

  --/ Landcover-level C flux from harvest associated with landcover change (donating landcover)
  acflux_landuse_change_lc : [NLANDCOVERTYPES]real,

  --/ Landcover-level N flux from slow harvested products (donating landcover)
  anflux_harvest_slow_lc : [NLANDCOVERTYPES]real,

  --/ Landcover-level N flux from harvest associated with landcover change (donating landcover)
  anflux_landuse_change_lc : [NLANDCOVERTYPES]real,

  --/ Which landcover types create new stands when area increases.
  expand_to_new_stand : [NLANDCOVERTYPES]bool,

  --/ Whether to pool all transferred land from a donor landcover (overrides different landcover targets of different stand types and stands in a landcover)
  pool_to_all_landcovers : [NLANDCOVERTYPES]bool,

  --/ Whether to pool transferred land to a receptor landcover (crop and pasture stands to new natural stand: pool!)
  pool_from_all_landcovers : [NLANDCOVERTYPES]bool
}



--/ State variables common to all individuals of a particular PFT in a GRIDCELL.
type Gridcellpft = {
  -- MEMBER VARIABLES
  --/ A number identifying this object within its list array
  gridcellpft_id : int,

  --/ A reference to the Pft object for this Gridcellpft
  pft_id : int,

  --/ annual degree day sum above threshold damaging temperature
  --  used in calculation of heat stess mortality, Sitch et al 2000, Eqn 55
  --/
  addtw : real,

  --/ Michaelis-Menten kinetic parameters
  --  Half saturation concentration for N uptake (Rothstein 2000, Macduff 2002)
  --/
  Km : real,

  --/Crop-specific variables:
  --/ whether the daily temperature has fallen below the autumn temperature limit (tempautumn) this year
  autumnoccurred : bool,
  --/ whether the daily temperature has risen above the spring temperature limit (tempspring) this year
  springoccurred : bool,
  --/ whether the daily temperature has fallen below the vernalization limit (trg) this year
  vernstartoccurred : bool,
  --/ whether the daily temperature rises over the vernalization limit (trg) this year
  vernendoccurred : bool,
  --/ first day when temperature fell below the autumn temperature limit (tempautumn) this year
  first_autumndate : int,
  --/ 20-year mean
  first_autumndate20 : int,
  --/ memory of the last 20 years' values
  first_autumndate_20 : [20]int,
  --/ last day when temperature rose above the spring temperature limit (tempspring) this year
  last_springdate : int,
  --/ 20-year mean
  last_springdate20 : int,
  --/ memory of the last 20 years' values
  last_springdate_20 : [20]int,
  --/ last day when temperature has fallen below the vernilisation temperature limit (trg) this year (if vernstartoccurred==true)
  last_verndate : int,
  --/ 20-year mean
  last_verndate20 : int,
  --/ memory of the last 20 years' values
  last_verndate_20 : [20]int,
  --/ default sowing date (pft.sdatenh/sdatesh)
  sdate_default : int,
  --/ calculated sowing date from temperature limits
  sdatecalc_temp : int,
  --/ calculated sowing date from precipitation limits
  sdatecalc_prec : int,
  --/ sowing date from input file
  sdate_force : int,
  --/ harvest date from input file
  hdate_force : int,
  --/ N fertilization from input file
  Nfert_read : real,
  --/ Manure N fertilization from input file
  Nfert_man_read : real,
  --/ default harvest date (pft.hlimitdatenh/hlimitdatesh)
  hlimitdate_default : int,
  --/ whether autumn sowing is either calculated or prescribed
  wintertype : bool,
  --/ first and last day of crop sowing window, calculated in calc_sowing_windows()
  swindow : [2]int,
  --/ first and last day of crop sowing window for irrigated crops, calculated in calc_sowing_windows()
  swindow_irr : [2]int,
  --/ temperature limits precludes crop sowing
  sowing_restriction : bool
}

--/ State variables common to all individuals of a particular STANDTYPE in a GRIDCELL.
type Gridcellst = {
  -- MEMBER VARIABLES
  --/ A number identifying this object within its list array
  gridcellst_id : int,

  --/ A reference to the StandType object for this Gridcellst
  standtype_id : int,

  --/ fraction of this stand type relative to the gridcell
  frac : real,
  --/ old fraction of this stand type relative to the gridcell before update
  frac_old : real,
  --/ original input value of old fraction of this stand type before rescaling
  frac_old_orig : real,
  --/ fraction unavailable for transfer to other stand types
  protected_frac : real,

  --/ net fraction change
  frac_change : real,
  --/ gross fraction increase
  gross_frac_increase : real,
  --/ gross fraction decrease
  gross_frac_decrease : real,

  -- current number of stands of this stand type
  nstands : int,

  nfert : real
}

--/ The Gridcell class corresponds to a modelled locality or grid cell.
--- Member variables include an object of type Climate (holding climate, insolation and
--  CO2 data), a object of type Soiltype (holding soil static parameters) and a list
--  array of Stand objects. Soil objects (holding soil state variables) are associated
--  with patches, not gridcells. A separate Gridcell object must be declared for each modelled
--  locality or grid cell.
--/
type Gridcell = {

  --inherited data: should be accessed in the global Stand array, using this Gridcell_id

  -- MEMBER VARIABLES
  --/ climate, insolation and CO2 for this grid cell
  climate_id : int,

  --/ soil static parameters for this grid cell
  soiltype_id : int,

  --/ landcover fractions and landcover-specific variables
  landcover_id : int,

  --/ list array [0...npft-1] of Gridcellpft (initialised in constructor)
  gridcellpft_id : int,

  --/ list array [0...nst-1] of Gridcellst (initialised in constructor)
  gridcellst_id : int,

  --/ object for keeping track of carbon and nitrogen balance
  massbalance_id : int,

  -- SIMFIRE
  --/ the region index to chosose from set of optimisations
  simfire_region : int,
  --/ timeseries of population density from the Hyde 3.1 dataset (inhabitants/ha)
  hyde31_pop_density : [57]real,
  --/ current year's population density (inhabitants/ha)
  pop_density : real,
  --/ tuning factor for available litter
  k_tun_litter : real,
  --/ maximum annual Nesterov Index
  max_nesterov : real,
  --/ current Nexterov index
  cur_nesterov : real,
  --/ Monthly max Nexterov index (to keep track of running year)
  monthly_max_nesterov : [12]real,
  --/ biome classification used in SIMFIRE
  simfire_biome : int,
  --/ Average maximum annual fAPAR (over avg_interv_fpar years)
  ann_max_fapar : real,
  --/ Average maximum annual fAPAR of recent years
  recent_max_fapar : [AVG_INTERVAL_FAPAR]real,
  --/ maximum fapar of running year so far
  cur_max_fapar : real,
  --/ monthly fire risk (factor describing local monthly fire climatology)
  monthly_fire_risk : [12]real,
  --/ current burned area from SIMFIRE (fract.)
  burned_area : real,
  --/ accumulated burned area from SIMFIRE for tstep < 1a (fract.)
  burned_area_accumulated : real,
  --/ Simple tracker to check whether at least one patch has enough fuel to burn
  can_burn : int,
  --/ annual burned area from SIMFIRE (fract.)
  annual_burned_area : real,
  --/ monthly burned area from SIMFIRE (fract.)
  monthly_burned_area : [12]real,

  -- Nitrogen deposition
  --/ annual NH4 deposition (kgN/m2/year)
  aNH4dep : real,
  --/ annual NO3 deposition (kgN/m2/year)
  aNO3dep : real,
  --/ daily NH4 deposition (kgN/m2)
  dNH4dep : real,
  --/ daily NO3 deposition (kgN/m2)
  dNO3dep : real,

  --/ Seed for generating random numbers within this Gridcell
  --  The reason why Gridcell has its own seed, rather than using for instance
  --  a single global seed is to make it easier to compare results when for
  --  instance changing the order in which the simulation proceeds. It also
  --  gets serialized together with the rest of the Gridcell state to make it
  --  possible to get exactly identical results after a restart.
  --
  --  \see randfrac()
  --/
  seed : int, -- long

  -- MEMBER FUNCTIONS

  --/ Creates a new Stand in this grid cell
  --Stand& create_stand(landcovertype lc, int no_patch = 0),

  --/ Creates new stand and initiates land cover settings when run_landcover==true
  --Stand& create_stand_lu(StandType& st, double fraction, int no_patch = 0),




  --/ Longitude for this grid cell
  lon : real,

  --/ Latitude for this grid cell
  lat : real

}
