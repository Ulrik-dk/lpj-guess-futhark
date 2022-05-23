-------------------------------------guess.h--------------------------------------
--------------------------------------------------------------------------------------/
-- #INCLUDES FOR LIBRARY HEADER FILES
-- C/C++ libraries required for member functions of classes defined in this file.
-- These libraries will also be available globally (so omit these #includes from source
-- files). In addition to various standard C/C++ runtime libraries, the framework
-- requires the following libraries (individual modules may use additional libraries)
--
-- GUTIL
--   Includes class xtring, providing functionality for pointer-free dynamic handling
--   of character strings wherever possible in LPJ-GUESS, strings are represented as
--   objects of type xtring rather than simple arrays of type char. GUTIL also provides
--   templates for dynamic collection classes (list arrays of various types), argument
--   processing for printf-style functions, timing functions and other utilities.

--#include "gutil.h"
--#include "shell.h"
open import "../framework/guessmath"
--#include "archive.h"
open import "../framework/parameters"
--#include "guesscontainer.h"
open import "../modules/soil"
open import "../futhark-extras"

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
-- Nitrogen preferance
type n_pref_type = enum_type
let NO : n_pref_type = 0
let NH4 : n_pref_type = 1
let NO3 : n_pref_type = 2
----------------------------------------------------------
-- GLOBAL CONSTANTS

-- number  of soil layers modelled
let NSOILLAYER_UPPER : int = 5
let NSOILLAYER_LOWER : int = NSOILLAYER - NSOILLAYER_UPPER

-- bvoc: number of monoterpene species used
let NMTCOMPOUNDS : int = NMTCOMPOUNDTYPES

-- SOIL DEPTH VALUES

-- soil upper layer depth (mm)
let SOILDEPTH_UPPER : real = 500.0
-- soil lower layer depth (mm)
let SOILDEPTH_LOWER : real = 1000.0

-- Depth of sublayer at top of upper soil layer, from which evaporation is
-- possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
-- Must be a multiple of Dz_soil
let SOILDEPTH_EVAP : real = 200.0

-- Year at which to calculate equilibrium soil carbon
let SOLVESOM_END : int = 400

-- Year at which to begin documenting means for calculation of equilibrium soil carbon
let SOLVESOM_BEGIN : int = 350

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


-- type Date
-- type Day
-- type MassBalance


-- This struct contains the result of a photosynthesis calculation.
-- see photosynthesis
type PhotosynthesisResult = {
  -- Constructs an empty result

  -- RuBisCO capacity (gC/m2/day)
  vm : real,

  -- gross daily photosynthesis (gC/m2/day)
  agd_g : real,

  -- leaf-level net daytime photosynthesis
  -- expressed in CO2 diffusion units (mm/m2/day) */
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

let net_assimilation(psr : PhotosynthesisResult) : real =
  (psr.agd_g - psr.rd_g) * KG_PER_G

--- Clears all members
--  This is returned by the photosynthesis function when no photosynthesis
--  takes place.
let PhotosynthesisResult() : PhotosynthesisResult =
  {
    agd_g       = 0,
    adtmm       = 0,
    rd_g        = 0,
    vm          = 0,
    je          = 0,
    nactive_opt = 0.0,
    vmaxnlim    = 1.0
  }


  -- type WeatherGenState


-- This struct contains the environmental input to a photosynthesis calculation.
-- \see photosynthesis */
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

-- Constructs an empty result
let PhotosynthesisEnvironment() : PhotosynthesisEnvironment =
  {
    co2 = 0,
    temp = 0,
    par = 0,
    fpar = 0,
    daylength = 0
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


let PhotosynthesisStresses() : PhotosynthesisStresses =
  {
    ifnlimvmax = false,
    moss_ps_limit = 1.0,
    graminoid_ps_limit = 1.0,
    inund_stress = 1.0
  }


-- type Climate


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
  -- For the fluxes only stored as totals for the whole patch */
  monthly_fluxes_patch : [12][NPERPATCHFLUXTYPES]real,

  -- Stores one flux value per month and flux type
  --- For the fluxes stored per pft for annual values */
  monthly_fluxes_pft : [12][NPERPFTFLUXTYPES]real,

  -- Stores one flux value per day and flux type
  daily_fluxes_patch : [365][NPERPATCHFLUXTYPES]real,

  -- Stores one flux value per day and flux type
  daily_fluxes_pft : [365][NPERPFTFLUXTYPES]real

}


--------------------------------------------------------------------------------
-- Implementation of Fluxes member functions
--------------------------------------------------------------------------------

-- TODO check if this constructor is correct
let Fluxes() : Fluxes = {
    --- Reference to patch to which this Fluxes object belongs
    --Patch& patch, --FUTHARK: ID instead of pointer? What is it used for?

    ---- Private members, mutable
    -- Stores one flux value per PFT and flux type
    annual_fluxes_per_pft = replicate npft (replicate NPERPFTFLUXTYPES realzero),

    -- Stores one flux value per month and flux type
    -- For the fluxes only stored as totals for the whole patch */
    monthly_fluxes_patch = replicate 12 (replicate NPERPATCHFLUXTYPES realzero),

    -- Stores one flux value per month and flux type
    --- For the fluxes stored per pft for annual values */
    monthly_fluxes_pft = replicate 12 (replicate NPERPFTFLUXTYPES realzero),

    -- Stores one flux value per day and flux type
    daily_fluxes_patch = replicate 365 (replicate NPERPATCHFLUXTYPES realzero),

    -- Stores one flux value per day and flux type
    daily_fluxes_pft = replicate 365 (replicate NPERPFTFLUXTYPES realzero)
  }

type Date = {
  month : int,
  day : int
}

-- Forgive me for this, there is no other way of updating arrays in records.
let report_flux_PerPFTFluxType(
    { annual_fluxes_per_pft
    , monthly_fluxes_patch
    , monthly_fluxes_pft
    , daily_fluxes_patch
    , daily_fluxes_pft }: *Fluxes
    , flux_type: PerPFTFluxType
    , pft_id: int
    , value: real
    , date: Date) : Fluxes =

    { annual_fluxes_per_pft =
        annual_fluxes_per_pft with [pft_id, flux_type] =
          annual_fluxes_per_pft[pft_id, flux_type] + value
    , monthly_fluxes_patch
    , monthly_fluxes_pft =
        monthly_fluxes_pft with [date.month, flux_type] =
          monthly_fluxes_pft[date.month, flux_type] + value
    , daily_fluxes_patch
    , daily_fluxes_pft = daily_fluxes_pft with [date.day,flux_type] = daily_fluxes_pft[date.day,flux_type] + value
    }

let report_flux_PerPatchFluxType(
    { annual_fluxes_per_pft
    , monthly_fluxes_patch
    , monthly_fluxes_pft
    , daily_fluxes_patch
    , daily_fluxes_pft } : *Fluxes
    , flux_type: PerPatchFluxType
    , value: real
    , date: Date) : Fluxes =

    { annual_fluxes_per_pft
    , monthly_fluxes_patch =
        monthly_fluxes_patch with [date.month, flux_type] =
          monthly_fluxes_patch[date.month, flux_type] + value
    , monthly_fluxes_pft
    , daily_fluxes_patch =
        daily_fluxes_patch with [date.day, flux_type] =
          daily_fluxes_patch[date.day, flux_type] + value
    , daily_fluxes_pft
    }


let get_daily_flux_PerPFTFluxType(this: Fluxes, flux_type: PerPFTFluxType, day: int) : real =
  this.daily_fluxes_pft[day,flux_type]

let get_daily_flux_PerPatchFluxType(this: Fluxes, flux_type: PerPatchFluxType, day: int) : real =
  this.daily_fluxes_patch[day,flux_type]


let get_monthly_flux_PerPFTFluxType(this: Fluxes, flux_type: PerPFTFluxType, month: int) : real =
  this.monthly_fluxes_pft[month,flux_type]

let get_monthly_flux_PerPatchFluxType(this: Fluxes, flux_type: PerPatchFluxType, month: int) : real =
  this.monthly_fluxes_patch[month,flux_type]

-- I dont trust this function -Ulrik
--let get_annual_flux_PerPFTFluxType_pft_id(this: Fluxes, flux_type: PerPFTFluxType, pft_id: int) : real =
--  this.annual_fluxes_per_pft[pft_id][flux_type]

-- its unfortunate that the indices are ordered as they are
let get_annual_flux_PerPFTFluxType(this: Fluxes, flux_type: PerPFTFluxType) : real =
  reduce (+) realzero <| (transpose this.annual_fluxes_per_pft)[flux_type]

let get_annual_flux_PerPatchFluxType(this: Fluxes, flux_type: PerPatchFluxType) : real =
  reduce (+) realzero <| (transpose this.annual_fluxes_per_pft)[flux_type]


-- type CropRotation
-- type StandType


-- Holds static functional parameters for a plant functional type (PFT).
-- There should be one Pft object for each potentially occurring PFT. The same Pft object
-- may be referenced (via the pft member of the Individual object see below) by different
-- average individuals. Member functions are included for initialising SLA given leaf
-- longevity, and for initialising sapling/regen characteristics (required for
-- population mode).

type Pft  = {
  -- MEMBER VARIABLES
  -- id code (should be zero based and sequential, 0...npft-1)
  id: int,
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
  -- Half saturation concentration for N uptake [kgN l-1] (Rothstein 2000) */
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
  -- J/m2/day, individual and cohort modes */
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
  -- \see ifbgestab */
  kest_bg: real,
  -- constant used in calculation of sapling establishment rate when spatial
  -- mass effect disabled (individual and cohort modes)

  kest_pres: real,
  -- expected longevity under non-stressed conditions (individual and cohort modes)
  longevity: real,
  -- threshold growth efficiency for imposition of growth suppression mortality
  -- kgC/m2 leaf/year, individual and cohort modes */
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

  ---- FIXME TODO HACK! comes from Date::MAX_YEAR_LENGTH+1
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
  -- Used to implement drought-limited establishment */
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
  -- kgC/m2 */

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

let avg_cton (min: real, max: real) : real =
  2.0 / (1.0 / min + 1.0 / max)

-- MEMBER FUNCTIONS
-- Constructor (initialises array gdd0)
let Pft() : Pft =
  {
    --std::fill_n(gdd0, Date::MAX_YEAR_LENGTH + 1, -1.0) -- value<0 signifies "unknown" see function phenology()
    gdd0 = replicate Date_MAX_YEAR_LENGTH_plusone (-1.0),
    nlim = false,

    root_beta = 0.0,

    drought_tolerance = 0.0, -- Default, means that the PFT will never be limited by drought.
    res_outtake = 0.0,
    harv_eff = 0.0,
    harv_eff_ic = 0.0,
    turnover_harv_prod = 1.0,  -- default 1 year turnover time

    isintercropgrass = false,
    ifsdautumn = false,
    maxtemp_sowing = 60,
    sdatenh = -1,
    sdatesh = -1,
    lgp_def = 190,
    hlimitdatenh = -1,
    hlimitdatesh = -1,
    tb = -999.9,
    trg = -999.9,
    pvd = -1,
    psens = -1.0,
    pb = -1.0,
    vern_lag=0,
    ps = -1.0,
    phu = -1.0,
    phu_red_spring_sow = 1.0,
    fphusen = -1.0,
    shapesenescencenorm = false,
    flaimaxharvest = -1.0,
    laimax = 0.0,
    aboveground_ho = true,
    frootstart = 0.0,
    frootend = 0.0,
    forceautumnsowing = 0,

    fertrate = (0.0,1.0),
    fertdates = (0,30),

    fert_stages = (0.5, 0.9),
    fertilised = (false, false),

    N_appfert = 0.0,

    T_vn_min=0.0,
    T_vn_opt=0.0,
    T_vn_max=0.0,
    T_veg_min=0.0,
    T_veg_opt=0.0,
    T_veg_max=0.0,
    T_rep_min=0.0,
    T_rep_opt=0.0,
    T_rep_max=0.0,

    a1=0.0,
    b1=0.0,
    c1=0.0,
    d1=0.0,
    a2=0.0,
    b2=0.0,
    c2=0.0,
    d2=0.0,
    a3=0.0,
    b3=0.0,
    c3=0.0,
    d3=0.0,

    photo = (0.0, 0.0, 0.0),

    ----- THE FOLLOWING WERE NOT DEFINED IN THE ORIGINAL CONSTRUCTOR

    alphar = nan,
    crownarea_max = nan,
    cton_leaf_avr = nan,
    cton_leaf_max = nan,

    cton_leaf_min = nan,
    cton_root = nan,
    cton_root_avr = nan,
    cton_root_max = nan,

    cton_sap = nan,
    cton_sap_avr = nan,
    cton_sap_max = nan,
    cton_stem_avr = nan,

    cton_stem_max = nan,
    dev_rate_rep = nan,
    dev_rate_veg = nan,
    emax = nan,
    eps_iso = nan,

    eps_mon = replicate NMTCOMPOUNDS nan,
    est_max = nan,
    fireresist = nan,
    fnstorage = nan,
    ga = nan,
    gdd0_max = nan,

    gdd0_min = nan,
    gdd5min_est = nan,
    gmin = nan,
    greff_min = nan,
    harvest_slow_frac = nan,

    has_aerenchyma = false,
    himin = nan,
    hiopt = nan,
    intc = nan,
    id = -1,
    inund_duration = -1,

    k_allom1 = nan,
    k_allom2 = nan,
    k_allom3 = nan,
    k_chilla = nan,
    k_chillb = nan,
    k_chillk = nan,

    k_latosa = nan,
    k_rp = nan,
    kest_bg = nan,
    kest_pres = nan,
    kest_repr = nan,
    km_volume = nan,

    lambda_max = nan,
    landcover = URBAN,
    leaflong = nan,
    leafphysiognomy = NOLEAFTYPE,
    lifeform = NOLIFEFORM,

    litterme = nan,
    longevity = nan,
    ltor_max = nan,
    max_snow = nan,
    min_snow = nan,
    name = -1,

    ndays_ramp_phu = nan,
    nupscoeff = nan,
    nuptoroot = nan,
    parff_min = nan,
    pathway = C4,

    phengdd5ramp = nan,
    phenology = CROPGREEN,
    phu_calc_lin = false,
    phu_calc_quad = false,

    phu_interc = nan,
    phu_max = nan,
    phu_min = nan,
    pstemp_high = nan,
    pstemp_low = nan,

    pstemp_max = nan,
    pstemp_min = nan,
    regen_cmass_heart = nan,
    regen_cmass_leaf = nan,

    regen_cmass_root = nan,
    regen_cmass_sap = nan,
    reprfrac = nan,
    respcoeff = nan,

    rootdist = replicate NSOILLAYER nan,
    sd_adjust = false,
    sd_adjust_par1 = nan,
    sd_adjust_par2 = nan,

    sd_adjust_par3 = nan,
    seas_iso = false,
    selection = -1,
    sla = nan,
    storfrac_mon = replicate NMTCOMPOUNDS nan,

    tcmax_est = nan,
    tcmin_est = nan,
    tcmin_surv = nan,
    tempautumn = nan,
    tempspring = nan,

    turnover_leaf = nan,
    turnover_root = nan,
    turnover_sap = nan,
    twmin_est = nan,

    twminusc = nan,
    wooddens = nan,
    wscal_min = nan,
    wtp_max = nan
  }

-- Calculates SLA given leaf longevity
let initsla(this: Pft) : Pft =
  -- SLA has to be supplied in the insfile for crops with N limitation
  if (!(this.phenology == CROPGREEN && this.nlim)) then
    -- Reich et al 1992, Table 1 (includes conversion x2.0 from m2/kg_dry_weight to
    -- m2/kgC)
    this with sla =
    if (this.leafphysiognomy == BROADLEAF) then
      0.2 * pow(10.0, 2.41 - 0.38 * log10(12.0 * this.leaflong))
    else if (this.leafphysiognomy == NEEDLELEAF) then
      0.2 * pow(10.0, 2.29 - 0.4 * log10(12.0 * this.leaflong))
    else this.sla
  else this


-- Calculates minimum leaf C:N ratio given leaf longevity
let init_cton_min(this: Pft) : Pft =
  -- cton_leaf_min has to be supplied in the insfile for crops with N limitation
  if (!(this.phenology == CROPGREEN && this.nlim)) then
    -- Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
    -- kgN/kgC)

    this with cton_leaf_min =
    if (this.leafphysiognomy == BROADLEAF) then
      500.0 / pow(10.0, 1.75 - 0.33 * log10(12.0 * this.leaflong))
    else if (this.leafphysiognomy == NEEDLELEAF) then
      500.0 / pow(10.0, 1.52 - 0.26 * log10(12.0 * this.leaflong))
    else this.cton_leaf_min
  else this

let init_cton_limits(this: Pft) : Pft =
  -- Fraction between min and max C:N ratio White et al. 2000
  let frac_mintomax = if (this.phenology == CROPGREEN && this.nlim) then 5.0 else 2.78  -- Use value also without nlim ?

  -- Fraction between leaf and root C:N ratio
  let frac_leaftoroot = 1.16 -- Friend et al. 1997

  -- Fraction between leaf and sap wood C:N ratio
  let frac_leaftosap = 6.9   -- Friend et al. 1997

  -- Max leaf C:N ratio
  let this = this with cton_leaf_max = this.cton_leaf_min * frac_mintomax

  -- Average leaf C:N ratio
  let this = this with cton_leaf_avr = avg_cton(this.cton_leaf_min, this.cton_leaf_max)

  -- Tighter C:N ratio range for roots and sapwood: picked out thin air
  let frac_maxtomin = 0.9

  -- Maximum fine root C:N ratio
  let this = this with cton_root_max = this.cton_leaf_max * frac_leaftoroot

  let cton_root_min = this.cton_root_max * frac_maxtomin

  -- Average fine root C:N ratio
  let this = this with cton_root_avr = avg_cton(cton_root_min, this.cton_root_max)

  -- Maximum sap C:N ratio
  let this = this with cton_sap_max  = this.cton_leaf_max * frac_leaftosap

  let cton_sap_min = this.cton_sap_max * frac_maxtomin

  -- Average sap C:N ratio
  let this = this with cton_sap_avr  = avg_cton(cton_sap_min, this.cton_sap_max)

  let this = this with respcoeff =
  if (this.lifeform == GRASS || this.lifeform == MOSS) then
    this.respcoeff / 2.0 * this.cton_root / (this.cton_root_avr + cton_root_min)
  else
    this.respcoeff / this.cton_root / (this.cton_root_avr + cton_root_min) +
                     this.cton_sap  / (this.cton_sap_avr  + cton_sap_min)
  let this = this with cton_stem_max = 1.0/(2.0*0.0034) --Maize params
  in this with cton_stem_avr = 1.0/(2.0*0.0068)

-- Calculates coefficient to compensate for different vertical distribution of fine root on nitrogen uptake
let init_nupscoeff(this: Pft) : Pft =
  -- Fraction fine root in upper soil layer should have higher possibility for mineralized nitrogen uptake
  -- Soil nitrogen profile is considered to have a exponential decline (Franzluebbers et al. 2009) giving
  -- an approximate advantage of 2 of having more roots in the upper soil layer
  let upper_adv = 2.0

  -- Simple solution until we get C and N in all soil layers.
  let (_, rootdist_upper, rootdist_lower) =
  loop (sl, rootdist_upper, rootdist_lower) = (0, 0.0, 0.0)
  while (sl < NSOILLAYER) do
    if (sl < NSOILLAYER_UPPER) then
      (sl+1, rootdist_upper + this.rootdist[sl], rootdist_lower) else
      (sl+1, rootdist_upper, rootdist_lower + this.rootdist[sl])
  in this with nupscoeff = rootdist_upper * upper_adv + rootdist_lower


-- Initialises sapling/regen characteristics in population mode following LPJF formulation
let initregen(this: Pft) : Pft =

  -- see function allometry in growth module.

  -- Note: primary PFT parameters, including SLA, must be set before this
  --       function is called

  let REGENLAI_TREE = 1.5
  let REGENLAI_GRASS = 0.001
  let SAPLINGHW = 0.2

  let this = if (this.lifeform == TREE) then

    -- Tree sapling characteristics

    let this = this with regen_cmass_leaf =
        pow(REGENLAI_TREE * this.k_allom1 * pow(1.0 + SAPLINGHW, this.k_rp) * pow(4.0 * this.sla / PI / this.k_latosa, this.k_rp * 0.5) / this.sla, 2.0 / (2.0 - this.k_rp))


    let this = this with regen_cmass_leaf =
        this.wooddens * this.k_allom2 * pow((1.0 + SAPLINGHW) * sqrt(4.0 * this.regen_cmass_leaf * this.sla / PI / this.k_latosa), this.k_allom3) * this.regen_cmass_leaf * this.sla / this.k_latosa

    in this with regen_cmass_heart = SAPLINGHW * this.regen_cmass_sap

  else if (this.lifeform == GRASS || this.lifeform == MOSS) then
    -- Grass regeneration characteristics
    this with regen_cmass_leaf = REGENLAI_GRASS / this.sla
    else this

  in this with regen_cmass_root = 1.0 / this.ltor_max * this.regen_cmass_leaf


-- Inits root fractions in each soil layer through a shape parameter beta (see Jackson et al., 1996)
let init_rootdist(this: Pft) : Pft =

  let depth = Dz_soil * CM_PER_MM

  let rootdist = copy this.rootdist with [0] = 1.0 - pow(this.root_beta, depth) -- init first layer

  let tot = rootdist[0]

  let (_, _, rootdist, tot) =
  loop (i, depth, rootdist, tot) = (1, depth, rootdist, tot) while (i<NSOILLAYER) do
    let depth = depth + (Dz_soil * CM_PER_MM)
    let rootdist[i] =  1.0 - pow(this.root_beta, depth) - (1.0 - pow(this.root_beta, depth - Dz_soil * CM_PER_MM))
    let tot = tot + rootdist[i]
    in (i+1, depth, rootdist, tot)

    -- Calibrated the root_beta for each PFT to match rootdist_upper from 'old' (pre LPJG 4.1) parameterisation.
    -- Sometimes the rootdist goes below our maximum soildepth. When that happens, put the residual fraction in lowest soil layer
    let rootdist[NSOILLAYER-1] = rootdist[NSOILLAYER-1] + 1.0 - tot
    in this with rootdist = rootdist

let ismoss(this: Pft) : bool = this.lifeform == MOSS
let isgrass(this: Pft) : bool = this.lifeform == GRASS
let istree(this: Pft) : bool = this.lifeform == TREE
let iswetlandspecies(this: Pft) : bool = (this.lifeform == MOSS || this.has_aerenchyma)

-- type cropindiv
-- type Individual
-- type Vegetation


--- Soiltype stores static parameters for soils and the snow pack.
--- One Soiltype object is defined for each Gridcell. State variables for soils
--  are held by objects of class Soil, of which there is one for each patch
--  (see below).

type Soiltype = {

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

let Soiltype() : Soiltype = {
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
let updateSolveSOMvalues(this : Soiltype, nyrspinup : int) =
  let this = this with solvesom_end = intFromReal (0.8 * (realFromInt nyrspinup))
  let this = this with solvesom_begin = intFromReal (0.7 * (realFromInt nyrspinup))
  in this

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


--- This struct contains litter for solving Century SOM pools.
-- We dont need getters, just read the field.
type LitterSolveSOM = {
  --- Carbon litter
  clitter : [NSOMPOOL]real,
  --- Nitrogen litter
  nlitter : [NSOMPOOL]real
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
  --- fraction of available water holding capacity */
  wcont_evap : real,

  --- reference to parent Patch object
  --patch : Patch, -- NOTE we dont do references in futhark
  --- reference to Soiltype object holding static parameters for this soil
  soiltype : Soiltype,
  --- the average wcont over the growing season, for each of the upper soil layers. Used in drought limited establishment.
  awcont_upper : real,
  --- daily water content in upper soil layer for each day of year
  dwcontupper : [Date_MAX_YEAR_LENGTH]real,
  --- mean water content in upper soil layer for last month
  --- (valid only on last day of month following call to daily_accounting_patch) */
  mwcontupper : real,
  --- stored snow as average over modelled area (mm rainfall equivalent)
  snowpack : real,
  --- total runoff today (mm/day)
  runoff : real,
  --- daily temperatures for the last month (deg C)
  --- (valid only on last day of month following call to daily_accounting_patch) */
  dtemp : [31]real,
  --- mean soil temperature for the last month (deg C)
  --- (valid only on last day of month following call to daily_accounting_patch) */
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
  --- (valid only on last day of month following call to daily_accounting_patch) */
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
  -- (awc + Fpwp) volumetric water content MINUS the ice fraction.
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
  --- gaseous CH4 concentration in each layer [g CH4-C layer-1] yesterday
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


let Soil(soiltype : Soiltype) : Soil = {
  soiltype = soiltype,

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


let cropphen_struct() : cropphen_struct = {
    sdate=(-1),
    sdate_harv=(-1),
    nsow=0,
    sownlastyear=false,
    sendate=(-1),
    hdate=(-1),
    hlimitdate=(-1),
    hucountend=(-1),
    nharv=0,
    tb=0.0,
    pvd=0,
    vdsum=0,
    vrf=1.0,
    phu=0.0,
    phu_old=0.0,
    husum_max=0.0,
    husum_sampled=0.0,
    husum_max_10=0.0,
    nyears_hu_sample = 0,
    prf=1.0,
    husum=0.0,
    fphu=0.0,
    fphu_harv=0.0,
    hu_samplingdays=0,
    hu_samplingperiod=false,

    hi=0.0,
    fhi=0.0,
    fhi_phen=0.0,
    fhi_water=1.0,
    fhi_harv=0.0,
    demandsum_crop=0.0,
    supplysum_crop=0.0,

    growingseason=false,  --Initialized to true for normal grass growth (CC3G & CC4G) in establishment
    growingseason_ystd=false,
    senescence=false,
    senescence_ystd=false,
    intercropseason=false,
    bicdate=(-1),
    eicdate=(-1),
    growingdays=0,
    growingdays_y=0,
    lgp=0,

    sdate_harvest = replicate 2 (-1),
    hdate_harvest = replicate 2 (-1),
    sdate_thisyear = replicate 2 (-1),

    vdsum_alloc=0.0,
    vd = 0.0,
    f_alloc_root=0.0,
    f_alloc_leaf=0.0,
    f_alloc_horg=0.0,
    f_alloc_stem=0.0,
    dev_stage = 0.0,

    fertilised = replicate 3 false
  }


--/ State variables common to all individuals of a particular PFT in a particular patch
-- Used in individual and cohort modes only. */
type Patchpft = {
  -- MEMBER VARIABLES:
  --/ id code (equal to value of member variable id in corresponding Pft object)
  id : int,
  --/ reference to corresponding Pft object in PFT list
  pft : Pft,
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


-- Constructor: initialises id, pft and data members
let Patchpft(i: int, p: Pft) : Patchpft = {
  id = i,
  pft = p,

  litter_leaf = 0.0,
  litter_root = 0.0,
  litter_sap   = 0.0,
  litter_heart = 0.0,
  litter_repr = 0.0,

  nmass_litter_leaf  = 0.0,
  nmass_litter_root  = 0.0,
  nmass_litter_sap   = 0.0,
  nmass_litter_heart = 0.0,

  wscal = 1.0,
  wscal_mean = 1.0,
  anetps_ff = 0.0,
  aphen = 0.0,
  phen = 0.0,
  wsupply = 0.0,
  wsupply_leafon = 0.0,
  anetps_ff_est = 0.0,
  anetps_ff_est_initial = 0.0,
  wscal_mean_est = 0.0,
  nsapling = 0,

  fwuptake = replicate NSOILLAYER 0.0,

  cropphen = cropphen_struct(),--NULL, -- We cant have this.
  harvested_products_slow = 0.0,
  harvested_products_slow_nmass = 0.0,

  swindow = replicate 2 (-1), -- TODO: this is inefficient

  --if(pft.landcover==CROPLAND)
  --{
  --  cropphen=new cropphen_struct,
  --}

  inund_count=0,
  inund_stress=1.0, -- No stress by default
  -- TODO: THE FOLLOWING FIELDS WERE UNINITIALIZED IN THE C++ CODE
  establish=false,
  gcbase=nan,
  gcbase_day=nan,
  water_deficit_d=nan,
  water_deficit_y=nan,
  wstress=false,
  wstress_day=false
}


-- Stores data for a patch.
-- In cohort and individual modes, replicate patches are
--  required in each stand to accomodate stochastic variation, in population mode there
--  should be just one Patch object, representing average conditions for the entire
--  stand. A reference to the parent Stand object (defined below) is included as a
--  member variable.

type Patch = {

    -- id code in range 0-npatch for patch
    id : int,
    -- reference to parent Stand object
    stand: Stand,
    -- list array [0...npft-1] of Patchpft objects (initialised in constructor)
    pfts: [npft]Patchpft,
    --ListArray_idin1<Patchpft,Pft> pft, -- In the c++ code, this was called "pft", even though it is multiple. Should be called "pfts".
    -- vegetation for this patch
    vegetation: Vegetation,
    -- soil for this patch
    soil: Soil,
    -- fluxes for this patch
    fluxes: Fluxes,
    -- FPAR at top of grass canopy today
    fpar_grass: real,
    -- FPAR at soil surface today
    fpar_ff: real,
    -- mean growing season PAR at top of grass canopy (J/m2/day)
    par_grass_mean: real,
    -- number of days in growing season, estimated from mean vegetation leaf-on fraction
    -- \see function fpar in canopy exchange module */
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
    -- i.e. those days for which daily temp > 5.0 degC */
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
    -- mm/day, patch vegetative area basis  */
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
