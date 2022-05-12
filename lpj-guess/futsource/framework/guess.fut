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
type lifeformtype = #NOLIFEFORM | #TREE | #GRASS | #MOSS

-- Phenology class for PFTs
type phenologytype = #NOPHENOLOGY | #EVERGREEN | #RAINGREEN | #SUMMERGREEN | #CROPGREEN | #ANY

-- Biochemical pathway for photosynthesis (C3 or C4)
type pathwaytype = #NOPATHWAY | #C3 | #C4

-- Leaf physiognomy types for PFTs
type leafphysiognomytype = #NOLEAFTYPE | #NEEDLELEAF | #BROADLEAF

-- The level of verbosity of LPJ-GUESS. Decides the amount of information that is written to the log-file.
type verbositylevel = #ERROR | #WARNING | #INFO | #DEBUG_WARNING

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
type insoltype =
  #NOINSOL -- No insolation type chosen
  | #SUNSHINE -- Percentage sunshine
  | #NETSWRAD -- Net shortwave radiation flux during daylight hours (W/m2)
  | #SWRAD -- Total shortwave radiation flux during daylight hours (W/m2)
  | #NETSWRAD_TS -- Net shortwave radiation flux during whole time step (W/m2)
  | #SWRAD_TS -- Total shortwave radiation flux during whole time step (W/m2)

-- CENTURY pool names, NSOMPOOL number of SOM pools
type pooltype =
  #SURFSTRUCT
  | #SOILSTRUCT
  | #SOILMICRO
  | #SURFHUMUS
  | #SURFMICRO
  | #SURFMETA
  | #SURFFWD
  | #SURFCWD
  | #SOILMETA
  | #SLOWSOM
  | #PASSIVESOM
  | #NSOMPOOL


-- Irrigation type for PFTs
type hydrologytype = #RAINFED | #IRRIGATED

-- Intercrop type for PFTs
type intercroptype = #NOINTERCROP | #NATURALGRASS

-- Seasonality type of gridcell
--  0:SEASONALITY_NO      No seasonality
--  1:SEASONALITY_PREC      Precipitation seasonality only
--  2:SEASONALITY_PRECTEMP    Both temperature and precipitation seasonality, but "weak" temperature seasonality (coldest month > 10degC)
--  3:SEASONALITY_TEMP      Temperature seasonality only
--  4:SEASONALITY_TEMPPREC    Both temperature and precipitation seasonality, but temperature most important (coldest month < 10degC)
--  5:    Temperature seasonality, always above 10 degrees (currently not used)
--let SEASONALITY_TEMPWARM : seasonality_type = 5 -- TODO: this one was in the original comment but not original source
type TYPENAME = #SEASONALITY_NO | #SEASONALITY_PREC | #SEASONALITY_PRECTEMP | #SEASONALITY_TEMP | #SEASONALITY_TEMPPREC

-- Precipitation seasonality type of gridcell
-- 0:DRY            (minprec_pet20<=0.5 && maxprec_pet20<=0.5)
--  1:DRY_INTERMEDIATE      (minprec_pet20<=0.5 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
--  2:DRY_WET          (minprec_pet20<=0.5 && maxprec_pet20>1.0)
--  3:INTERMEDIATE        (minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
--  4:INTERMEDIATE_WET      (minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>1.0)
--  5:WET            (minprec_pet20>1.0 && maxprec_pet20>1.0)
--/
type prec_seasonality_type = #DRY | #DRY_INTERMEDIATE | #DRY_WET | #INTERMEDIATE | #INTERMEDIATE_WET | #WET

-- Temperature seasonality type of gridcell
-- 0:COLD            (mtemp_max20<=10)
-- 1:COLD_WARM          (mtemp_min20<=10 && mtemp_max20>10 && mtemp_max20<=30)
-- 2:COLD_HOT          (mtemp_min20<=10 && mtemp_max20>30)
-- 3:WARM            (mtemp_min20>10 && mtemp_max20<=30)
-- 4:WARM_HOT          (mtemp_min20>10 && mtemp_max20>30)
-- 5:HOT            (mtemp_min20>30)
--/
type temp_seasonality_type = #COLD | #COLD_WARM | #COLD_HOT | #WARM | #WARM_HOT | #HOT | #SIX | #SEVEN

-- Gas type (used in methane code)
type gastype = #O2gas | #CO2gas | #CH4gas

-- Nitrogen preferance
type n_pref_type = #NO | #NH4 | #NO3

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
    landcover = #URBAN,
    leaflong = nan,
    leafphysiognomy = #NOLEAFTYPE,
    lifeform = #NOLIFEFORM,

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
    pathway = #C4,

    phengdd5ramp = nan,
    phenology = #CROPGREEN,
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
  if (!(this.phenology == #CROPGREEN && this.nlim)) then
    -- Reich et al 1992, Table 1 (includes conversion x2.0 from m2/kg_dry_weight to
    -- m2/kgC)
    this with sla =
    if (this.leafphysiognomy == #BROADLEAF) then
      0.2 * pow(10.0, 2.41 - 0.38 * log10(12.0 * this.leaflong))
    else if (this.leafphysiognomy == #NEEDLELEAF) then
      0.2 * pow(10.0, 2.29 - 0.4 * log10(12.0 * this.leaflong))
    else this.sla
  else this


-- Calculates minimum leaf C:N ratio given leaf longevity
let init_cton_min(this: Pft) : Pft =
  -- cton_leaf_min has to be supplied in the insfile for crops with N limitation
  if (!(this.phenology == #CROPGREEN && this.nlim)) then
    -- Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
    -- kgN/kgC)

    this with cton_leaf_min =
    if (this.leafphysiognomy == #BROADLEAF) then
      500.0 / pow(10.0, 1.75 - 0.33 * log10(12.0 * this.leaflong))
    else if (this.leafphysiognomy == #NEEDLELEAF) then
      500.0 / pow(10.0, 1.52 - 0.26 * log10(12.0 * this.leaflong))
    else this.cton_leaf_min
  else this

let init_cton_limits(this: Pft) : Pft =
  -- Fraction between min and max C:N ratio White et al. 2000
  let frac_mintomax = if (this.phenology == #CROPGREEN && this.nlim) then 5.0 else 2.78  -- Use value also without nlim ?

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
  if (this.lifeform == #GRASS || this.lifeform == #MOSS) then
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

  let this = if (this.lifeform == #TREE) then

    -- Tree sapling characteristics

    let this = this with regen_cmass_leaf =
        pow(REGENLAI_TREE * this.k_allom1 * pow(1.0 + SAPLINGHW, this.k_rp) * pow(4.0 * this.sla / PI / this.k_latosa, this.k_rp * 0.5) / this.sla, 2.0 / (2.0 - this.k_rp))


    let this = this with regen_cmass_leaf =
        this.wooddens * this.k_allom2 * pow((1.0 + SAPLINGHW) * sqrt(4.0 * this.regen_cmass_leaf * this.sla / PI / this.k_latosa), this.k_allom3) * this.regen_cmass_leaf * this.sla / this.k_latosa

    in this with regen_cmass_heart = SAPLINGHW * this.regen_cmass_sap

  else if (this.lifeform == #GRASS || this.lifeform == #MOSS) then
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

let isgrass(this: Pft) : bool = this.lifeform == #GRASS
let istree(this: Pft) : bool = this.lifeform == #TREE
let iswetlandspecies(this: Pft) : bool = (this.lifeform == #MOSS || this.has_aerenchyma)
