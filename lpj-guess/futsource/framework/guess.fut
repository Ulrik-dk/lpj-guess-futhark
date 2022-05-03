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
open import "../framework/guessmath.h"
--#include "archive.h"
open import "../framework/parameters.h"
--#include "guesscontainer.h"
open import "../modules/soil.h"


----------------------------------------------------------
-- GLOBAL ENUMERATED TYPE DEFINITIONS

-- Life form class for PFTs (trees, grasses)
type lifeformtype = i8
let NOLIFEFORM : lifeformtype = 0
let TREE : lifeformtype = 1
let GRASS : lifeformtype = 2
let MOSS : lifeformtype = 3


-- Phenology class for PFTs
type phenologytype = i8
let NOPHENOLOGY : phenologytype = 0
let EVERGREEN : phenologytype = 1
let RAINGREEN : phenologytype = 2
let SUMMERGREEN : phenologytype = 3
let CROPGREEN : phenologytype = 4
let ANY : phenologytype = 5

-- Biochemical pathway for photosynthesis (C3 or C4)
type pathwaytype = i8
let NOPATHWAY : pathwaytype = 0
let C3 : pathwaytype = 1
let C4 : pathwaytype = 2

-- Leaf physiognomy types for PFTs
type leafphysiognomytype = i8
let NOLEAFTYPE : leafphysiognomytype = 0
let NEEDLELEAF : leafphysiognomytype = 1
let BROADLEAF : leafphysiognomytype = 2

-- The level of verbosity of LPJ-GUESS. Decides the amount of information that is written to the log-file.
type verbositylevel = i8
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
type insoltype = i8
let NOINSOL : insoltype = 0 -- No insolation type chosen
let SUNSHINE : insoltype = 1 -- Percentage sunshine
let NETSWRAD : insoltype = 2 -- Net shortwave radiation flux during daylight hours (W/m2)
let SWRAD : insoltype = 3 -- Total shortwave radiation flux during daylight hours (W/m2)
let NETSWRAD_TS : insoltype = 4 -- Net shortwave radiation flux during whole time step (W/m2)
let SWRAD_TS : insoltype = 5 -- Total shortwave radiation flux during whole time step (W/m2)

-- CENTURY pool names, NSOMPOOL number of SOM pools
type pooltype = i8
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
type hydrologytype = i8
let RAINFED : hydrologytype = 0
let IRRIGATED : hydrologytype = 1

-- Intercrop type for PFTs
type intercroptype = i8
let NOINTERCROP : intercroptype = 0
let NATURALGRASS : intercroptype = 1

-- Seasonality type of gridcell
--  0:SEASONALITY_NO			No seasonality
--  1:SEASONALITY_PREC			Precipitation seasonality only
--  2:SEASONALITY_PRECTEMP		Both temperature and precipitation seasonality, but "weak" temperature seasonality (coldest month > 10degC)
--  3:SEASONALITY_TEMP			Temperature seasonality only
--  4:SEASONALITY_TEMPPREC		Both temperature and precipitation seasonality, but temperature most important (coldest month < 10degC)
--  5:		Temperature seasonality, always above 10 degrees (currently not used)
type seasonality_type = i8
let SEASONALITY_NO : seasonality_type = 0
let SEASONALITY_PREC : seasonality_type = 1
let SEASONALITY_PRECTEMP : seasonality_type = 2
let SEASONALITY_TEMP : seasonality_type = 3
let SEASONALITY_TEMPPREC : seasonality_type = 4
--let SEASONALITY_TEMPWARM : seasonality_type = 5 -- TODO: this one was in the original comment but not source

-- Precipitation seasonality type of gridcell
-- 0:DRY						(minprec_pet20<=0.5 && maxprec_pet20<=0.5)
--  1:DRY_INTERMEDIATE			(minprec_pet20<=0.5 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
--  2:DRY_WET					(minprec_pet20<=0.5 && maxprec_pet20>1.0)
--  3:INTERMEDIATE				(minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
--  4:INTERMEDIATE_WET			(minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>1.0)
--  5:WET						(minprec_pet20>1.0 && maxprec_pet20>1.0)
--/
type prec_seasonality_type = i8
let DRY : prec_seasonality_type = 0
let DRY_INTERMEDIATE : prec_seasonality_type = 1
let DRY_WET : prec_seasonality_type = 2
let INTERMEDIATE : prec_seasonality_type = 3
let INTERMEDIATE_WET : prec_seasonality_type = 4
let WET : prec_seasonality_type = 5


-- Temperature seasonality type of gridcell
-- 0:COLD						(mtemp_max20<=10)
-- 1:COLD_WARM					(mtemp_min20<=10 && mtemp_max20>10 && mtemp_max20<=30)
-- 2:COLD_HOT					(mtemp_min20<=10 && mtemp_max20>30)
-- 3:WARM						(mtemp_min20>10 && mtemp_max20<=30)
-- 4:WARM_HOT					(mtemp_min20>10 && mtemp_max20>30)
-- 5:HOT						(mtemp_min20>30)
--/
type temp_seasonality_type = i8
let COLD : temp_seasonality_type = 0
let COLD_WARM : temp_seasonality_type = 1
let COLD_HOT : temp_seasonality_type = 2
let WARM : temp_seasonality_type = 3
let WARM_HOT : temp_seasonality_type = 4
let HOT : temp_seasonality_type = 5

-- Gas type (used in methane code)
--
type gastype = i8
let O2gas : gastype = 0
let CO2gas : gastype = 1
let CH4gas : gastype = 2

-- Nitrogen preferance
type n_pref_type = i8
let NO : n_pref_type = 0
let NH4 : n_pref_type = 1
let NO3 : n_pref_type = 2

----------------------------------------------------------
-- GLOBAL CONSTANTS

-- number  of soil layers modelled
let NSOILLAYER_UPPER : i64 = 5
let NSOILLAYER_LOWER : i64 = NSOILLAYER - NSOILLAYER_UPPER

-- bvoc: number of monoterpene species used
let NMTCOMPOUNDS : i64 = NMTCOMPOUNDTYPES

-- SOIL DEPTH VALUES

-- soil upper layer depth (mm)
let SOILDEPTH_UPPER : f64 = 500.0
-- soil lower layer depth (mm)
let SOILDEPTH_LOWER : f64 = 1000.0

-- Depth of sublayer at top of upper soil layer, from which evaporation is
-- possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
-- Must be a multiple of Dz_soil
let SOILDEPTH_EVAP : f64 = 200.0

-- Year at which to calculate equilibrium soil carbon
let SOLVESOM_END : i64 = 400

-- Year at which to begin documenting means for calculation of equilibrium soil carbon
let SOLVESOM_BEGIN : i64 = 350

-- Number of years to average growth efficiency over in function mortality
let NYEARGREFF : i64 = 5

-- Coldest day in N hemisphere (January 15)
-- Used to decide when to start counting GDD's and leaf-on days
--  for summergreen phenology.
--/
let COLDEST_DAY_NHEMISPHERE : i64 = 14

-- Coldest day in S hemisphere (July 15)
-- Used to decide when to start counting GDD's and leaf-on days
--  for summergreen phenology.
--/
let COLDEST_DAY_SHEMISPHERE : i64 = 195

-- Warmest day in N hemisphere (same as COLDEST_DAY_SHEMISPHERE)
let WARMEST_DAY_NHEMISPHERE : i64 = COLDEST_DAY_SHEMISPHERE

-- Warmest day in S hemisphere (same as COLDEST_DAY_NHEMISPHERE)
let WARMEST_DAY_SHEMISPHERE : i64 = COLDEST_DAY_NHEMISPHERE

-- number of years to average aaet over in function soilnadd
let NYEARAAET : i64 = 5

-- number of years to average max snow depth over in function soilnadd
let NYEARMAXSNOW : i64 = 20

-- Priestley-Taylor coefficient (conversion factor from equilibrium evapotranspiration to PET)
let PRIESTLEY_TAYLOR : f64 = 1.32

-- Solving Century SOM pools

-- fraction of nyear_spinup minus freenyears at which to begin documenting for calculation of Century equilibrium
let SOLVESOMCENT_SPINBEGIN  : f64 = 0.1
-- fraction of nyear_spinup minus freenyears at which to end documentation and start calculation of Century equilibrium
let SOLVESOMCENT_SPINEND : f64 = 0.3

-- Kelvin to deg C conversion
let K2degC : f64 = 273.15

-- Maximum number of crop rotation items
let NROTATIONPERIODS_MAX : i64 = 3

-- Conversion factor for CO2 from ppmv to mole fraction
let CO2_CONV : f64 = 1.0e-6

-- Initial carbon allocated to crop organs at sowing, kg m-2
let CMASS_SEED : f64 = 0.01

-- Precision in land cover fraction input
let INPUT_PRECISION : f64 = 1.0e-14
let INPUT_ERROR : f64 = 0.5e-6
let INPUT_RESOLUTION : f64 = INPUT_PRECISION - INPUT_PRECISION * INPUT_ERROR

-- Averaging interval for average maximum annual fapar (SIMFIRE)
let AVG_INTERVAL_FAPAR : i64 = 3

-- Averaging interval for biome averaging (SIMFIRE)
let N_YEAR_BIOMEAVG : i64 = 3

















-- This struct contains the result of a photosynthesis calculation.
-- see photosynthesis
type PhotosynthesisResult = {
	-- Constructs an empty result

	-- RuBisCO capacity (gC/m2/day)
	vm : f64,

	-- gross daily photosynthesis (gC/m2/day)
	agd_g : f64,

	-- leaf-level net daytime photosynthesis
	-- expressed in CO2 diffusion units (mm/m2/day) */
	adtmm : f64,

	-- leaf respiration (gC/m2/day)
	rd_g : f64,

	-- PAR-limited photosynthesis rate (gC/m2/h)
	je : f64,

	-- optimal leaf nitrogen associated with photosynthesis (kgN/m2)
	nactive_opt : f64,

	-- nitrogen limitation on vm
	vmaxnlim : f64

	-- net C-assimilation (gross photosynthesis minus leaf respiration) (kgC/m2/day)
}

let net_assimilation(psr : PhotosynthesisResult) : f64 =
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
	co2 : f64,
	-- mean air temperature today (deg C)
	temp : f64,
	-- total daily photosynthetically-active radiation today (J / m2 / day) (ALPHAA not yet accounted for)
	par : f64,
	-- fraction of PAR absorbed by foliage
	fpar : f64,
	-- day length, must equal 24 in diurnal mode(h)
	daylength : f64
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

-- This struct contains the stresses used in a photosynthesis calculation.
-- \see photosynthesis */
type PhotosynthesisStresses = {
  -- whether nitrogen should limit Vmax
	ifnlimvmax : bool,
	--  limit to moss photosynthesis. [0,1], where 1 means no limit
	moss_ps_limit : f64,
	-- limit to graminoid photosynthesis. [0,1], where 1 means no limit
	graminoid_ps_limit : f64,
	-- limit to photosynthesis due to inundation, where 1 means no limit
	inund_stress : f64
}

let PhotosynthesisStresses() : PhotosynthesisStresses =
  -- All members set to no stress values
	-- Default values indicating no stress
  {
    ifnlimvmax = false,
    moss_ps_limit = 1.0,
    graminoid_ps_limit = 1.0,
    inund_stress = 1.0
  }


-- Holds static functional parameters for a plant functional type (PFT).
-- There should be one Pft object for each potentially occurring PFT. The same Pft object
-- may be referenced (via the pft member of the Individual object see below) by different
-- average individuals. Member functions are included for initialising SLA given leaf
-- longevity, and for initialising sapling/regen characteristics (required for
-- population mode).

--module Pft = {
  type Pft  = {
    -- MEMBER VARIABLES
    -- id code (should be zero based and sequential, 0...npft-1)
    id: i64,
    -- name of PFT
    name: xtring,
    -- life form (tree or grass)
    lifeform: lifeformtype,
    -- leaf phenology (raingreen, summergreen, evergreen, rain+summergreen, cropgreen)
    phenology: phenologytype,
    -- leaf physiognomy (needleleaf, broadleaf)
    leafphysiognomy: leafphysiognomytype,
    -- growing degree sum on 5 degree base required for full leaf cover
    phengdd5ramp: f64,
    -- water stress threshold for leaf abscission (range 0-1 raingreen PFTs)
    wscal_min: f64,
    -- biochemical pathway for photosynthesis (C3 or C4)
    pathway: pathwaytype,
    -- approximate low temperature limit for photosynthesis (deg C)
    pstemp_min: f64,
    -- approximate lower range of temperature optimum for photosynthesis (deg C)
    pstemp_low: f64,
    -- approximate upper range of temperature optimum for photosynthesis (deg C)
    pstemp_high: f64,
    -- maximum temperature limit for photosynthesis (deg C)
    pstemp_max: f64,
    -- non-water-stressed ratio of i64ercellular to ambient CO2 partial pressure
    lambda_max: f64,
    -- vegetation root profile in an array containing fraction of roots in each soil layer, [0=upper layer]
    rootdist: [NSOILLAYER]f64,
    -- shape parameter for initialisation of root distribtion
    root_beta: f64,
    -- canopy conductance component not associated with photosynthesis (mm/s)
    gmin: f64,
    -- maximum evapotranspiration rate (mm/day)
    emax: f64,
    -- mai64enance respiration coefficient (0-1)
    respcoeff: f64,

    -- minimum leaf C:N mass ratio allowed when nitrogen demand is determined
    cton_leaf_min: f64,
    -- maximum leaf C:N mass ratio  allowed when nitrogen demand is determined
    cton_leaf_max: f64,
    -- average leaf C:N mass ratio (between min and max)
    cton_leaf_avr: f64,
    -- average fine root C:N mass ratio (connected cton_leaf_avr)
    cton_root_avr: f64,
    -- maximum fine root C:N mass ratio (used when mass is negligible)
    cton_root_max: f64,
    -- average sapwood C:N mass ratio (connected cton_leaf_avr)
    cton_sap_avr: f64,
    -- maximum sapwood C:N mass ratio (used when mass is negligible)
    cton_sap_max: f64,
    -- reference fine root C:N mass ratio
    cton_root: f64,
    -- reference sapwood C:N mass ratio
    cton_sap: f64,
    -- Maximum nitrogen (NH4+ and NO3- seperatly) uptake per fine root [kgN kgC-1 day-1]
    nuptoroot: f64,
    -- coefficient to compensate for vertical distribution of fine root on nitrogen uptake
    nupscoeff: f64,
    -- fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen longterm storage scalar
    fnstorage: f64,

    -- Michaelis-Menten kinetic parameters
    -- Half saturation concentration for N uptake [kgN l-1] (Rothstein 2000) */
    km_volume: f64,

    -- fraction of NPP allocated to reproduction
    reprfrac: f64,
    -- annual leaf turnover as a proportion of leaf C biomass
    turnover_leaf: f64,
    -- annual fine root turnover as a proportion of fine root C biomass
    turnover_root: f64,
    -- annual sapwood turnover as a proportion of sapwood C biomass
    turnover_sap: f64,
    -- sapwood and heartwood density (kgC/m3)
    wooddens: f64,
    -- maximum tree crown area (m2)
    crownarea_max: f64,
    -- constant in allometry equations
    k_allom1: f64,
    -- constant in allometry equations
    k_allom2: f64,
    -- constant in allometry equations
    k_allom3: f64,
    -- constant in allometry equations
    k_rp: f64,
    -- tree leaf to sapwood area ratio
    k_latosa: f64,
    -- specific leaf area (m2/kgC)
    sla: f64,
    -- leaf longevity (years)
    leaflong: f64,
    -- leaf to root mass ratio under non-water-stressed conditions
    ltor_max: f64,
    -- litter moisture flammability threshold (fraction of AWC)
    litterme: f64,
    -- fire resistance (0-1)
    fireresist: f64,
    -- minimum forest-floor PAR level for growth (grasses) or establishment (trees)
    -- J/m2/day, individual and cohort modes */
    parff_min: f64,
    -- parameter capturing non-linearity in recruitment rate relative to
    -- understorey growing conditions for trees (Fulton 1991) (individual and
    -- cohort modes)

    alphar: f64,
    -- maximum sapling establishment rate (saplings/m2/year) (individual and cohort modes)
    est_max: f64,
    -- constant used in calculation of sapling establishment rate when spatial
    -- mass effect enabled (individual and cohort modes)

    kest_repr: f64,
    -- constant affecting amount of background establishment
    -- \see ifbgestab */
    kest_bg: f64,
    -- constant used in calculation of sapling establishment rate when spatial
    -- mass effect disabled (individual and cohort modes)

    kest_pres: f64,
    -- expected longevity under non-stressed conditions (individual and cohort modes)
    longevity: f64,
    -- threshold growth efficiency for imposition of growth suppression mortality
    -- kgC/m2 leaf/year, individual and cohort modes */
    greff_min: f64,

    -- Bioclimatic limits (all temperatures deg C)

    -- minimum 20-year coldest month mean temperature for survival
    tcmin_surv: f64,
    -- maximum 20-year coldest month mean temperature for establishment
    tcmax_est: f64,
    -- minimum degree day sum on 5 deg C base for establishment
    gdd5min_est: f64,
    -- minimum 20-year coldest month mean temperature for establishment
    tcmin_est: f64,
    -- minimum warmest month mean temperature for establishment
    twmin_est: f64,
    -- continentality parameter for boreal summergreen trees
    twminusc: f64,
    -- constant in equation for budburst chilling time requirement (Sykes et al 1996)
    k_chilla: f64,
    -- coefficient in equation for budburst chilling time requirement
    k_chillb: f64,
    -- exponent in equation for budburst chilling time requirement
    k_chillk: f64,
    -- array containing values for GDD0(c) given c=number of chill days
    -- Sykes et al 1996, Eqn 1
    -- gdd0 has one element for each possible value for number of chill days

    ---- FIXME TODO HACK! comes from Date::MAX_YEAR_LENGTH+1
    gdd0: [Date_MAX_YEAR_LENGTH_plusone]f64,

    -- i64erception coefficient (unitless)
    i64c: f64,

    -- the amount of N that is applied (kg N m-2)
    N_appfert: f64,
    -- 0 - 1 how much of the fertiliser is applied the first date, default 1.
    fertrate: (f64, f64),
    -- dates relative to sowing date
    fertdates: (i64, i64),
    fert_stages: (f64, f64),
    fertilised: (bool, bool),

    T_vn_min: f64,
    T_vn_opt: f64,
    T_vn_max: f64,

    T_veg_min: f64,
    T_veg_opt: f64,
    T_veg_max: f64,

    T_rep_min: f64,
    T_rep_opt: f64,
    T_rep_max: f64,

    photo: (f64, f64, f64),

    dev_rate_veg: f64,
    dev_rate_rep: f64,

    a1: f64, b1: f64, c1: f64, d1: f64, a2: f64, b2: f64, c2: f64, d2: f64, a3: f64, b3: f64, c3: f64, d3: f64,
    cton_stem_avr: f64,
    cton_stem_max: f64,

    -- Drought tolerance level (0 = very -> 1 = not at all) (unitless)
    -- Used to implement drought-limited establishment */
    drought_tolerance: f64,

    -- bvoc

    -- aerodynamic conductance (m s-1)
    ga: f64,
    -- isoprene emission capacity (ug C g-1 h-1)
    eps_iso: f64,
    -- whether (1) or not (1) isoprene emissions show a seasonality
    seas_iso: bool,
    -- monoterpene emission capacity (ug C g-1 h-1) per monoterpene species
    eps_mon: [NMTCOMPOUNDS]f64,
    -- fraction of monoterpene production that goes i64o storage pool (-) per monoterpene species
    storfrac_mon: [NMTCOMPOUNDS]f64,

    -- Bioclimatic limits parameters from Wolf et al. 2008

    -- snow max [mm]
    max_snow: f64,
    -- snow min [mm]
    min_snow: f64,
    -- GDD0 min
    gdd0_min: f64,
    -- GDD0 max
    gdd0_max: f64,

    -- New parameters from parameters from Wania et al. (2009a, 2009b, 2010)

    -- Days per month for which inundation is tolerated
    inund_duration: i64,
    -- Inundation stress is felt when the water table (mm) is above wtp_max
    wtp_max: f64,
    -- Whether this PFT has aerenchyma through which O2 and CH4 can be transported (Wania et al. 2010 - Sec 2.6)
    has_aerenchyma: bool,

    -- Sapling/regeneration characteristics (used only in population mode)
    -- For trees, on sapling individual basis (kgC) for grasses, on stand area basis,
    -- kgC/m2 */

    -- leaf C biomass
    regen_cmass_leaf: f64,
    -- fine root C biomass
    regen_cmass_root: f64,
    -- sapwood C biomass
    regen_cmass_sap: f64,
    -- heartwood C biomass
    regen_cmass_heart: f64,

    -- specifies type of landcover pft is allowed to grow in (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND)
    landcover: landcovertype,
    -- pft selection
    selection: xtring,
    -- fraction of residue outtake at harvest
    res_outtake: f64,
    -- harvest efficiencytype leafphysiognomy = i64

    harv_eff: f64,
    -- harvest efficiency of i64ercrop grass
    harv_eff_ic: f64,
    -- fraction of harvested products that goes i64o patchpft.harvested_products_slow
    harvest_slow_frac: f64,
    -- yearly turnover fraction of patchpft.harvested_products_slow (goes to gridcell.acflux_harvest_slow)
    turnover_harv_prod: f64,
    -- whether pft may grow as cover crop
    isi64ercropgrass: bool,
    -- whether autumn temperature dependent sowing date is calculated
    ifsdautumn: bool,
    -- upper temperature limit for autumn sowing
    tempautumn: f64,
    -- lower temperature limit for spring sowing
    tempspring: f64,
    -- default length of growing period
    lgp_def: i64,
    -- upper minimum temperature limit for crop sowing
    maxtemp_sowing: f64,
    -- default sowing date in the northern hemisphere (julian day)
    sdatenh: i64,
    -- default sowing date in the southern hemisphere
    sdatesh: i64,
    -- whether sowing date adjusting equation is used
    sd_adjust: bool,
    -- parameter 1 in sowing date adjusting equation
    sd_adjust_par1: f64,
    -- parameter 2 in sowing date adjusting equation
    sd_adjust_par2: f64,
    -- parameter 3 in sowing date adjusting equation
    sd_adjust_par3: f64,
    -- latest date for harvesting in the northern hemisphere
    hlimitdatenh: i64,
    -- latest date for harvesting in the southern hemisphere
    hlimitdatesh: i64,
    -- default base temperature (°C) for heat unit (hu) calculation
    tb: f64,
    -- temperature under which vernalisation is possible (°C)
    trg: f64,
    -- default number of vernalising days required
    pvd: i64,
      -- sensitivity to the photoperiod effect [0-1]
    psens: f64,
    -- basal photoperiod (h) (pb<ps for longer days plants)
    pb: f64,
    -- lag in days after sowing before vernalization starts
    vern_lag: i64,
    -- saturating photoperiod (h) (ps<pb for shorter days plants)
    ps: f64,
    -- default potential heat units required for crop maturity (degree-days)
    phu: f64,
    -- whether quadratic equation used for calculating potential heat units (Bondeau method)
    phu_calc_quad: bool,
    -- whether linear equation used for calculating potential heat units (Bondeau method)
    phu_calc_lin: bool,
    -- minimum potential heat units required for crop maturity (Bondeau method) (degree-days)
    phu_min: f64,
    -- maximum potential heat units required for crop maturity (Bondeau method) (degree-days)
    phu_max: f64,
    -- reduction factor of potential heat units in spring crops (Bondeau method) (degree-days)
    phu_red_spring_sow: f64,
    -- number of days of phu decrease in the linear phu equation (Bondeau method)
    ndays_ramp_phu: f64,
    -- i64ercept for the linear phu equation (Bondeau method)
    phu_i64erc: f64,
    -- fraction of growing season (phu) at which senescence starts [0-1]
    fphusen: f64,
    -- type of senescence curve (see Bondeau et al. 2007)
    shapesenescencenorm: bool,
    -- fraction of maximal LAI still present at harvest [0-1]
    flaimaxharvest: f64,
    -- default maximum LAI (only used for i64ercrop grass in the case where no pasture is present in any stand)
    laimax: f64,
    -- whether harvestable organs are above ground
    aboveground_ho: bool,
    -- optimum harvest index
    hiopt: f64,
    -- minimum harvest index
    himin: f64,
    -- initial fraction of growing season's npp allocated to roots
    frootstart: f64,
    -- final fraction of growing season's npp allocated to roots
    frootend: f64,
    -- autumn/spring sowing of pft:s with tempautumn = 1
    forceautumnsowing: i64,  --0 = NOFORCING,  1 = AUTUMNSOWING, 2 = SPRINGSOWING
    -- N limited version of pft
    nlim: bool

  }

  let avg_cton (min: f64, max: f64) : f64 =
    2.0 / (1.0 / min + 1.0 / max)
  -- MEMBER FUNCTIONS

  -- Constructor (initialises array gdd0)
  let Pft() : Pft =
    let pft : Pft = {
      --std::fill_n(gdd0, Date::MAX_YEAR_LENGTH + 1, -1.0) -- value<0 signifies "unknown" see function phenology()
      gdd0 = replicate Date_MAX_YEAR_LENGTH_plusone (-1.0),
      nlim = false,

      root_beta = 0.0,

      drought_tolerance = 0.0, -- Default, means that the PFT will never be limited by drought.
      res_outtake = 0.0,
      harv_eff = 0.0,
      harv_eff_ic = 0.0,
      turnover_harv_prod = 1.0,  -- default 1 year turnover time

      isi64ercropgrass = false,
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
      i64c = nan,
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
      landcover = -1,
      leaflong = nan,
      leafphysiognomy = -1,
      lifeform = -1,

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
      pathway = -1,

      phengdd5ramp = nan,
      phenology = -1,
      phu_calc_lin = false,
      phu_calc_quad = false,

      phu_i64erc = nan,
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

    } in pft

  -- Calculates SLA given leaf longevity
  let initsla(pft: Pft) : Pft =

    -- SLA has to be supplied in the insfile for crops with N limitation
    if (!(pft.phenology == CROPGREEN && pft.nlim)) then
      -- Reich et al 1992, Table 1 (includes conversion x2.0 from m2/kg_dry_weight to
      -- m2/kgC)
      pft with sla =
      if (pft.leafphysiognomy == BROADLEAF) then
        0.2 * pow(10.0, 2.41 - 0.38 * f64.log10(12.0 * pft.leaflong))
      else if (pft.leafphysiognomy == NEEDLELEAF) then
        0.2 * pow(10.0, 2.29 - 0.4 * f64.log10(12.0 * pft.leaflong))
      else pft.sla
    else pft

  -- Calculates minimum leaf C:N ratio given leaf longevity
  let init_cton_min(pft: Pft) : Pft =
    -- cton_leaf_min has to be supplied in the insfile for crops with N limitation
    if (!(pft.phenology == CROPGREEN && pft.nlim)) then
      -- Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
      -- kgN/kgC)
      pft with cton_leaf_min =
      if (pft.leafphysiognomy == BROADLEAF) then
        500.0 / pow(10.0, 1.75 - 0.33 * f64.log10(12.0 * pft.leaflong))
      else if (pft.leafphysiognomy == NEEDLELEAF) then
        500.0 / pow(10.0, 1.52 - 0.26 * f64.log10(12.0 * pft.leaflong))
      else pft.cton_leaf_min
    else pft

  let init_cton_limits(pft: Pft) : Pft =
    -- Fraction between min and max C:N ratio White et al. 2000
    let frac_mi64omax = if (pft.phenology == CROPGREEN && pft.nlim) then 5.0 else 2.78  -- Use value also without nlim ?

    -- Fraction between leaf and root C:N ratio
    let frac_leaftoroot = 1.16 -- Friend et al. 1997

    -- Fraction between leaf and sap wood C:N ratio
    let frac_leaftosap = 6.9   -- Friend et al. 1997

    -- Max leaf C:N ratio
    let pft = pft with cton_leaf_max = pft.cton_leaf_min * frac_mi64omax

    -- Average leaf C:N ratio
    let pft = pft with cton_leaf_avr = avg_cton(pft.cton_leaf_min, pft.cton_leaf_max)

    -- Tighter C:N ratio range for roots and sapwood: picked out thin air
    let frac_maxtomin = 0.9

    -- Maximum fine root C:N ratio
    let pft = pft with cton_root_max = pft.cton_leaf_max * frac_leaftoroot

    let cton_root_min = pft.cton_root_max * frac_maxtomin

    -- Average fine root C:N ratio
    let pft = pft with cton_root_avr = avg_cton(cton_root_min, pft.cton_root_max)

    -- Maximum sap C:N ratio
    let pft = pft with cton_sap_max  = pft.cton_leaf_max * frac_leaftosap

    let cton_sap_min = pft.cton_sap_max * frac_maxtomin

    -- Average sap C:N ratio
    let pft = pft with cton_sap_avr  = avg_cton(cton_sap_min, pft.cton_sap_max)

    let pft = pft with respcoeff =
    if (pft.lifeform == GRASS || pft.lifeform == MOSS) then
      pft.respcoeff / 2.0 * pft.cton_root / (pft.cton_root_avr + cton_root_min)
    else
      pft.respcoeff / pft.cton_root / (pft.cton_root_avr + cton_root_min) +
                       pft.cton_sap  / (pft.cton_sap_avr  + cton_sap_min)
    let pft = pft with cton_stem_max = 1.0/(2.0*0.0034) --Maize params
    in pft with cton_stem_avr = 1.0/(2.0*0.0068)

  -- Calculates coefficient to compensate for different vertical distribution of fine root on nitrogen uptake
  let init_nupscoeff(pft: Pft) : Pft =
    -- Fraction fine root in upper soil layer should have higher possibility for mineralized nitrogen uptake
    -- Soil nitrogen profile is considered to have a exponential decline (Franzluebbers et al. 2009) giving
    -- an approximate advantage of 2 of having more roots in the upper soil layer
    let upper_adv = 2.0

    -- Simple solution until we get C and N in all soil layers.
    let (_, rootdist_upper, rootdist_lower) =
    loop (sl, rootdist_upper, rootdist_lower) = (0, 0.0, 0.0)
    while (sl < NSOILLAYER) do
      if (sl < NSOILLAYER_UPPER) then
        (sl+1, rootdist_upper+pft.rootdist[sl], rootdist_lower) else
        (sl+1, rootdist_upper, rootdist_lower+pft.rootdist[sl])
    in pft with nupscoeff = rootdist_upper * upper_adv + rootdist_lower


  -- Initialises sapling/regen characteristics in population mode following LPJF formulation
  let initregen(pft: Pft) : Pft =

    -- see function allometry in growth module.

    -- Note: primary PFT parameters, including SLA, must be set before this
    --       function is called

    let REGENLAI_TREE = 1.5
    let REGENLAI_GRASS = 0.001
    let SAPLINGHW = 0.2

    in if (pft.lifeform == TREE) then
      -- Tree sapling characteristics
      let pft = pft with regen_cmass_leaf =
        pow(REGENLAI_TREE * pft.k_allom1 * pow(1.0 + SAPLINGHW, pft.k_rp) *
          pow(4.0 * pft.sla / PI / pft.k_latosa, pft.k_rp * 0.5) / pft.sla, 2.0 / (2.0 - pft.k_rp))
      let pft = pft with regen_cmass_leaf =
        pft.wooddens * pft.k_allom2 * pow((1.0 + SAPLINGHW) *
        f64.sqrt(4.0 * pft.regen_cmass_leaf * pft.sla / PI / pft.k_latosa), pft.k_allom3) *
        pft.regen_cmass_leaf * pft.sla / pft.k_latosa
      in pft with regen_cmass_heart = SAPLINGHW * pft.regen_cmass_sap

    else if (pft.lifeform == GRASS || pft.lifeform==MOSS) then

      -- Grass regeneration characteristics
      pft with regen_cmass_leaf = REGENLAI_GRASS / pft.sla

    else
      pft with regen_cmass_root = 1.0 / pft.ltor_max * pft.regen_cmass_leaf

  -- Inits root fractions in each soil layer through a shape parameter beta (see Jackson et al., 1996)
  let init_rootdist(pft: Pft) : Pft =

    let depth = Dz_soil * CM_PER_MM

    let rootdist = copy pft.rootdist with [0] = 1.0 - pow(pft.root_beta, depth) -- init first layer

    let tot = rootdist[0]

    let (_, _, rootdist, tot) =
    loop (i, depth, rootdist, tot) =
         (1, depth, rootdist, tot) while (i<NSOILLAYER) do
         let i2 = i+1
         let depth2 = depth + (Dz_soil * CM_PER_MM)
         let rootdist2 = rootdist with [i] =  1.0 - pow(pft.root_beta, depth2) - (1.0 - pow(pft.root_beta, depth2 - Dz_soil * CM_PER_MM))
         let tot2 = tot + rootdist2[i]
         in (i2, depth2, rootdist2, tot2)

    -- Calibrated the root_beta for each PFT to match rootdist_upper from 'old' (pre LPJG 4.1) parameterisation.
    -- Sometimes the rootdist goes below our maximum soildepth. When that happens, put the residual fraction in lowest soil layer
    let rootdist[NSOILLAYER-1] = rootdist[NSOILLAYER-1] + 1.0 - tot
    in pft with rootdist = rootdist

  let ismoss(pft: Pft) : bool = pft.lifeform == MOSS
  let isgrass(pft: Pft) : bool = pft.lifeform == GRASS
  let istree(pft: Pft) : bool = pft.lifeform == TREE
  let iswetlandspecies(pft: Pft) : bool = (ismoss(pft) || pft.has_aerenchyma)
--}
