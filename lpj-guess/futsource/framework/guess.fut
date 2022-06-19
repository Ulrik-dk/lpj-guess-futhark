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
open import "guess_datatypes"




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

let MassBalance() : MassBalance = {
  start_year = nyear_spinup,
  ccont = 0.0,
  ccont_zero = 0.0,
  ccont_zero_scaled = 0.0,
  cflux = 0.0,
  cflux_zero = 0.0,
  ncont = 0.0,
  ncont_zero = 0.0,
  ncont_zero_scaled = 0.0,
  nflux = 0.0,
  nflux_zero = 0.0
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


let PhotosynthesisStresses() : PhotosynthesisStresses =
  {
    ifnlimvmax = false,
    moss_ps_limit = 1.0,
    graminoid_ps_limit = 1.0,
    inund_stress = 1.0
  }


--- FIXME: This code was taken from weathergen.cpp!
let WeatherGenState() : WeatherGenState =
  {
    q = replicate QSIZ 0,
    carry = 362,
    xcng  = 1236789,
    xs    = 521288629, --default seed
    indx  = QSIZ+1,
    have  = false,
    gamma_vals = replicate 2 realzero,
    pday = replicate 2 false,
    resid = replicate 4 realzero
  }


  --- constructor function: initialises gridcell member
  -- takes Gridcell& gc
let Climate(latitude : real, climate_id : int, gridcell_id : int) : Climate = {
  gridcell_id = gridcell_id,
  climate_id = climate_id,
  aprec = 0.0,
  aprec_lastyear = 0.0,
  mtemp20 = replicate 12 realzero,
  mprec20 = replicate 12 realzero,
  mpet20 = replicate 12 realzero,
  mpet_year = replicate 12 realzero,
  mprec_pet20 = replicate 12 realzero,

  mtemp_20 = replicate 20 (replicate 12 realzero),
  mprec_20 = replicate 20 (replicate 12 realzero),
  mpet_20 = replicate 20 (replicate 12 realzero),
  mprec_pet_20 = replicate 20 (replicate 12 realzero),

  mprec_petmin_20 = replicate 20 realzero,
  mprec_petmax_20 = replicate 20 realzero,

  mprec_petmin20=0.0,
  mprec_petmax20=0.0,

  seasonality=SEASONALITY_NO,
  seasonality_lastyear=SEASONALITY_NO,
  prec_seasonality=DRY,
  prec_seasonality_lastyear=DRY,
  prec_range=DRY,
  prec_range_lastyear=DRY,
  temp_seasonality=COLD,
  temp_seasonality_lastyear=COLD,
  biseasonal=false,

  eet=0.0,

  gdd0 = 0.0,
  agdd0 = 0.0,

  --- Initialises certain member variables
  --- Should be called before Climate object is applied to a new grid cell */
  --void initdrivers(double latitude) {
  -- Futhark: Fuzed "initializer" with constructor
  -- FIXME this might not be correct in practice, but its quick

  mtemp_min_20 = replicate 20 realzero,
  mtemp_max_20 = replicate 20 realzero,

  mtemp_min20 = 0.0,
  mtemp_max20 = 0.0,
  mtemp = 0.0,
  maxtemp = 0.0,
  gdd5 = 0.0,
  chilldays = 0,
  ifsensechill = true,
  atemp_mean = 0.0,

  lat = latitude,

  doneday = replicate Date_MAX_YEAR_LENGTH false,
  sinelat = sin(latitude * DEGTORAD),
  cosinelat = cos(latitude * DEGTORAD),

  testday_temp = if (latitude >= 0) then 180 else 364,
  testday_prec = if (latitude >= 0) then 364 else 180,
  coldestday = if (latitude >= 0) then COLDEST_DAY_NHEMISPHERE else COLDEST_DAY_SHEMISPHERE,
  adjustlat = if (latitude >= 0) then 0 else 181,

  -- UNINITIALIZED
  temps = replicate Date_subdaily 0.0,
  insols = replicate Date_subdaily 0.0,
  pars = replicate Date_subdaily 0.0,
  rads = replicate Date_subdaily 0.0,
  gtemps = replicate Date_subdaily 0.0,
  agdd5=nan,
  avg_annual_rainfall=nan,
  co2=nan,
  cur_rainfall=nan,
  daylength=nan,
  daylength_save=replicate Date_MAX_YEAR_LENGTH nan,
  days_since_last_rainfall=nan,
  dprec_10=replicate 10 nan,
  dtr=nan,
  gtemp=nan,
  hh=replicate Date_MAX_YEAR_LENGTH nan,
  insol=nan,
  instype=intnan,
  kbdi=nan,
  last_rainfall=nan,
  mcarthur_forest_fire_index=nan,
  months_ffdi=replicate 30 nan,
  mtemp_max=nan,
  mtemp_min=nan,
  par=nan,
  prec=nan,
  qo=replicate Date_MAX_YEAR_LENGTH nan,
  rad=nan,
  relhum=nan,
  sinehh=replicate Date_MAX_YEAR_LENGTH nan,
  sprec_2=replicate 2 nan,
  temp=nan,
  tmax=nan,
  tmin=nan,
  u=replicate Date_MAX_YEAR_LENGTH nan,
  u10=nan,
  v=replicate Date_MAX_YEAR_LENGTH nan,
  var_prec=nan,
  var_temp=nan,
  weathergenstate=WeatherGenState()
}


let Date() : Date =
  let data = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  let scanned = scan (+) 0 data
  let middaymonth = map (\i -> scanned[i] + data[i] / 2) <| iota 12 --TODO FIXME: chech that this is right, it probably isnt
  in
  {
  ndaymonth = copy data,
  middaymonth = middaymonth,
  subdaily = 1,
  first_calendar_year = 0,
  ---
  MAX_YEAR_LENGTH = intnan,
  day = intnan,
  dayofmonth = intnan,
  islastday = false,
  islastmonth = false,
  islastyear = false,
  ismidday = false,
  month = intnan,
  nyear = intnan,
  year = intnan
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
    -- For the fluxes only stored as totals for the whole patch--
    monthly_fluxes_patch = replicate 12 (replicate NPERPATCHFLUXTYPES realzero),

    -- Stores one flux value per month and flux type
    --- For the fluxes stored per pft for annual values--
    monthly_fluxes_pft = replicate 12 (replicate NPERPFTFLUXTYPES realzero),

    -- Stores one flux value per day and flux type
    daily_fluxes_patch = replicate 365 (replicate NPERPATCHFLUXTYPES realzero),

    -- Stores one flux value per day and flux type
    daily_fluxes_pft = replicate 365 (replicate NPERPFTFLUXTYPES realzero)
  }


let CropRotation() : CropRotation = {
  ncrops = 0,
  firstrotyear = 0
}




let ManagementType(managementtype_id : int) : ManagementType = {
  managementtype_id = managementtype_id,
  planting_system = planting_system_NONE,
  harvest_system = harvest_system_NONE,
--  pftname = "",
--  selection = "",
  nyears = 1.0,
  hydrology = RAINFED,
--    firr = 0.0,
  sdate = -1,
  hdate = -1,
  nfert = (-1.0),
  fallow = false,
  multicrop = false
}


let StandType(standtype_id : int) : StandType = {
  standtype_id = standtype_id,
  intercrop = NOINTERCROP,
  naturalveg = naturalvegNONE,
  restrictpfts = false,
  reestab = reestabALL,
  firstmanageyear = 100000,
  management = ManagementType(0),
  landcover = NATURAL, --TODO FIXME: READ FROM FILE
  rotation = CropRotation()
}

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

let istruecrop_or_intercropgrass(indiv: Individual, pft : Pft) : bool =
  (pft.landcover==CROPLAND && (pft.phenology==CROPGREEN || indiv.cropindiv.isintercropgrass))


let cmass_leaf_today(this : Individual, pft : Pft, patchpft: Patchpft) : real =
  if (istruecrop_or_intercropgrass(this, pft)) then
    if patchpft.cropphen.growingseason then this.cropindiv.grs_cmass_leaf else 0
  else
    this.cmass_leaf * this.phen



let lai_nitrogen_today(this : Individual, patchpft : Patchpft, pft: Pft) : real =
	if (pft.phenology==CROPGREEN) then
		if (patchpft.cropphen.growingseason && cmass_leaf_today(this, pft, patchpft) > 0.0) then
			let k = 0.5
			let ktn = 0.52*k + 0.01 -- Yin et al 2003
			let nb = 1/(pft.cton_leaf_max*pft.sla)
			in (1/ktn) * log(1+ktn*this.nmass_leaf/nb)
		else 0.0
  else 1.0


--- Gets the individual's daily cmass_root value
let cmass_root_today(this : Individual, pft: Pft, patchpft: Patchpft) : real =
  if (istruecrop_or_intercropgrass(this, pft)) then
    if patchpft.cropphen.growingseason then this.cropindiv.grs_cmass_root else 0
  else
    this.cmass_root * this.phen

let Individual_report_flux_PerPFTFluxType(this: *Fluxes, alive: bool, istruecrop_or_intercropgrass: bool, flux_type : PerPFTFluxType, value : real, date : Date, pft_id : int) =
  if (alive || istruecrop_or_intercropgrass) then
    report_flux_PerPFTFluxType(copy this, flux_type, pft_id, value, date)
  else this

let Individual_report_flux_PerPatchFluxType(this: *Fluxes, alive: bool, istruecrop_or_intercropgrass: bool, flux_type : PerPFTFluxType, value : real, date : Date) =
  if (alive || istruecrop_or_intercropgrass) then
    report_flux_PerPatchFluxType(copy this, flux_type, value, date)
  else this

--let get_annual_flux_PerPFTFluxType_pft_id(this: Fluxes, flux_type: PerPFTFluxType, pft_id: int) : real =
--  this.annual_fluxes_per_pft[pft_id][flux_type]

let get_annual_flux_PerPFTFluxType(this: Fluxes, flux_type: PerPFTFluxType) : real =
  reduce (+) realzero <| (transpose this.annual_fluxes_per_pft)[flux_type]

let get_annual_flux_PerPatchFluxType(this: Fluxes, flux_type: PerPatchFluxType) : real =
  reduce (+) realzero <| (transpose this.annual_fluxes_per_pft)[flux_type]

let avg_cton (min: real, max: real) : real =
  2.0 / (1.0 / min + 1.0 / max)



-- MEMBER FUNCTIONS
-- Constructor (initialises array gdd0)
let Pft(pft_id : int) : Pft = {
  pft_id = pft_id,
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



let cropindiv_struct() : cropindiv_struct = {
  cmass_ho=0.0,
  cmass_agpool=0.0,
  cmass_stem = 0.0,
  cmass_leaf_max=0.0,
  cmass_leaf_sen=0.0,
  yield=0.0,
  yield_harvest=replicate 2 realzero,
  dcmass_leaf=0.0,
  dcmass_root=0.0,
  dcmass_plant=0.0,
  dcmass_ho=0.0,
  dcmass_agpool=0.0,
  grs_cmass_leaf=0.0,
  grs_cmass_root=0.0,
  grs_cmass_plant=0.0,
  grs_cmass_ho=0.0,
  grs_cmass_agpool=0.0,
  grs_cmass_stem = 0.0,
  grs_cmass_dead_leaf = 0.0,
  grs_cmass_leaf_luc=0.0,
  grs_cmass_root_luc=0.0,
  grs_cmass_ho_luc=0.0,
  grs_cmass_agpool_luc=0.0,
  grs_cmass_dead_leaf_luc = 0.0,
  grs_cmass_stem_luc = 0.0,
  nmass_ho=0.0,
  nmass_agpool=0.0,
  nmass_dead_leaf = 0.0,
  ycmass_leaf=0.0,
  ycmass_root=0.0,
  ycmass_plant=0.0,
  ycmass_ho=0.0,
  ycmass_agpool=0.0,
  ycmass_stem = 0.0,
  ycmass_dead_leaf = 0.0,
  harv_cmass_leaf=0.0,
  harv_cmass_root=0.0,
  harv_cmass_ho=0.0,
  harv_yield=0.0,
  harv_cmass_agpool=0.0,
  harv_cmass_stem = 0.0,
  cmass_ho_harvest = replicate 2 realzero,
  --Nitrogen
  dnmass_leaf=0.0,
  dnmass_root=0.0,
  dnmass_ho=0.0,
  dnmass_agpool=0.0,
  ynmass_leaf=0.0,
  ynmass_root=0.0,
  ynmass_ho=0.0,
  ynmass_agpool=0.0,
  ynmass_dead_leaf = 0.0,
  harv_nmass_leaf=0.0,
  harv_nmass_root=0.0,
  harv_nmass_ho=0.0,
  harv_nmass_agpool=0.0,
  nmass_dead_leaf_luc = 0.0,
  nmass_ho_harvest = replicate 2 realzero,

  isprimarycrop=false,
  isprimarycovegetation=false,
--    issecondarycrop=false,
  isintercropgrass=false,

  -- UNDEFINED FIELDS:
  dcmass_stem=nan,
  harv_cmass_plant=nan,
  nmass_agpool_luc=nan,
  nmass_ho_luc=nan
}


-- Returns true if stand is true high-latitude peatland stand, as opposed to a wetland < PEATLAND_WETLAND_LATITUDE_LIMIT N
let is_highlatitude_peatland_stand(this: Stand, gridcell : Gridcell) : bool =
  let lat : real = gridcell.lat
  in this.landcover==PEATLAND && lat >= PEATLAND_WETLAND_LATITUDE_LIMIT

-- Returns true if stand is wetland stand, as opposed to a peatland >= PEATLAND_WETLAND_LATITUDE_LIMIT N
let is_true_wetland_stand(this: Stand, gridcell : Gridcell) : bool =
  let lat : real = gridcell.lat
  in this.landcover==PEATLAND && lat < PEATLAND_WETLAND_LATITUDE_LIMIT


let Individual(individual_id : int
              ,gridcell_id : int
              ,patch_id : int
              ,stand_id : int
              ,pft_id : int
              ,stand_pftid : int
              ,stand_hasgrassintercrop : bool
              ,pft_isintercropgrass : bool
            --,stand: Stand [num_patches]
              ) : Individual =

  let ps_res = PhotosynthesisResult() in
  {
  individual_id = individual_id,
  pft_id = pft_id,
  gridcell_id = gridcell_id,
  patch_id = patch_id,
  stand_id = stand_id,
  --vegetation = v,
  anpp              = 0.0,
  fpc               = 0.0,
  fpc_daily      = 0.0,
  densindiv         = 0.0,
  cmass_leaf        = 0.0,
  cmass_root        = 0.0,
  cmass_sap         = 0.0,
  cmass_heart       = 0.0,
  cmass_debt        = 0.0,
  cmass_leaf_post_turnover      = 0.0,
  cmass_root_post_turnover      = 0.0,
  cmass_tot_luc     = 0.0,
  phen              = 0.0,
  aphen             = 0.0,
  deltafpc          = 0.0,

  nmass_leaf        = 0.0,
  nmass_root        = 0.0,
  nmass_sap         = 0.0,
  nmass_heart       = 0.0,
  cton_leaf_aopt    = 0.0,
  cton_leaf_aavr    = 0.0,
  cton_status       = 0.0,
  cmass_veg         = 0.0,
  nmass_veg         = 0.0,
  nmass_tot_luc     = 0.0,

  nactive           = 0.0,
  nextin            = 1.0,
  nstore_longterm   = 0.0,
  nstore_labile     = 0.0,
  ndemand           = 0.0,
  fnuptake          = 1.0,
  anuptake          = 0.0,
  max_n_storage     = 0.0,
  scale_n_storage   = 0.0,

  leafndemand       = 0.0,
  rootndemand       = 0.0,
  sapndemand        = 0.0,
  storendemand      = 0.0,
  leaffndemand      = 0.0,
  rootfndemand      = 0.0,
  sapfndemand       = 0.0,
  storefndemand     = 0.0,
  leafndemand_store = 0.0,
  rootndemand_store = 0.0,

  nstress           = false,

  -- additional initialisation
  age               = 0.0,
  fpar              = 0.0,
  aphen_raingreen   = 0i64,
  intercep          = 0.0,
  phen_mean         = 0.0,
  wstress           = false,
  lai               = 0.0,
  lai_layer         = 0.0,
  lai_indiv         = 0.0,
  lai_daily         = 0.0,
  lai_indiv_daily   = 0.0,
  alive             = false,

  mlai = replicate 12 realzero,
  mlai_max = replicate 12 realzero,

  -- bvoc
  iso               = 0.0,
  fvocseas          = 1.0,
  mon = replicate NMTCOMPOUNDS realzero,
  monstor = replicate NMTCOMPOUNDS realzero,
  dnpp              = 0.0,
  last_turnover_day = (-1i64),
  gpterms = replicate Date_subdaily realzero,

  --Stand& stand = vegetation.patch.stand,

  -- there is a case where it is not initialized, but we cant have that
  cropindiv = (let newone = cropindiv_struct()
              in if (stand_pftid == pft_id)
              then newone with isprimarycrop = true
              else if (stand_hasgrassintercrop && pft_isintercropgrass)
              then newone with isintercropgrass = true
              else newone),

  --uninitialized:
  phots=replicate Date_subdaily ps_res,
  aaet=nan,
  aet=nan,
  avmaxnlim=nan,
  boleht=nan,
  crownarea=nan,
  daily_cmass_leafloss=nan,
  daily_cmass_rootloss=nan,
  daily_nmass_leafloss=nan,
  daily_nmass_rootloss=nan,
  fpar_leafon=nan,
  gpterm=nan,
  height=nan,
  hondemand=nan,
  lai_leafon_layer=nan,
  ltor=nan,
  nday_leafon=intnan,
  nmass_heart_luc=nan,
  nmass_leaf_luc=nan,
  nmass_root_luc=nan,
  nmass_sap_luc=nan,
  nstore_labile_luc=nan,
  nstore_longterm_luc=nan,
  photosynthesis_result=PhotosynthesisResult()
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


-- Constructor: initialises id, pft and data members
let Patchpft(gridcell_id: int, patch_id: int, stand_id: int, patchpft_id: int, pft_id: int) : Patchpft = {

  gridcell_id = gridcell_id,
  stand_id = stand_id,
  patch_id = patch_id,

  patchpft_id = patchpft_id,
  pft_id = pft_id,

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

let diurnal(this : Date) : bool =
  this.subdaily > 1


--- Constructs beginning of the day period (the only one in daily mode)
let Day(date : Date) : Day = {
  isstart = true,
  isend = !diurnal(date),
  period = 0
}

--- Gets the growingseason status for crop individual. Non-crop individuals always return true.
let growingseason(patchpft: Patchpft) : bool =
	patchpft.cropphen.growingseason -- FIXME is this right?

--- Advances to the next sub-daily period
let next(this : Day, date: Date) : Day = {
  period = this.period+1,
  isstart = false,
  isend = this.period+1 == date.subdaily - 1
}

--- Gets the individual's daily fpc value
let fpc_today(this: Individual, pft: Pft, ppft: Patchpft) : real =
  if (pft.phenology == CROPGREEN) then
    if ppft.cropphen.growingseason then this.fpc_daily else 0.0
  else
    this.fpc * this.phen

--- Gets the individual's daily lai_indiv value
let lai_indiv_today(this: Individual, pft: Pft, ppft: Patchpft) : real =
	if (pft.phenology == CROPGREEN) then
		if ppft.cropphen.growingseason then this.lai_indiv_daily else 0
	else
		this.lai_indiv * this.phen
--
let is_true_crop_stand(stand: Stand, pft: Pft) : bool =
  stand.landcover==CROPLAND && pft.phenology==CROPGREEN--	-- OK also for fallow (pftid always cropgreen)

let cton_leaf(this : Individual, use_phen : bool, stand : Stand, standpft: Standpft, pfts : [npft]Pft, patchpft: Patchpft) : real =
	if (is_true_crop_stand(stand, pfts[standpft.pft_id]) && !negligible(cmass_leaf_today(this, pfts[standpft.pft_id], patchpft)) && !negligible(this.nmass_leaf))
  then cmass_leaf_today(this, pfts[standpft.pft_id], patchpft) / this.nmass_leaf
	else
    if (!is_true_crop_stand(stand, pfts[standpft.pft_id]) && !negligible(this.cmass_leaf) && !negligible(this.nmass_leaf)) then
  		if (use_phen) then
  			if (!negligible(this.phen)) then
  				cmass_leaf_today(this, pfts[standpft.pft_id], patchpft) / this.nmass_leaf
  			else pfts[standpft.pft_id].cton_leaf_avr
  		else this.cmass_leaf / this.nmass_leaf
  	else pfts[standpft.pft_id].cton_leaf_max


let cton_sap(this : Individual, pft : Pft) : real =
	if (pft.lifeform == TREE) then
		if (!negligible(this.cmass_sap) && !negligible(this.nmass_sap)) then
			max(pft.cton_sap_avr * pft.cton_leaf_min / pft.cton_leaf_avr, this.cmass_sap / this.nmass_sap)
		else
			pft.cton_sap_max
	else
		1.0
--
let cton_root(this : Individual, use_phen : bool, pft : Pft, patchpft: Patchpft) : real =
	if (!negligible(this.cmass_root) && !negligible(this.nmass_root)) then
		if (use_phen) then
			if (!negligible(cmass_root_today(this, pft, patchpft))) then
				max(pft.cton_root_avr * pft.cton_leaf_min / pft.cton_leaf_avr, cmass_root_today(this, pft, patchpft) / this.nmass_root)
			else
				pft.cton_root_avr
		else
			max(pft.cton_root_avr * pft.cton_leaf_min / pft.cton_leaf_avr, this.cmass_root / this.nmass_root)
	else
		pft.cton_root_max



--- Total storage of nitrogen
let nstore(individual : Individual) : real =
	individual.nstore_longterm + individual.nstore_labile

--- Total carbon wood biomass
let cmass_wood(individual : Individual) : real =
	individual.cmass_sap + individual.cmass_heart - individual.cmass_debt

--- Total nitrogen wood biomass
let nmass_wood(individual : Individual) : real =
	individual.nmass_sap + individual.nmass_heart

let ndemand_storage(this: Individual, cton_leaf_opt : real, stand : Stand, pfts: [npft]Pft, standpft: Standpft, patchpft : Patchpft) : real =
	if (is_true_crop_stand(stand, pfts[standpft.pft_id]) && ifnlim) then	-- only CROPGREEN, only ifnlim ?
		-- analogous with root demand
		max(0.0, this.cropindiv.grs_cmass_stem / (cton_leaf_opt * pfts[standpft.pft_id].cton_stem_avr / pfts[standpft.pft_id].cton_leaf_avr) - this.cropindiv.nmass_agpool)
	else
		max(0.0, min(this.anpp * this.scale_n_storage / cton_leaf(this, true, stand, standpft, pfts, patchpft), this.max_n_storage) - nstore(this))

--- Gets the individual's daily lai value
let lai_today(this : Individual, pft: Pft, ppft: Patchpft) : real =
  if (pft.phenology == CROPGREEN) then
    if ppft.cropphen.growingseason then this.lai_daily else 0
  else
    this.lai * this.phen


let Patch(gridcell_id: int
         ,stand_id: int
         ,patch_id : int
         ,stand_pftid : int
         ,stand_hasgrassintercrop : bool
         ,st : Soiltype)
         : Patch =
  let (_, _, _) = unzip3 <| -- TODO FIXME return these for global the masterobject
    map (\i ->
        let pft_id = 0
        let pft = Pft(pft_id)
        let patchpft = Patchpft(gridcell_id, patch_id, stand_id, i, pft_id)
        let individual = Individual(i, gridcell_id, patch_id, stand_id, pft_id, stand_pftid, stand_hasgrassintercrop, pft.isintercropgrass)
        in (pft, patchpft, individual)
      ) <| iota npft
  in
{
  gridcell_id = gridcell_id,
  stand_id = stand_id,
  patch_id = patch_id,

  --pfts = patchpfts,

  --vegetation = individuals,

  --soil = Soil(st),
  fluxes = Fluxes(),

  age = 0,
  disturbed = false,
  managed = false,
  man_strength = 0.0,
  managed_this_year = false,
  plant_this_year = false,
  wdemand = 0.0,
  wdemand_leafon = 0.0,

  growingseasondays = 0,

  fireprob = 0.0,
  ndemand = 0.0,
  dnfert = 0.0,
  anfert = 0.0,
  nharv = 0,

  -- UNUSED? aaet_5 = replicate NYEARAAET realzero,

  avg_fbrlt = replicate N_YEAR_BIOMEAVG realzero,
  avg_fgrass = replicate N_YEAR_BIOMEAVG realzero,
  avg_fndlt = replicate N_YEAR_BIOMEAVG realzero,
  avg_fshrb = replicate N_YEAR_BIOMEAVG realzero,
  avg_ftot = replicate N_YEAR_BIOMEAVG realzero,

  -- FIXME unshared fields
  aaet = nan,
  abaserunoff = nan,
  adrainrunoff = nan,
  aevap = nan,
  aintercep = nan,
  apet = nan,
  arunoff = nan,
  asurfrunoff = nan,
  awetland_water_added = nan,
  eet_net_veg = nan,
  fire_line_intensity = nan,
  fpar_ff = nan,
  fpar_grass = nan,
  fpc_rescale = nan,
  fpc_total = nan,
  intercep = nan,
  irrigation_d = nan,
  irrigation_y = nan,
  is_litter_day = false,
  isharvestday = false,
  lcwd_to_atm = nan,
  leaf_to_atm = nan,
  leaf_to_lit = nan,
  lfwd_to_atm = nan,
  litf_to_atm = nan,
  maet = replicate 12 nan,
  mevap = replicate 12 nan,
  mintercep = replicate 12 nan,
  mpet = replicate 12 nan,
  mrunoff = replicate 12 nan,
  nday_growingseason = intnan,
  par_grass_mean = nan,
  wdemand_day = nan,
  wetland_water_added_today = nan,
  wood_to_atm = nan,
  wood_to_cwd = nan,
  wood_to_fwd = nan,
  wood_to_str = nan
  }


let Standpft(gridcell_id: int, stand_id: int, standpft_id: int, pft_id : int) : Standpft = {
  gridcell_id = gridcell_id,
  stand_id = stand_id,

  standpft_id = standpft_id,
  pft_id = pft_id,

  anetps_ff_max = 0.0,
  active = !run_landcover,
  plant = false,
  reestab = false,
  irrigated = false,
  sdate_force = -1,
  hdate_force = -1,
  -- uninitialized
  cmass_repr=nan,
  fpc_total=nan,
  photosynthesis_result = PhotosynthesisResult()
}


let npatch : int = 15 --from global.ins -- TODO FIXME move this somewhere reasonable

--------------------------------------------------------------------------------
-- Implementation of Stand member functions
--------------------------------------------------------------------------------

let Stand(stand_id : int,
          gridcell_id: int,
          soiltype_id : int,
          standtype_id : int,
          --witness : [num_patches](),
          --num_patches : int,
          landcover : landcovertype,
          npatch_l : int,
          date : Date) : ([]Patch, []Standpft, Stand) =
  let num_patches =
    if (landcover == FOREST || landcover == NATURAL || (disturb_pasture && landcover == PASTURE))
      then npatch -- use the global variable npatch for stands with stochastic events
    else if npatch_l > 0 then npatch_l -- use patch number provided by calling function
    else 1
  let st = Soiltype(soiltype_id)
  let pft_id = (-1)
  let stand_pftid = pft_id
  let patches = map (\patch_id -> Patch(gridcell_id, stand_id, patch_id, stand_pftid, false, st)) <| iota num_patches
  let standpfts = map (\standpft_id -> Standpft(gridcell_id, stand_id, standpft_id, standpft_id)) <| iota npft
  in (patches, standpfts, {
    pft_id = pft_id,
    num_patches = num_patches,
    standtype_id = standtype_id,
    gridcell_id = gridcell_id,
    stand_id = stand_id,
    soiltype_id=soiltype_id,
    landcover=landcover,
    original=landcover,
    frac=1.0,
    --standpft = standpfts,
    --data=patches,
    first_year = date.year,
    clone_year = (-1),
    transfer_area_st = replicate nst realzero,
    seed = 12345678,
    --stid = 0,
    current_rot = 0,
    ndays_inrotation = 0,
    infallow = false,
    isrotationday = false,
    isirrigated = false,
    hasgrassintercrop = false,
    gdd5_intercrop = 0.0,
    frac_old = 0.0,
    frac_temp = 0.0,
    protected_frac = 0.0,
    frac_change = 0.0,
    gross_frac_increase = 0.0,
    gross_frac_decrease = 0.0,
    cloned_fraction = 0.0,
    cloned = false,
    anpp = 0.0,
    cmass = 0.0,
    scale_LC_change = 1.0,
    -- UNINITIALIZED
    origin=intnan
  })


let Landcover() : Landcover = {

  updated = false,

  acflux_harvest_slow = 0.0,
  acflux_landuse_change = 0.0,
  anflux_harvest_slow = 0.0,
  anflux_landuse_change = 0.0,

  frac = replicate NLANDCOVERTYPES realzero,
  frac_old = replicate NLANDCOVERTYPES realzero,
  frac_change = replicate NLANDCOVERTYPES realzero,
  acflux_harvest_slow_lc = replicate NLANDCOVERTYPES realzero,
  acflux_landuse_change_lc = replicate NLANDCOVERTYPES realzero,
  anflux_harvest_slow_lc = replicate NLANDCOVERTYPES realzero,
  anflux_landuse_change_lc = replicate NLANDCOVERTYPES realzero,

  frac_transfer = replicate NLANDCOVERTYPES (replicate NLANDCOVERTYPES realzero),
  primary_frac_transfer = replicate NLANDCOVERTYPES (replicate NLANDCOVERTYPES realzero),

  expand_to_new_stand = map (\i -> (i == NATURAL || i == FOREST)) <| iota NLANDCOVERTYPES,

  pool_to_all_landcovers = replicate NLANDCOVERTYPES false, --from a donor landcover, alt.c
  pool_from_all_landcovers = replicate NLANDCOVERTYPES false --to a receptor landcover, alt.a
}

--/ Constructs a Gridcellpft object
--  \param i   The id for this object
--  \param p   A reference to the Pft for this Gridcellpft
--/
let Gridcellpft(gridcellpft_id: int, gridcell_id : int, pft_id : int) : Gridcellpft = {
  gridcellpft_id = gridcellpft_id,
  pft_id = pft_id,
  addtw = 0.0,
  Km = 0.0,

  autumnoccurred=false,
  springoccurred=false,
  vernstartoccurred=false,
  vernendoccurred=false,
  first_autumndate=(-1),
  first_autumndate20=(-1),
  last_springdate=(-1),
  last_springdate20=(-1),
  last_verndate=(-1),
  last_verndate20=(-1),
  first_autumndate_20 = replicate 20 (-1),
  last_springdate_20 = replicate 20 (-1),
  last_verndate_20 = replicate 20 (-1),
  sdate_default=(-1),
  sdate_force=(-1),
  hdate_force=(-1),
  Nfert_read=(-1),
  Nfert_man_read=(-1),
  sdatecalc_temp=(-1),
  sdatecalc_prec=(-1),
  hlimitdate_default=(-1),
  wintertype=false,
  swindow = replicate 2 (-1),
  sowing_restriction = false,
  --
  swindow_irr = replicate 2 (-1)
}

--/ Constructs a Gridcellst object
--  \param i   The id for this object
--  \param s   A reference to the StandType for this Gridcellst
--/
let Gridcellst(gridcellst_id : int, standtype_id : int) : Gridcellst = {
  gridcellst_id = gridcellst_id,
  standtype_id = standtype_id,
  frac = 1.0,
  frac_old = 0.0,
  frac_old_orig = 0.0,
  protected_frac = 0.0,
  frac_change = 0.0,
  gross_frac_increase = 0.0,
  gross_frac_decrease = 0.0,
  nstands = 0,
  nfert = (-1.0)
}
type Vegetation = [npft]Individual

let Gridcell() : Gridcell = {
  --TODO this constructor is wrong
	seed = 12345678,
  aNH4dep = nan,
  aNO3dep = nan,
  ann_max_fapar = nan,
  annual_burned_area = nan,
  burned_area = nan,
  burned_area_accumulated = nan,
  can_burn = intnan,
  climate_id = intnan,
  cur_max_fapar = nan,
  cur_nesterov = nan,
  dNH4dep = nan,
  dNO3dep = nan,
  gridcellpft_id = intnan,
  gridcellst_id = intnan,
  hyde31_pop_density = replicate 57 nan,
  k_tun_litter = nan,
  landcover_id = intnan,
  lat = nan,
  lon = nan,
  massbalance_id = intnan,
  max_nesterov = nan,
  monthly_burned_area = replicate 12 nan,
  monthly_fire_risk = replicate 12 nan,
  monthly_max_nesterov = replicate 12 nan,
  pop_density = nan,
  recent_max_fapar = replicate AVG_INTERVAL_FAPAR nan,
  simfire_biome = intnan,
  simfire_region = intnan,
  soiltype_id = intnan
}
