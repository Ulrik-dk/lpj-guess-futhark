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
let Climate(latitude : real) : Climate = {
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

let Individual(i : int, p : Pft, stand: Stand) : Individual = {
  id = i,
  pft = p,
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

  --Stand& stand = vegetation.patch.stand,

  -- there is a case where it is not initialized, but we cant have that
  cropindiv = (let newone = cropindiv_struct()
              in if (stand.pftid == p.id)
              then newone with isprimarycrop = true
              else if (stand.hasgrassintercrop && p.isintercropgrass)
              then newone with isintercropgrass = true
              else newone),

  --uninitialized:
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

let Patch(i : int, s : Stand, st : Soiltype) : Patch =
  let (_, patchpfts, individuals) = unzip3 <|
    map (\i ->
        let p = Pft()
        let patchp = Patchpft(i, p)
        let indv = Individual(i, p, s)
        in (p, patchp, indv)
      ) <| iota npft
  in
  {
  id = i,
  --stand = s,

  pfts = patchpfts,

  vegetation = individuals,

  soil = Soil(st),
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


let Standpft(i : int, p : Pft) : Standpft = {
  idx = i,
  pft = p,
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
let Stand(i : int,
          --Gridcell* gc,
          st : Soiltype,
          landcover : landcovertype,
          npatch_l : int,
          date : Date) : Stand =
  let num_patches =
    if (landcover == FOREST || landcover == NATURAL || (disturb_pasture && landcover == PASTURE))  -- TODO: obey this comment! \/ make it global
      then npatch -- use the global variable npatch for stands with stochastic events
    else if npatch_l > 0 then npatch_l -- use patch number provided by calling function
    else 1
  let st = Soiltype()
  let patches = map (\j -> Patch(j, i, st) ) <| iota num_patches 
  let p = Pft()
  in
  {id = i
--,gridcell=gc
  ,soiltype=st
  ,landcover=landcover
  ,original=landcover
  ,frac=1.0
  ,pft = replicate npft p
  ,data=patches
  first_year = date.year,
  clone_year = (-1),
  transfer_area_st = replicate nst realzero,
  seed = 12345678,
  stid = 0,
  pftid = (-1),
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
  origin=nan
}
