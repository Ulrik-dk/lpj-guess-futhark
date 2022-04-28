-- Holds static functional parameters for a plant functional type (PFT).
-- There should be one Pft object for each potentially occurring PFT. The same Pft object
-- may be referenced (via the pft member of the Individual object see below) by different
-- average individuals. Member functions are included for initialising SLA given leaf
-- longevity, and for initialising sapling/regen characteristics (required for
-- population mode).

let Date_MAX_YEAR_LENGTH_plusone = 367

module Pft = {
  type~ Pft  = {
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
    phengdd5ramp: double,
    -- water stress threshold for leaf abscission (range 0-1 raingreen PFTs)
    wscal_min: double,
    -- biochemical pathway for photosynthesis (C3 or C4)
    pathway: pathwaytype,
    -- approximate low temperature limit for photosynthesis (deg C)
    pstemp_min: double,
    -- approximate lower range of temperature optimum for photosynthesis (deg C)
    pstemp_low: double,
    -- approximate upper range of temperature optimum for photosynthesis (deg C)
    pstemp_high: double,
    -- maximum temperature limit for photosynthesis (deg C)
    pstemp_max: double,
    -- non-water-stressed ratio of intercellular to ambient CO2 partial pressure
    lambda_max: double,
    -- vegetation root profile in an array containing fraction of roots in each soil layer, [0=upper layer]
    rootdist: double[NSOILLAYER],
    -- shape parameter for initialisation of root distribtion
    root_beta: double,
    -- canopy conductance component not associated with photosynthesis (mm/s)
    gmin: double,
    -- maximum evapotranspiration rate (mm/day)
    emax: double,
    -- maintenance respiration coefficient (0-1)
    respcoeff: double,

    -- minimum leaf C:N mass ratio allowed when nitrogen demand is determined
    cton_leaf_min: double,
    -- maximum leaf C:N mass ratio  allowed when nitrogen demand is determined
    cton_leaf_max: double,
    -- average leaf C:N mass ratio (between min and max)
    cton_leaf_avr: double,
    -- average fine root C:N mass ratio (connected cton_leaf_avr)
    cton_root_avr: double,
    -- maximum fine root C:N mass ratio (used when mass is negligible)
    cton_root_max: double,
    -- average sapwood C:N mass ratio (connected cton_leaf_avr)
    cton_sap_avr: double,
    -- maximum sapwood C:N mass ratio (used when mass is negligible)
    cton_sap_max: double,
    -- reference fine root C:N mass ratio
    cton_root: double,
    -- reference sapwood C:N mass ratio
    cton_sap: double,
    -- Maximum nitrogen (NH4+ and NO3- seperatly) uptake per fine root [kgN kgC-1 day-1]
    nuptoroot: double,
    -- coefficient to compensate for vertical distribution of fine root on nitrogen uptake
    nupscoeff: double,
    -- fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen longterm storage scalar
    fnstorage: double,

    -- Michaelis-Menten kinetic parameters
    -- Half saturation concentration for N uptake [kgN l-1] (Rothstein 2000) */
    km_volume: double,

    -- fraction of NPP allocated to reproduction
    reprfrac: double,
    -- annual leaf turnover as a proportion of leaf C biomass
    turnover_leaf: double,
    -- annual fine root turnover as a proportion of fine root C biomass
    turnover_root: double,
    -- annual sapwood turnover as a proportion of sapwood C biomass
    turnover_sap: double,
    -- sapwood and heartwood density (kgC/m3)
    wooddens: double,
    -- maximum tree crown area (m2)
    crownarea_max: double,
    -- constant in allometry equations
    k_allom1: double,
    -- constant in allometry equations
    k_allom2: double,
    -- constant in allometry equations
    k_allom3: double,
    -- constant in allometry equations
    k_rp: double,
    -- tree leaf to sapwood area ratio
    k_latosa: double,
    -- specific leaf area (m2/kgC)
    sla: double,
    -- leaf longevity (years)
    leaflong: double,
    -- leaf to root mass ratio under non-water-stressed conditions
    ltor_max: double,
    -- litter moisture flammability threshold (fraction of AWC)
    litterme: double,
    -- fire resistance (0-1)
    fireresist: double,
    -- minimum forest-floor PAR level for growth (grasses) or establishment (trees)
    -- J/m2/day, individual and cohort modes */
    parff_min: double,
    -- parameter capturing non-linearity in recruitment rate relative to
    -- understorey growing conditions for trees (Fulton 1991) (individual and
    -- cohort modes)

    alphar: double,
    -- maximum sapling establishment rate (saplings/m2/year) (individual and cohort modes)
    est_max: double,
    -- constant used in calculation of sapling establishment rate when spatial
    -- mass effect enabled (individual and cohort modes)

    kest_repr: double,
    -- constant affecting amount of background establishment
    -- \see ifbgestab */
    kest_bg: double,
    -- constant used in calculation of sapling establishment rate when spatial
    -- mass effect disabled (individual and cohort modes)

    kest_pres: double,
    -- expected longevity under non-stressed conditions (individual and cohort modes)
    longevity: double,
    -- threshold growth efficiency for imposition of growth suppression mortality
    -- kgC/m2 leaf/year, individual and cohort modes */
    greff_min: double,

    -- Bioclimatic limits (all temperatures deg C)

    -- minimum 20-year coldest month mean temperature for survival
    tcmin_surv: double,
    -- maximum 20-year coldest month mean temperature for establishment
    tcmax_est: double,
    -- minimum degree day sum on 5 deg C base for establishment
    gdd5min_est: double,
    -- minimum 20-year coldest month mean temperature for establishment
    tcmin_est: double,
    -- minimum warmest month mean temperature for establishment
    twmin_est: double,
    -- continentality parameter for boreal summergreen trees
    twminusc: double,
    -- constant in equation for budburst chilling time requirement (Sykes et al 1996)
    k_chilla: double,
    -- coefficient in equation for budburst chilling time requirement
    k_chillb: double,
    -- exponent in equation for budburst chilling time requirement
    k_chillk: double,
    -- array containing values for GDD0(c) given c=number of chill days
    -- Sykes et al 1996, Eqn 1
    -- gdd0 has one element for each possible value for number of chill days

    ---- FIXME TODO HACK! comes from Date::MAX_YEAR_LENGTH+1
    gdd0: double[Date_MAX_YEAR_LENGTH_plusone],

    -- interception coefficient (unitless)
    intc: double,

    -- the amount of N that is applied (kg N m-2)
    N_appfert: double,
    -- 0 - 1 how much of the fertiliser is applied the first date, default 1.
    fertrate: (double, double),
    -- dates relative to sowing date
    fertdates: (int, int),
    fert_stages: (double, double),
    fertilised: (bool, bool),

    T_vn_min: double,
    T_vn_opt: double,
    T_vn_max: double,

    T_veg_min: double,
    T_veg_opt: double,
    T_veg_max: double,

    T_rep_min: double,
    T_rep_opt: double,
    T_rep_max: double,

    photo: (double, double, double),

    dev_rate_veg: double,
    dev_rate_rep: double,

    a1: double, b1: double, c1: double, d1: double, a2: double, b2: double, c2: double, d2: double, a3: double, b3: double, c3: double, d3: double,
    cton_stem_avr: double,
    cton_stem_max: double,

    -- Drought tolerance level (0 = very -> 1 = not at all) (unitless)
    -- Used to implement drought-limited establishment */
    drought_tolerance: double,

    -- bvoc

    -- aerodynamic conductance (m s-1)
    ga: double,
    -- isoprene emission capacity (ug C g-1 h-1)
    eps_iso: double,
    -- whether (1) or not (1) isoprene emissions show a seasonality
    seas_iso: bool,
    -- monoterpene emission capacity (ug C g-1 h-1) per monoterpene species
    eps_mon: double[NMTCOMPOUNDS],
    -- fraction of monoterpene production that goes into storage pool (-) per monoterpene species
    storfrac_mon: double[NMTCOMPOUNDS],

    -- Bioclimatic limits parameters from Wolf et al. 2008

    -- snow max [mm]
    max_snow: double,
    -- snow min [mm]
    min_snow: double,
    -- GDD0 min
    gdd0_min: double,
    -- GDD0 max
    gdd0_max: double,

    -- New parameters from parameters from Wania et al. (2009a, 2009b, 2010)

    -- Days per month for which inundation is tolerated
    inund_duration: int,
    -- Inundation stress is felt when the water table (mm) is above wtp_max
    wtp_max: double,
    -- Whether this PFT has aerenchyma through which O2 and CH4 can be transported (Wania et al. 2010 - Sec 2.6)
    has_aerenchyma: bool,

    -- Sapling/regeneration characteristics (used only in population mode)
    -- For trees, on sapling individual basis (kgC) for grasses, on stand area basis,
    -- kgC/m2 */

    -- leaf C biomass
    regen_cmass_leaf: double,
    -- fine root C biomass
    regen_cmass_root: double,
    -- sapwood C biomass
    regen_cmass_sap: double,
    -- heartwood C biomass
    regen_cmass_heart: double,

    -- specifies type of landcover pft is allowed to grow in (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND)
    landcover: landcovertype,
    -- pft selection
    selection: xtring,
    -- fraction of residue outtake at harvest
    res_outtake: double,
    -- harvest efficiency
    harv_eff: double,
    -- harvest efficiency of intercrop grass
    harv_eff_ic: double,
    -- fraction of harvested products that goes into patchpft.harvested_products_slow
    harvest_slow_frac: double,
    -- yearly turnover fraction of patchpft.harvested_products_slow (goes to gridcell.acflux_harvest_slow)
    turnover_harv_prod: double,
    -- whether pft may grow as cover crop
    isintercropgrass: bool,
    -- whether autumn temperature dependent sowing date is calculated
    ifsdautumn: bool,
    -- upper temperature limit for autumn sowing
    tempautumn: double,
    -- lower temperature limit for spring sowing
    tempspring: double,
    -- default length of growing period
    lgp_def: int,
    -- upper minimum temperature limit for crop sowing
    maxtemp_sowing: double,
    -- default sowing date in the northern hemisphere (julian day)
    sdatenh: int,
    -- default sowing date in the southern hemisphere
    sdatesh: int,
    -- whether sowing date adjusting equation is used
    sd_adjust: bool,
    -- parameter 1 in sowing date adjusting equation
    sd_adjust_par1: double,
    -- parameter 2 in sowing date adjusting equation
    sd_adjust_par2: double,
    -- parameter 3 in sowing date adjusting equation
    sd_adjust_par3: double,
    -- latest date for harvesting in the northern hemisphere
    hlimitdatenh: int,
    -- latest date for harvesting in the southern hemisphere
    hlimitdatesh: int,
    -- default base temperature (°C) for heat unit (hu) calculation
    tb: double,
    -- temperature under which vernalisation is possible (°C)
    trg: double,
    -- default number of vernalising days required
    pvd: int,
      -- sensitivity to the photoperiod effect [0-1]
    psens: double,
    -- basal photoperiod (h) (pb<ps for longer days plants)
    pb: double,
    -- lag in days after sowing before vernalization starts
    vern_lag: int,
    -- saturating photoperiod (h) (ps<pb for shorter days plants)
    ps: double,
    -- default potential heat units required for crop maturity (degree-days)
    phu: double,
    -- whether quadratic equation used for calculating potential heat units (Bondeau method)
    phu_calc_quad: bool,
    -- whether linear equation used for calculating potential heat units (Bondeau method)
    phu_calc_lin: bool,
    -- minimum potential heat units required for crop maturity (Bondeau method) (degree-days)
    phu_min: double,
    -- maximum potential heat units required for crop maturity (Bondeau method) (degree-days)
    phu_max: double,
    -- reduction factor of potential heat units in spring crops (Bondeau method) (degree-days)
    phu_red_spring_sow: double,
    -- number of days of phu decrease in the linear phu equation (Bondeau method)
    ndays_ramp_phu: double,
    -- intercept for the linear phu equation (Bondeau method)
    phu_interc: double,
    -- fraction of growing season (phu) at which senescence starts [0-1]
    fphusen: double,
    -- type of senescence curve (see Bondeau et al. 2007)
    shapesenescencenorm: bool,
    -- fraction of maximal LAI still present at harvest [0-1]
    flaimaxharvest: double,
    -- default maximum LAI (only used for intercrop grass in the case where no pasture is present in any stand)
    laimax: double,
    -- whether harvestable organs are above ground
    aboveground_ho: bool,
    -- optimum harvest index
    hiopt: double,
    -- minimum harvest index
    himin: double,
    -- initial fraction of growing season's npp allocated to roots
    frootstart: double,
    -- final fraction of growing season's npp allocated to roots
    frootend: double,
    -- autumn/spring sowing of pft:s with tempautumn = 1
    forceautumnsowing: int,  --0 = NOFORCING,  1 = AUTUMNSOWING, 2 = SPRINGSOWING
    -- N limited version of pft
    nlim: bool

  }

  let avg_cton (min: double, max: double) : double =
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
      shapesenescencenorm = 0,
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

      photo = (0.0, 0.0, 0.0)

    } in pft

  -- Calculates SLA given leaf longevity
  let initsla(Pft: pft) : Pft =

    -- SLA has to be supplied in the insfile for crops with N limitation
    if (!(phenology == CROPGREEN && ifnlim)) then
      -- Reich et al 1992, Table 1 (includes conversion x2.0 from m2/kg_dry_weight to
      -- m2/kgC)
      pft with pft.sla =
      if (leafphysiognomy == BROADLEAF) then
        0.2 * pow(10.0, 2.41 - 0.38 * log10(12.0 * pft.leaflong))
      else if (leafphysiognomy == NEEDLELEAF) then
        0.2 * pow(10.0, 2.29 - 0.4 * log10(12.0 * pft.leaflong))
      else pft.sla
    else pft

  -- Calculates minimum leaf C:N ratio given leaf longevity
  let init_cton_min(Pft: pft) : Pft =
    -- cton_leaf_min has to be supplied in the insfile for crops with N limitation
    if (!(phenology == CROPGREEN && ifnlim)) then
      -- Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
      -- kgN/kgC)
      pft with pft.cton_leaf_min =
      if (leafphysiognomy == BROADLEAF) then
        500.0 / pow(10.0, 1.75 - 0.33 * log10(12.0 * leaflong))
      else if (leafphysiognomy == NEEDLELEAF) then
        500.0 / pow(10.0, 1.52 - 0.26 * log10(12.0 * leaflong))
      else cton_leaf_min
    else pft

  let init_cton_limits(Pft: pft) : Pft =
    -- Fraction between min and max C:N ratio White et al. 2000
    let frac_mintomax = if (phenology == CROPGREEN && ifnlim) then 5.0 else 2.78  -- Use value also without nlim ?

    -- Fraction between leaf and root C:N ratio
    let frac_leaftoroot = 1.16 -- Friend et al. 1997

    -- Fraction between leaf and sap wood C:N ratio
    let frac_leaftosap = 6.9   -- Friend et al. 1997

    -- Max leaf C:N ratio
    let pft = pft with pft.cton_leaf_max = pft.cton_leaf_min * frac_mintomax

    -- Average leaf C:N ratio
    let pft = pft with pft.cton_leaf_avr = pft.avg_cton(pft.cton_leaf_min, pft.cton_leaf_max)

    -- Tighter C:N ratio range for roots and sapwood: picked out thin air
    let frac_maxtomin = 0.9

    -- Maximum fine root C:N ratio
    let pft = pft with pft.cton_root_max = pft.cton_leaf_max * frac_leaftoroot

    let cton_root_min = pft.cton_root_max * pft.frac_maxtomin

    -- Average fine root C:N ratio
    let pft = pft with pft.cton_root_avr = avg_cton(cton_root_min, cton_root_max)

    -- Maximum sap C:N ratio
    let pft = pft with pft.cton_sap_max  = cton_leaf_max * frac_leaftosap

    let cton_sap_min = cton_sap_max * frac_maxtomin

    -- Average sap C:N ratio
    let pft = pft with pft.cton_sap_avr  = avg_cton(cton_sap_min, cton_sap_max)

    let pft = pft with pft.respcoeff =
    if (lifeform == GRASS || lifeform == MOSS) then
      respcoeff /= 2.0 * cton_root / (cton_root_avr + cton_root_min)
    else
      respcoeff /= cton_root / (cton_root_avr + cton_root_min) +
                   cton_sap  / (cton_sap_avr  + cton_sap_min)
    let pft = pft with pft.cton_stem_max = 1.0/(2.0*0.0034) --Maize params
    let pft = pft with pft.cton_stem_avr = 1.0/(2.0*0.0068)
    in ptf

  -- Calculates coefficient to compensate for different vertical distribution of fine root on nitrogen uptake
  let init_nupscoeff(Pft: pft) : Pft =
    -- Fraction fine root in upper soil layer should have higher possibility for mineralized nitrogen uptake
    -- Soil nitrogen profile is considered to have a exponential decline (Franzluebbers et al. 2009) giving
    -- an approximate advantage of 2 of having more roots in the upper soil layer
    let upper_adv = 2.0

    -- Simple solution until we get C and N in all soil layers.
    let rootdist_upper = 0.0
    let rootdist_lower = 0.0

    let (sl, rootdist_upper, rootdist_lower) =
    loop (sl, rootdist_upper, rootdist_lower) = (0, 0.0, 0.0)
    while (sl < NSOILLAYER) do
      if (sl < NSOILLAYER_UPPER) then
        (sl+1, rootdist_upper+=pft.rootdist[sl], rootdist_lower) else
        (sl+1, rootdist_upper, rootdist_lower+=pft.rootdist[sl])
    in pft with pft.nupscoeff = rootdist_upper * upper_adv + rootdist_lower


  -- Initialises sapling/regen characteristics in population mode following LPJF formulation
  let initregen(Pft: pft) : Pft =

    -- see function allometry in growth module.

    -- Note: primary PFT parameters, including SLA, must be set before this
    --       function is called

    let REGENLAI_TREE = 1.5
    let REGENLAI_GRASS = 0.001
    let SAPLINGHW = 0.2

    in if (lifeform == TREE) then
      -- Tree sapling characteristics
      let pft = pft with pft.regen_cmass_leaf =
        pow(REGENLAI_TREE * k_allom1 * pow(1.0 + SAPLINGHW, k_rp) *
          pow(4.0 * sla / PI / k_latosa, k_rp * 0.5) / sla, 2.0 / (2.0 - k_rp))
      let pft = pft with pft.regen_cmass_leaf =
        wooddens * k_allom2 * pow((1.0 + SAPLINGHW) *
        sqrt(4.0 * regen.cmass_leaf * sla / PI / k_latosa), k_allom3) *
        regen.cmass_leaf * sla / k_latosa
      in pft with pft.regen_cmass_heart = SAPLINGHW * regen_cmass_sap

    else if (lifeform == GRASS || lifeform==MOSS) then

      -- Grass regeneration characteristics
      pft with pft.regen_cmass_leaf = REGENLAI_GRASS / sla

    else
      pft with pft.regen_cmass_root = 1.0 / ltor_max * pft.regen_cmass_leaf

  -- Inits root fractions in each soil layer through a shape parameter beta (see Jackson et al., 1996)
  let init_rootdist(Pft: pft) : Pft =

    let depth = Dz_soil * CM_PER_MM
    let rootdist = pft.rootdist with [0] = 1.0 - pow(root_beta, depth) -- init first layer
    let tot = rootdist[0]

    let (i, depth, rootdist, tot) =
    loop (i, depth, rootdist, tot) =
         (1, depth, rootdist, tot) while (i<NSOILLAYER) do
      ( i+1
      , depth + Dz_soil * CM_PER_MM
      , rootdist with [i] =  1.0 - pow(root_beta, depth) - (1.0 - pow(root_beta, depth - Dz_soil * CM_PER_MM))
      , tot + rootdist[i])

    -- Calibrated the root_beta for each PFT to match rootdist_upper from 'old' (pre LPJG 4.1) parameterisation.
    -- Sometimes the rootdist goes below our maximum soildepth. When that happens, put the residual fraction in lowest soil layer
    in pft with pft.rootdist = rootdist with [NSOILLAYER-1] = [NSOILLAYER-1] + 1.0 - tot

  let ismoss(pft: Pft) : bool = pft.lifeform == MOSS
  let isgrass(pft: Pft) : bool = pft.lifeform == GRASS
  let istree(pft: Pft) : bool = pft.lifeform == TREE
  let iswetlandspecies(pft: Pft) : bool = (ismoss(pft) || pft.has_aerenchyma)
  let Pft
