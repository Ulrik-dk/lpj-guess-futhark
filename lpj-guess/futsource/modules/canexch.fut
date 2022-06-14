-------------------------------------Combined imports-----------------------------------
open import "../framework/guess"
--#include "config.h"
--#include "canexch.h"
--#include "driver.h"
open import "../modules/q10"
open import "../futhark-extras"
--#include "bvoc.h"
open import "../modules/ncompete"
--#include <assert.h>
-------------------------------------PART FROM .h--------------------------------------
-- Constants for photosynthesis calculations
-- conversion factor for solar radiation at 550 nm from J/m2 to mol_quanta/m2 (E=mol quanta) mol J-1
let CQ : real = 4.6e-6
-- intrinsic quantum efficiency of CO2 uptake, C3 plants
let ALPHA_C3 : real = 0.08
-- intrinsic quantum efficiency of CO2 uptake, C4 plants
let ALPHA_C4 : real = 0.053
-- O2 partial pressure (Pa)
let PO2 : real = 2.09e4
-- colimitation (shape) parameter
let THETA : real = 0.7
-- 'saturation' ratio of intercellular to ambient CO2 partial pressure for C4 plants
let LAMBDA_SC4 : real = 0.4
-- leaf respiration as fraction of maximum rubisco, C3 plants
let BC3 : real = 0.015
-- leaf respiration as fraction of maximum rubisco, C4 plants
let BC4 : real = 0.02
-- leaf respiration as fraction of maximum rubisco, mosses
-- see Wania et al. (2009b)
let BC_moss : real = 0.03
let CMASS : real = 12.0    -- atomic mass of carbon
let ALPHAA : real = 0.5    -- value chosen to give global carbon pool and flux values that
                -- agree with published estimates.
                -- scaling factor for PAR absorption from leaf to plant projective area level
                -- alias "twigloss". Should normally be in the range 0-1
let ALPHAA_NLIM : real = 0.65 -- Same as ALPHAA above but chosen to give pools and flux values
                 -- that agree with published estimates when Nitrogen limitation is
                 -- switched on.
let ALPHAA_CROP : real = 0.7      -- Value for crops without N limitation.
let ALPHAA_CROP_NLIM : real = 0.9  -- Value for crops with N limitation
-- Lambert-Beer extinction law (Prentice et al 1993 Monsi & Saeki 1953)
let lambertbeer (lai : real) : real = exp(-0.5 * lai)
-------------------------------------PART FROM .cpp------------------------------------
-- leaf nitrogen (kgN/kgC) not associated with photosynthesis
-- (value given by Haxeltine & Prentice 1996b)---
let N0 : real = 7.15 * G_PER_MG
-- Lookup tables for parameters with Q10 temperature responses
-- lookup table for Q10 temperature response of CO2/O2 specificity ratio
let lookup_tau = LookupQ10(0.57, 2600.0)
-- lookup table for Q10 temperature response of Michaelis constant for O2
let lookup_ko = LookupQ10(1.2, 3.0e4)
-- lookup table for Q10 temperature response of Michaelis constant for CO2
let lookup_kc = LookupQ10(2.1, 30.0)
--void interception(Patch& patch,Climate& climate) {
--------------------------------------------------------------------------------------/

-- replacing the respective call with these will result in futhark being unable to typecheck it
let type_checker_helper1 (i: Individual, ppfts: [npft]Patchpft) : Individual =
    if ppfts[i.pft_id].cropphen.growingseason then
      let i = i with fpar_leafon = 0.0
      let i = i with fpar = 0.0
      in i
    else i

let type_checker_helper2 (i: Individual, value: real) : Individual =
  i with lai_leafon_layer = value

let type_checker_helper3 (indiv: Individual, plai_leafon_layer: real, fpar_uptake_leafon_layer : real, plai_layer : real, fpar_uptake_layer : real) : Individual =
  let indiv : Individual =
  if (!negligible(plai_leafon_layer)) then
    -- FPAR partitioned according to the relative amount
    -- of leaf area in this layer for this individual
    indiv with fpar_leafon = indiv.fpar_leafon + fpar_uptake_leafon_layer*indiv.lai_leafon_layer/plai_leafon_layer
  else indiv
  let indiv =
    if (!negligible(plai_layer)) then
      indiv with fpar = indiv.fpar + fpar_uptake_layer* (indiv.lai_leafon_layer*indiv.phen)/plai_layer
    else indiv
  in indiv

let type_checker_helper4 (spft: Standpft, ps_res: PhotosynthesisResult) : Standpft =
  spft with photosynthesis_result = ps_res

let type_checker_helper5 (i: Individual, ps_res: PhotosynthesisResult) : Individual =
  i with photosynthesis_result = ps_res

let type_checker_helper6 (i: Individual, subdaily : int) : Individual =
  let new_gpterms = copy i.gpterms
  let new_gpterms[subdaily] = realzero
  in i with gpterms = new_gpterms

let type_checker_helper7 (i: Individual, ci: NCompetingIndividual) : Individual =
  i with fnuptake = ci.fnuptake

let type_checker_helper8 (i: Individual) : Individual =
  i


-- Internal function - not intended to be called by framework
let fpar(patch : Patch, vegetation : [npft]Individual, climate : Climate, pfts : [npft]Pft, date : Date, ppfts : [npft]Patchpft)
  : (Patch, Vegetation) =
  -- DESCRIPTION
  -- Calculates daily fraction of incoming PAR (FPAR) taken up by individuals in a
  -- particular patch over their projective areas, given current leaf phenological
  -- status. Calculates PAR and FPAR at top of grass canopy (individual and cohort
  -- modes). Calculates fpar assuming leaf-on (phen=1) for all vegetation.
  --
  -- Note: In order to compensate for the additional computational cost of
  --       calculating fpar_leafon in cohort/individual mode, the grain of the
  --       integration of FPAR through the canopy has been increased from 1 to 2 m
  --
  -- NEW ASSUMPTIONS CONCERNING FPC AND FPAR (Ben Smith 2002-02-20)
  -- FPAR = average individual fraction of PAR absorbed on patch basis today,
  --        including effect of current leaf phenology (this differs from previous
  --        versions of LPJ-GUESS in which FPAR was on an FPC basis)
  -- FPC =  PFT population (population mode), cohort (cohort mode) or individual
  --        (individual mode) fractional projective cover as a fraction of patch area
  --        (in population mode, corresponds to LPJF variable fpc_grid). Updated
  --        annually based on leaf-out LAI (see function allometry in growth module).
  --        (FPC was previously equal to summed crown area as a fraction of patch
  --        area in cohort/individual mode)
  --
  -- Population mode: FPAR on patch (grid cell) area basis assumed to be equal to fpc
  -- under full leaf cover i.e.
  --     (1) fpar        = fpc*phen
  --     (2) fpar_leafon = fpc
  --
  -- Individual and cohort modes: FPAR calculated assuming trees shade themselves
  --   and all individuals below them according to the Lambert-Beer law (Prentice
  --   et al 1993, Eqn 27 Monsi & Saeki 1953):
  --     (3) fpar = integral [0-tree height] exp ( -k * plai(z) )
  --   where
  --       k       = extinction coefficient
  --       plai(z) = summed leaf-area index for leaves of all individuals, above
  --                 canopy depth z, taking account of current phenological status
  let VSTEP : real = 2.0 -- width of vertical layers for canopy-area integration (m)
  let PHEN_GROWINGSEASON : real = 0.5
    -- minimum expected vegetation leaf-on fraction for growing season
  --double plai -- cumulative leaf-area index (LAI) for patch (m2 leaf/m2 ground)
  --double plai_leafon
    -- cumulative LAI for patch assuming full leaf cover for all individuals
  --double plai_layer -- summed LAI by layer for patch
  --double plai_leafon_layer
    -- summed LAI by layer for patch assuming full leaf cover for all individuals
  --double plai_grass -- summed LAI for grasses
  --double plai_leafon_grass
    -- summed LAI for grasses assuming full leaf cover for all individuals
  --double flai -- fraction of total grass LAI represented by a particular grass
  --double fpar_layer_top -- FPAR by layer
  --double fpar_leafon_layer_top
    -- FPAR by layer assuming full leaf cover for all individuals
  --double fpar_layer_bottom
  --double fpar_leafon_layer_bottom
  --double fpar_grass -- FPAR at top of grass canopy
  --double fpar_leafon_grass
    -- FPAR at top of grass canopy assuming full leaf cover for all individuals
  --double fpar_ff -- FPAR at forest floor (beneath grass canopy)
  --double fpar_leafon_ff
    -- FPAR at forest floor assuming full leaf cover for all individuals
  --double frac
    -- vertical fraction of layer occupied by crown cylinder(s) of a particular
    -- individual or cohort
  --double atoh -- term in calculating LAI sum for a given layer
  --double height_veg -- maximum vegetation height (m)
  --int toplayer -- number of vertical layers of width VSTEP in vegetation (minus 1)
  --int layer -- layer number (0=lowest)
  --double lowbound -- lower bound of current layer (m)
  --double highbound -- upper bound of current layer (m)
  --double fpar_min -- minimum FPAR required for grass growth
  --double par_grass -- PAR reaching top of grass canopy (J/m2/day)
  --double phen_veg -- LAI-weighted mean fractional leaf-out for vegetation
  --variables needed for "S�kes" FPAR scheme
  --double fpar_uptake_layer
  --double fpar_uptake_leafon_layer
  -- Obtain reference to Vegetation object
  --Vegetation& vegetation=patch.vegetation
  -- And to Climate object
  --const Climate& climate = patch.get_climate()
  in
    if (vegmode==POPULATION) then
      -- POPULATION MODE
      -- Loop through individuals
      (patch, map (\indiv ->
          let indiv = indiv with fpar = fpc_today(indiv, pfts[indiv.pft_id], ppfts[indiv.pft_id])
          let indiv = indiv with fpar_leafon = if ppfts[indiv.pft_id].cropphen.growingseason then indiv.fpar else 0.0
          in indiv
        ) vegetation)
        -- For this individual ...
    else
      -- INDIVIDUAL OR COHORT MODE
      -- Initialise individual FPAR, find maximum height of vegetation, calculate
      -- individual LAI given current phenology, calculate summed LAI for grasses
      let plai : real = 0.0
      let plai_leafon : real = 0.0
      let plai_grass : real = 0.0
      let plai_leafon_grass : real = 0.0


      let height_veg : real =
        reduce (curry max) 0.0 <|
        map (\i -> if ppfts[i.pft_id].cropphen.growingseason then i.height else 0.0) vegetation
      let vegetation = map (\i -> type_checker_helper1(i, ppfts)) vegetation
      let plai_leafon = reduce (+) 0.0  <| map (\i -> if ppfts[i.pft_id].cropphen.growingseason then i.lai else 0.0) vegetation
      let (plai_grass, plai_leafon_grass) =
            reduce (\(a1,b1) (a2,b2) -> (a1+a2,b1+b2)) (0.0, 0.0)
            <| map (\i -> if pfts[i.pft_id].lifeform == GRASS || pfts[i.pft_id].lifeform == MOSS then (i.lai, lai_today(i, pfts[i.pft_id], ppfts[i.pft_id])) else (0.0, 0.0)) vegetation
      -- Accumulate LAI-weighted sum of individual leaf-out fractions
      let phen_veg : real = reduce (+) 0.0 <| map (\i -> lai_today(i, pfts[i.pft_id], ppfts[i.pft_id])) vegetation



      -- Calculate LAI-weighted mean leaf-out fraction for vegetation
      -- guess2008 - bugfix - was: if (!negligible(plai))
      let phen_veg = if (!negligible(plai_leafon)) then phen_veg / plai_leafon else 1.0

      -- Calculate number of layers (minus 1) from ground surface to top of canopy
      let toplayer = intFromReal (height_veg/VSTEP-0.0001)

      -- Calculate FPAR by integration from the top of the canopy (Eqn 2)
      let plai : real = 0.0
      let plai_leafon : real = 0.0

      -- Set FPAR for bottom of layer above (initially 1 at top of canopy)
      let fpar_layer_bottom = 1.0
      let fpar_leafon_layer_bottom = 1.0

      let  (layer,plai,plai_leafon,vegetation) =
      loop (layer,plai,plai_leafon,vegetation)
         = (toplayer,plai,plai_leafon,vegetation)
      while (layer >=0) do

        let lowbound = realFromInt layer*VSTEP
        let highbound = lowbound+VSTEP

        -- FPAR at top of this layer = FPAR at bottom of layer above
        let fpar_layer_top = fpar_layer_bottom
        let fpar_leafon_layer_top = fpar_leafon_layer_bottom

        -- Loop through individuals
        let (vegetation, pls, plls) = unzip3 <| map (\indiv ->
          if (pfts[indiv.pft_id].lifeform==TREE) then
            if (indiv.height>lowbound && indiv.boleht<highbound &&
              !negligible(indiv.height-indiv.boleht)) then
              -- Calculate vertical fraction of current layer occupied by
              -- crown cylinders of this cohort
              let frac = 1.0
              let frac =
              if (indiv.height<highbound) then
                frac-(highbound-indiv.height)/VSTEP
              else frac
              let frac =
              if (indiv.boleht>lowbound) then
                frac-(indiv.boleht-lowbound)/VSTEP
              else frac
              -- Calculate summed LAI of this cohort in this layer
              let atoh = indiv.lai/(indiv.height-indiv.boleht)
              let indiv = type_checker_helper2 (indiv,atoh*frac*VSTEP)
              let pl = indiv.lai_leafon_layer * indiv.phen
              let pll = indiv.lai_leafon_layer
              in (indiv, pl, pll)
            else
              let indiv = indiv with lai_layer = 0.0
              let indiv = indiv with lai_leafon_layer = 0.0
              in (indiv, 0.0, 0.0)
          else (indiv, 0.0, 0.0)) vegetation
        let plai_layer : real = reduce (+) 0.0 pls
        let plai_leafon_layer : real = reduce (+) 0.0 plls

        -- Calculate FPAR at bottom of this layer
        -- Eqn 27, Prentice et al 1993
        let fpar_layer_bottom = lambertbeer(plai)
        let fpar_leafon_layer_bottom = lambertbeer(plai_leafon)
        -- Total PAR uptake in this layer
        let fpar_uptake_layer = fpar_layer_top-fpar_layer_bottom
        let fpar_uptake_leafon_layer = fpar_leafon_layer_top-fpar_leafon_layer_bottom
        -- Partition PAR for this layer among trees,
        let vegetation = map (\indiv ->
          if (pfts[indiv.pft_id].lifeform==TREE) then
            type_checker_helper3 (indiv, plai_leafon_layer, fpar_uptake_leafon_layer, plai_layer, fpar_uptake_layer)
          else indiv
        ) vegetation
        in (layer-1,plai + plai_layer,plai_leafon + plai_leafon_layer,vegetation)

      -- FPAR reaching grass canopy
      let fpar_grass = lambertbeer(plai)
      let fpar_leafon_grass = lambertbeer(plai_leafon)
      -- Add grass LAI to calculate PAR reaching forest floor
      -- BLARP: Order changed Ben 050301 to overcome optimisation bug in pgCC
      --plai+=plai_grass
      let fpar_ff = lambertbeer(plai+plai_grass)
      let plai = plai + plai_grass
      -- Save this
      let patch = patch with fpar_ff = fpar_ff
      let plai_leafon = plai_leafon + plai_leafon_grass
      let fpar_leafon_ff = lambertbeer(plai_leafon)
      -- FPAR for grass PFTs is difference between relative PAR at top of grass canopy
      -- canopy and at forest floor, or lower if FPAR at forest floor below threshold
      -- for grass growth. PAR reaching the grass canopy is partitioned among grasses
      -- in proportion to their LAI (a somewhat simplified assumption)
      -- Loop through individuals
      let (vegetation) = map (\indiv ->
        if (pfts[indiv.pft_id].lifeform==GRASS || pfts[indiv.pft_id].lifeform==MOSS) then
          -- Calculate minimum FPAR for growth of this grass
          -- Fraction of total grass LAI represented by this grass
          let flai =
          if (!negligible(plai_grass)) then
            lai_today(indiv, pfts[indiv.pft_id], ppfts[indiv.pft_id])/plai_grass
          else 1.0
          let fpar_min =
          if (!negligible(climate.par)) then
            min(pfts[indiv.pft_id].parff_min/climate.par,1.0)
          else 1.0
          let indiv = indiv with fpar = max(0.0,fpar_grass*flai-max(fpar_ff*flai,fpar_min))
          -- Repeat assuming full leaf cover for all individuals
          let flai =
          if (!negligible(plai_leafon_grass)) then
            indiv.lai/plai_leafon_grass
          else 1.0
          let indiv = indiv with fpar_leafon=max(0.0,fpar_leafon_grass*flai-
            max(fpar_leafon_ff*flai,fpar_min))
          in indiv
        else (indiv)
        ) vegetation
      let fpar_tree_total : real =
        reduce (+) 0.0 <| map
          (\i -> if pfts[i.pft_id].lifeform==TREE then i.fpar else 0.0
          ) vegetation
      -- Save grass canopy FPAR and update mean growing season grass canopy PAR
      -- Growing season defined here as days when mean vegetation leaf-on fraction
      -- exceeds 50% and we're in the light half of the year (daylength >= 11).
      --
      -- The daylength condition was added because sites with evergreens can have
      -- a mean vegetation leaf-on fraction over 50% even during polar night.
      -- 11 hours was chosen because some sites never reach exactly 12 hours, the
      -- exact limit shouldn't matter much.
      let patch = patch with fpar_grass=fpar_grass
      let par_grass = fpar_grass * climate.par
      let patch =
          if (date.day==0) then
            let patch = patch with par_grass_mean = 0.0
            let patch = patch with nday_growingseason = 0
            in patch
          else patch
      let patch =
          if (phen_veg > PHEN_GROWINGSEASON && climate.daylength >= 11.0) then
            let patch = patch with par_grass_mean = patch.par_grass_mean + par_grass
            let patch = patch with nday_growingseason = patch.nday_growingseason + 1
            in patch
          else patch
      -- Convert from sum to mean on last day of year
      let patch =
          if (date.islastday && date.islastmonth && (patch.nday_growingseason > 0)) then
            patch with par_grass_mean = patch.par_grass_mean / realFromInt patch.nday_growingseason
          else patch
      in (patch, vegetation)

let alphaa(pft : Pft) =
  if (pft.phenology == CROPGREEN) then
    if ifnlim then ALPHAA_CROP_NLIM else ALPHAA_CROP
  else
    if ifnlim then ALPHAA_NLIM else ALPHAA
--- Non-water stressed rubisco capacity, with or without nitrogen limitation
let vmax(b : real,
         c1 : real,
         c2 : real,
         apar : real,
         tscal : real,
         daylength : real,
         temp : real,
         nactive : real,
         ifnlimvmax : bool)
         : (real, real, real) =
  -- Calculation of non-water-stressed rubisco capacity assuming leaf nitrogen not
  -- limiting (Eqn 11, Haxeltine & Prentice 1996a)
  -- Calculation of sigma is based on Eqn 12 (same source)
  let s : real =  24.0 / daylength * b
  let sigma : real = sqrt(max(0.0, 1.0 - (c2 - s) / (c2 - THETA * s)))
  let vm = 1 / b * CMASS * CQ * c1 / c2 * tscal * apar *
              (2.0 * THETA * s * (1.0 - sigma) - s + c2 * sigma)
  -- Calculate nitrogen-limited Vmax for current leaf nitrogen
  -- Haxeltine & Prentice 1996b Eqn 28
  let M : real = 25.0 -- corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996b
  -- Conversion factor in calculation of leaf nitrogen: includes conversion of:
  --    - Vm from gC/m2/day to umolC/m2/sec
  --      - nitrogen from mg/m2 to kg/m2
  let CN : real = 1.0 / (3600 * daylength * CMASS)
  let tfac : real = exp(-0.0693 * (temp - 25.0))
  let vm_max : real = nactive / (M * CN * tfac)
  -- Calculate optimal leaf nitrogen based on [potential] Vmax (Eqn 28 Haxeltine & Prentice 1996b)
  let nactive_opt : real = M * vm * CN * tfac
  in if (vm > vm_max && ifnlimvmax) then
      let vmaxnlim = vm_max / vm  -- Save vmax nitrogen limitation
      let vm = vm_max
      in (vm, vmaxnlim, nactive_opt)
    else
      (vm, 1.0, nactive_opt)
--- Total daily gross photosynthesis
--- Calculation of total daily gross photosynthesis and leaf-level net daytime
--  photosynthesis given degree of stomatal closure (as parameter lambda).
--  Includes implicit scaling from leaf to plant projective area basis.
--  Adapted from Farquhar & von Caemmerer (1982) photosynthesis model, as simplified
--  by Collatz et al (1991), Collatz et al (1992), Haxeltine & Prentice (1996a,b)
--  and Sitch et al. (2000).
--
--  To calculate vmax call w/ daily averages of temperature and par.
--  Vmax is to be calculated daily and only with lambda == lambda_max.
--  lambda values greater than lambda_max are forbidden.
--  In sub-daily mode daylength should be 24 h, to obtain values in daily units.
--
--  INPUT PARAMETERS
--
--  \param PhotosynthesisEnvironment struct containing the following public members:
--   - co2        atmospheric ambient CO2 concentration (ppmv)
--   - temp       mean air temperature today (deg C)
--   - par        total daily photosynthetically-active radiation today (J/m2/day)
--   - daylength  day length, must equal 24 in diurnal mode (h)
--   - fpar       fraction of PAR absorbed by foliage
--
--  \param PhotosynthesisStresses struct containing the following members:
--   - ifnlimvmax      - whether nitrogen should limit Vmax
--   - moss_ps_limit    - limit to moss photosynthesis. [0,1], where 1 means no limit
--   - graminoid_ps_limit  - limit to graminoid photosynthesis. [0,1], where 1 means no limit
--   - inund_stress      - limit to photosynthesis due to inundation, where 1 means no limit
--
--  \param pft        Pft object containing the following public members:
--   - pathway         biochemical pathway for photosynthesis (C3 or C4)
--   - pstemp_min      approximate low temperature limit for photosynthesis (deg C)
--   - pstemp_low      approximate lower range of temperature optimum for
--                     photosynthesis (deg C)
--   - pstemp_high     approximate upper range of temperature optimum for photosynthesis
--                     (deg C)
--   - pstemp_max      maximum temperature limit for photosynthesis (deg C)
--   - lambda_max      non-water-stressed ratio of intercellular to ambient CO2 pp
--
--  \param lambda     ratio of intercellular to ambient partial pressure of CO2
--
--  \param nactive    nitrogen available for photosynthesis
--
--  \param vm         pre-calculated value of Vmax for this stand for this day if
--                    available, otherwise calculated
--
-- OUTPUT PARAMETERS
--
-- \param ps_result      see documentation of PhotosynthesisResult struct
--
--
-- IMPORTANT for users adding new call parameters to the list above:
--
-- Never place new call parameters in the proper photosynthesis() function header, instead
-- place new parameters insided the structs PhotosynthesisEnvironment and PhotosynthesisStresses,
-- or in PhotosynthesisResult if it is a result.
let photosynthesis(ps_env : PhotosynthesisEnvironment,
                   ps_stresses : PhotosynthesisStresses,
                   pft: Pft,
                   lambda : real,
                   nactive : real,
                   vm : real)
                   : PhotosynthesisResult =
  -- NOTE: This function is identical to LPJF subroutine "photosynthesis" except for
  -- the formulation of low-temperature inhibition coefficient tscal (tstress LPJF).
  -- The function adopted here draws down metabolic activity in approximately the
  -- temperature range pstemp_min-pstemp_low but does not affect photosynthesis
  -- at high temperatures.
  -- HISTORY
  -- Ben Smith 18/1/2001: Tested in comparison to LPJF subroutine "photosynthesis":
  -- function showed identical behaviour except at temperatures >= c. 35 deg C where
  -- LPJF temperature inhibition function results in lower photosynthesis.
  -- Make sure that only two alternative modes are possible:
  --  * daily non-water stressed (forces Vmax calculation)
  --  * with pre-calculated Vmax (sub-daily and water-stressed)
  --assert(vm >= 0 || lambda == pft.lambda_max) --TODO
  --assert(lambda <= pft.lambda_max) --TODO
  let PATMOS : real = 1e5  -- atmospheric pressure (Pa)
  -- Get the environmental variables
  let temp : real = ps_env.temp
  let co2 : real = ps_env.co2
  let fpar : real = ps_env.fpar
  let par : real = ps_env.par
  let daylength : real = ps_env.daylength
  -- Get the stresses
  let ifnlimvmax : bool = ps_stresses.ifnlimvmax
  let pstemp_max = pft.pstemp_max
  let pstemp_high = pft.pstemp_high
  let pstemp_low = pft.pstemp_low
  let pstemp_min = pft.pstemp_min
  let pathway = pft.pathway
  let lifeform = pft.lifeform
  -- No photosynthesis during polar night, outside of temperature range or no RuBisCO activity
  in if (negligible(daylength) || negligible(fpar) || temp > pstemp_max || temp < pstemp_min || !(boolFromReal vm)) then PhotosynthesisResult() else
  -- Scale fractional PAR absorption at plant projective area level (FPAR) to
  -- fractional absorption at leaf level (APAR)
  -- Eqn 4, Haxeltine & Prentice 1996a
  let apar : real = par * fpar * alphaa(pft)
  -- Calculate temperature-inhibition coefficient
  -- This function (tscal) is mathematically identical to function tstress in LPJF.
  -- In contrast to earlier versions of modular LPJ and LPJ-GUESS, it includes both
  -- high- and low-temperature inhibition.
  let k1 : real = (pstemp_min+pstemp_low) / 2.0
  let tscal : real = (1.0 - 0.01*exp(4.6/(pstemp_max-pstemp_high)*(temp-pstemp_high)))/
                    (1.0+exp((k1-temp)/(k1-pstemp_min)*4.6))
  let (b, c1, c2) : (real, real, real) =
    if (pathway == C3) then      -- C3 photosynthesis
      -- Calculate CO2 compensation point (partial pressure)
      -- Eqn 8, Haxeltine & Prentice 1996a
      let gammastar : real = PO2 / 2.0 / readQ10(temp, lookup_tau)
      -- Intercellular partial pressure of CO2 given stomatal opening (Pa)
      -- Eqn 7, Haxeltine & Prentice 1996a
      let pi_co2 : real = lambda * co2 * PATMOS * CO2_CONV
      -- Calculation of C1_C3, Eqn 4, Haxeltine & Prentice 1996a
      -- High-temperature inhibition modelled by suppression of LUE by decreased
      -- relative affinity of rubisco for CO2 with increasing temperature (Table 3.7,
      -- Larcher 1983)
      -- Notes: - there is an error in Eqn 4, Haxeltine & Prentice 1996a (missing
      --          2.0* in denominator) which is fixed here (see Eqn A2, Collatz
      --          et al 1991)
      --        - the explicit low temperature inhibition function has been removed
      --          and replaced by a temperature-dependent upper limit on V_m, see
      --          below
      --        - the reduction in maximum photosynthesis due to leaf age (phi_c)
      --          has been removed
      --        - alpha_a, accounting for reduction in PAR utilisation efficiency
      --          from the leaf to ecosystem level, appears in the calculation of
      --          apar (above) instead of here
      --        - C_mass, the atomic weight of carbon, appears in the calculation
      --          of V_m instead of here
      let c1 : real = (pi_co2 - gammastar) / (pi_co2 + 2.0 * gammastar) * ALPHA_C3
      -- Calculation of C2_C3, Eqn 6, Haxeltine & Prentice 1996a
      let c2 : real = (pi_co2 - gammastar) / (pi_co2 + readQ10(temp, lookup_kc) * (1.0 + PO2/readQ10(temp, lookup_ko)))
      let b : real = if lifeform == MOSS then BC_moss else BC3
      -- see Wania et al. 2009b
      in (b, c1, c2)
    else               -- C4 photosynthesis
      -- Calculation of C1_C4 given actual pi (lambda)
      -- C1_C4 incorporates term accounting for effect of intercellular CO2
      -- concentration on photosynthesis (Eqn 14, 16, Haxeltine & Prentice 1996a)
      let c1 = min(lambda/LAMBDA_SC4, 1.0) * ALPHA_C4
      let c2 = 1.0
      let b = BC4
      in (b, c1, c2)
  -- Calculation of non-water-stressed rubisco capacity (Eqn 11, Haxeltine & Prentice 1996a)
  let ps_result = PhotosynthesisResult()
  let (vm, vmaxnlim, nactive_opt) =
    if (vm < 0) then vmax(b, c1, c2, apar, tscal, daylength, temp, nactive, ifnlimvmax)
    else (vm, ps_result.vmaxnlim, ps_result.nactive_opt)
  -- Calculation of daily leaf respiration
  -- Eqn 10, Haxeltine & Prentice 1996a
  let rd_g = vm * b
  -- PAR-limited photosynthesis rate (gC/m2/h)
  -- Eqn 3, Haxeltine & Prentice 1996a
  let je = (c1 * tscal * apar * CMASS * CQ / daylength)
  -- Rubisco-activity limited photosynthesis rate (gC/m2/h)
  -- Eqn 5, Haxeltine & Prentice 1996a
  let jc : real = c2 * vm / 24.0
  -- Calculation of daily gross photosynthesis
  -- Eqn 2, Haxeltine & Prentice 1996a
  -- Notes: - there is an error in Eqn 2, Haxeltine & Prentice 1996a (missing
  --       theta in 4*theta*je*jc term) which is fixed here
  let agd_g = ((je + jc - sqrt((je + jc) * (je + jc) - 4.0 * THETA * je * jc)) /(2.0 * THETA) * daylength)
  let (agd_g, rd_g) = if (!iftwolayersoil) then
    -- LIMITS TO PHOTOSYNTHESIS
    -- On wetlands, both agd_g and rd_g are scaled in the event of inundation (all PFTS),
    -- or, for mosses and graminoids only, in the event of the water table dropping below
    -- the PFT's optimal water table depth (dessication). See Wania et al. (2009b)
    -- Get the stresses
    let moss_ps_limit : real = ps_stresses.moss_ps_limit
    let graminoid_ps_limit : real = ps_stresses.graminoid_ps_limit
    let inund_stress : real = ps_stresses.inund_stress
    -- 1) Inundation stress
    -- Reduce GPP if there is inundation stress
    -- (possibility of) inund stress (i.e. values < 1) only on PEATLAND stands and when ifinundationstress is true
    let agd_g = agd_g * inund_stress
    let rd_g = rd_g * inund_stress
    -- 2a) Moss dessication
    let (agd_g, rd_g) =
      -- Reduce agd_g using moss_wtp_limit (lies between [0.3, 1.0])
      if (lifeform == MOSS) then (agd_g * moss_ps_limit, rd_g * moss_ps_limit)
      else (agd_g, rd_g)
    -- 2b) Graminoid dessication (NB! all other PFTs have graminoid_ps_limit == 1)
    let agd_g = agd_g * graminoid_ps_limit
    let rd_g = rd_g * graminoid_ps_limit
    in (agd_g, rd_g)
  else (agd_g, rd_g)
    -- Leaf-level net daytime photosynthesis (gC/m2/day)
  -- Based on Eqn 19, Haxeltine & Prentice 1996a
  let adt : real = agd_g - daylength / 24.0 * rd_g
  -- Convert to CO2 diffusion units (mm/m2/day) using ideal gas law
  let adtmm = (adt / 12 * 8.314 * (temp + 273.15) / 100000 * 1000)
  in {agd_g = agd_g
     ,adtmm = adtmm
     ,je = je
     ,rd_g = rd_g
     ,vm = vm
     ,nactive_opt = nactive_opt
     ,vmaxnlim = vmaxnlim}
-- Calculate value for canopy conductance component associated with photosynthesis (mm/s)
-- Eqn 21, Haxeltine & Prentice 1996
-- includes conversion of daylight from hours to seconds
let gpterm (adtmm: real, co2: real, lambda: real, daylength: real) : real =
  if (adtmm <= 0)
    then 0.0
    else 1.6 / CO2_CONV / 3600 * adtmm / co2 / (1 - lambda) / daylength
--- Determine CO2 in peatland watindiv with fnuptake = 0.0er
---
-- Updated daily
---
let get_co2(p : Patch, climate : Climate, pft : Pft, stand: Stand, soil : Soil, g: Gridcell) : real =
  if (is_highlatitude_peatland_stand(stand, g) && ismoss(pft))
  then soil.acro_co2 -- override for peat mosses
  else climate.co2
--- Determine inundation stress
---
-- Updated daily
-- --Was given a Patch, and accessed stand by parental-reference, rather just ask the stand directly
let get_inund_stress(ppft : Patchpft, stand : Stand, g: Gridcell) : real =
  if (ifinundationstress && is_highlatitude_peatland_stand(stand, g))
  then ppft.inund_stress -- Determine inundation stress for this Patchpft
  else 1.0 -- No stress by default
--- Determine limits on graminoid photosynthesis due to WTP
---
---Updated daily
let get_graminoid_wtp_limit(p : Patch, pft : Pft, stand : Stand, soil : Soil, gridcell : Gridcell) : real =
  if (!ismoss(pft) && iswetlandspecies(pft) && is_highlatitude_peatland_stand(stand, gridcell))
  then soil.dgraminoid_wtp_limit -- Update limit if this is a graminoid
  else 1.0 -- No limit by default
--- Determine limits on moss photosynthesis due to WTP
--
-- Updated daily
---
let get_moss_wtp_limit(pft : Pft, stand : Stand, soil : Soil, gridcell : Gridcell) : real =
  -- Update limit if this is a moss
  if (ismoss(pft) && is_highlatitude_peatland_stand(stand, gridcell))
  then soil.dmoss_wtp_limit
  else 1.0 -- No limit by default
--/ Pre-calculate Vmax and no-stress assimilation and canopy conductance
---
---Vmax is calculated on a daily scale (w/ daily averages of temperature and par)
---Subdaily values calculated if needed
 --
let photosynthesis_nostress(vegetation : Vegetation, patch : Patch, climate : Climate, spfts : [npft]Standpft, ppfts : [npft]Patchpft, pfts : [npft]Pft, stand : Stand, soil : Soil, gridcell : Gridcell, date : Date)
 : (Vegetation, [npft]Standpft, [npft]Patchpft) =
 --ps_stress.no_stress() --TODO: is this just the default constructor?
 -- If this is the first patch, calculate no-stress assimilation for
 -- each Standpft, assuming FPAR=1. This is then later used in
 -- forest_floor_conditions.
 let pftco2 : real = climate.co2 -- will override for peat mosses
 let (spfts, ppfts) =
    if patch.patch_id == 0 then
     unzip <| map2 (\spft ppft ->
       if (spft.active) then
         let pftco2 = get_co2(patch, climate, pfts[spft.pft_id], stand, soil, gridcell)
         let ps_env = {co2=pftco2, temp=climate.temp, par=climate.par, fpar=1.0, daylength=climate.daylength}
         let ps_stress = {ifnlimvmax=false, moss_ps_limit=get_moss_wtp_limit(pfts[spft.pft_id], stand, soil, gridcell), graminoid_ps_limit=get_graminoid_wtp_limit(patch, pfts[spft.pft_id], stand, soil, gridcell), inund_stress=get_inund_stress(ppft, stand, gridcell)}
         -- Call photosynthesis assuming stomates fully open (lambda = lambda_max)
         let ps_result = photosynthesis(ps_env, ps_stress, pfts[spft.pft_id],
                                        pfts[spft.pft_id].lambda_max, 1.0, -1)
         let spft = type_checker_helper4(spft, ps_result)
         in (spft, ppft)
       else (spft, ppft)
     ) spfts ppfts
   else (spfts, ppfts)
 -- Pre-calculation of no-stress assimilation for each individual
 let (vegetation) =
   map (\indiv ->
     let pft : Pft = pfts[indiv.pft_id]
     let ppft : Patchpft = ppfts[pft.pft_id]
     let pftco2 = get_co2(patch, climate, pft, stand, soil, gridcell)
     let ps_env = {co2=pftco2, temp=climate.temp, par=climate.par, fpar=indiv.fpar, daylength=climate.daylength}
     let ps_stress = {ifnlimvmax=false, moss_ps_limit=get_moss_wtp_limit(pft, stand, soil, gridcell), graminoid_ps_limit=get_graminoid_wtp_limit(patch, pft, stand, soil, gridcell), inund_stress=get_inund_stress(ppft, stand, gridcell)}
     -- Individual photosynthesis with no nitrogen limitation
     let ps_result = photosynthesis(ps_env, ps_stress, pft,
                                    pft.lambda_max, 1.0, -1)
     let indiv = type_checker_helper5(indiv, ps_result)
     let indiv = indiv with gpterm = gpterm(indiv.photosynthesis_result.adtmm, pftco2, pft.lambda_max, climate.daylength)
     let indiv =
     if (diurnal(date)) then
       let indiv = type_checker_helper6(indiv, date.subdaily)
       --FIXME is this correct?
       let res = PhotosynthesisResult()
       let phots = copy indiv.phots
       let phots[date.subdaily] = res
       let indiv = indiv with phots = phots
       let (indiv, _) =
       loop (indiv, i) =
            (indiv, 0)
       while (i < date.subdaily) do
         let ps_result : PhotosynthesisResult = copy indiv.phots[i]
         -- Update temperature and PAR
         let ps_env : PhotosynthesisEnvironment = {co2=pftco2, temp=climate.temps[i], par=climate.pars[i], fpar=indiv.fpar, daylength=24}
         let ps_result = photosynthesis(ps_env, ps_stress, pft, pft.lambda_max, 1.0, indiv.photosynthesis_result.vm)
         let gpterms = copy indiv.gpterms
         let gpterms[i] = gpterm(ps_result.adtmm, climate.co2, pft.lambda_max, 24)
         let indiv = indiv with gpterms = gpterms
         in (indiv, i+1)
         in indiv
       else indiv
       in indiv
   ) vegetation
 in (vegetation, spfts, ppfts)
--- Calculates individual fnuptake based on surface of fine root
--- Calculates individual fraction nitrogen uptake based on surface of fine root
--- Roots are cone formed with height == radius.
--- V = PI * r^3 / 3
--- A = (2^1/2 + 1) * PI * r^2
--- -> A = const * cmass_root^2/3
---
let nitrogen_uptake_strength(indiv : Individual, pfts: [npft]Pft, patchpft: Patchpft) : real =
 pow(max(0.0, cmass_root_today(indiv, pfts[indiv.pft_id], patchpft)) * pfts[indiv.pft_id].nupscoeff * indiv.cton_status / indiv.densindiv, 2.0 / 3.0) * indiv.densindiv
--- Individual nitrogen uptake fraction
--- Determining individual nitrogen uptake as a fraction of its nitrogen demand.
--- \see ncompete
--- Function nitrogen_uptake_strength() determines how good individuals are at
--- acquiring nitrogen.
---

let fnuptake [n] (vegetation : [n]Individual, nmass_avail : real, pfts: [npft]Pft, patchpft: Patchpft) =
  -- Create vector describing the individuals to ncompete()
  let competingIndividuals = map (\i -> {ndemand = vegetation[i].ndemand, strength = nitrogen_uptake_strength(vegetation[i], pfts, patchpft), fnuptake = 0}) <| iota n
  -- Let ncompete() do the actual distribution
  let individuals : [n]NCompetingIndividual = ncompete(competingIndividuals, nmass_avail)

  -- Get the results, nitrogen uptake fraction for each individual
  let vegetation = map (type_checker_helper7) <| zip vegetation individuals
  in vegetation
--- Use nitrogen storage to limit stress

--- Retranslocated nitrogen from last year is used to
--- limit nitrogen stress in leaves, roots, and sap wood
let nstore_usage(vegetation : [npft]Individual, pfts: [npft]Pft, patchpft: Patchpft) : [npft]Individual =
  map (\indiv ->
    -- individual excess nitrogen demand after uptake
    let indiv = type_checker_helper8(indiv)
    let excess_ndemand : real = (indiv.leafndemand + indiv.rootndemand) * (1.0 - indiv.fnuptake) + indiv.leafndemand_store + indiv.rootndemand_store
    in
    -- if individual is in need of using its labile nitrogen storage
    if (!negligible(excess_ndemand)) then
      -- if labile nitrogen storage is larger than excess nitrogen demand
      if (excess_ndemand <= indiv.nstore_labile) then
        -- leaf nitrogen demand
        let leaf_ndemand : real = indiv.leafndemand * (1.0 - indiv.fnuptake) + indiv.leafndemand_store
        let indiv = indiv with nmass_leaf = indiv.nmass_leaf + leaf_ndemand
        let indiv = indiv with nstore_labile = indiv.nstore_labile - leaf_ndemand

        -- root nitrogen demand
        let root_ndemand : real = indiv.rootndemand * (1.0 - indiv.fnuptake) + indiv.rootndemand_store
        let indiv = indiv with nmass_root = indiv.nmass_root + root_ndemand
        let indiv = indiv with nstore_labile = indiv.nstore_labile - root_ndemand
        let indiv = indiv with nstress = false
        in indiv
      else
        if (!negligible(indiv.nstore_labile)) then
          -- calculate total nitrogen mass
          let tot_nmass : real = indiv.nmass_leaf + indiv.nmass_root + indiv.fnuptake * (indiv.leafndemand + indiv.rootndemand) + indiv.nstore_labile

          -- new leaf C:N ratio
          let cton_leaf : real = (cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) + cmass_root_today(indiv, pfts[indiv.pft_id], patchpft) * (pfts[indiv.pft_id].cton_leaf_avr / pfts[indiv.pft_id].cton_root_avr)) / tot_nmass

          -- nitrogen added to leaf from storage
          let labile_nto_leaf : real = cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) / cton_leaf - (indiv.nmass_leaf + indiv.fnuptake * indiv.leafndemand)

          -- new leaf nitrogen
          let indiv = indiv with nmass_leaf = indiv.nmass_leaf + labile_nto_leaf

          -- new root nitrogen
          let indiv = indiv with nmass_root = indiv.nmass_root + indiv.nstore_labile - labile_nto_leaf
          let indiv = indiv with nstore_labile = 0.0
          in indiv
        -- nitrogen stressed photosynthesis is allowed only when nitrogen limitation is turned on
        else
          let indiv = indiv with nstress = ifnlim
          in indiv
    else
      let indiv = indiv with nstress = false
      in indiv
      -- photosynthesis will not be nitrogen stresses
  ) vegetation
--/ Nitrogen demand

--- Determines nitrogen demand based on vmax for leaves.
--  Roots and sap wood nitrogen concentration follows leaf
--  nitrogen concentration.
--  Also determines individual nitrogen uptake capability
--/
let ndemand(patch : Patch, vegetation : [npft]Individual, gridcell : Gridcell, soil : Soil, pfts : [npft]Pft, patchpft : Patchpft)
  : (Patch, Vegetation) =
  --  Gridcell& gridcell = patch.stand.get_gridcell()
  --  Soil& soil = patch.soil
  --/ daily nitrogen demand for patch (kgN/m2)
  let patch = patch with ndemand = 0.0
  -- soil available mineral nitrogen (kgN/m2)
  let nmin_avail : real = nmass_avail(soil,NO)
  -- Scalar to soil temperature (Eqn A9, Comins & McMurtrie 1993) for nitrogen uptake
  let soilT = get_soil_temp_25(soil)
  let temp_scale = if soilT > 0.0 then max(0.0, 0.0326 + 0.00351 * pow(soilT, 1.652) - pow(soilT / 41.748, 7.19)) else 0.0
  --/ Rate of nitrogen uptake not associated with Michaelis-Menten Kinetics (Zaehle and Friend 2010)
  let kNmin : real = 0.05
  let vegetation =
  map (\indiv ->
    let indiv = type_checker_helper8(indiv)
    -- Rescaler of nitrogen uptake
    let indiv = indiv with fnuptake = 1.0
    -- Starts with no nitrogen stress
    let indiv = indiv with nstress = false
    -- Optimal leaf nitrogen content
    --double leafoptn
    -- Optimal leaf C:N ratio
    --double cton_leaf_opt
    -- Calculate optimal leaf nitrogen content and demand
    let (leafoptn, cton_leaf_opt, indiv) =
    if (!negligible(indiv.phen)) then
      let indiv = indiv with nday_leafon = indiv.nday_leafon + 1
      -- Added a scalar depending on individual lai to slow down light optimization of newly shaded leafs
      -- Peltoniemi et al. 2012
      let indiv = indiv with nextin = exp(0.12 * min(10.0*indiv.phen, lai_indiv_today(indiv, pfts[indiv.pft_id], patchpft)))
      -- Calculate optimal leaf nitrogen associated with photosynthesis and none photosynthetic
      -- active nitrogen (Haxeltine et al. 1996 eqn 27/28)
      let leafoptn = indiv.photosynthesis_result.nactive_opt * indiv.nextin + N0 * cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft)
      -- Can not have higher nitrogen concentration than minimum leaf C:N ratio

      let leafoptn =
      if (cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) / leafoptn < pfts[indiv.pft_id].cton_leaf_min) then
        cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) / pfts[indiv.pft_id].cton_leaf_min
      -- Can not have lower nitrogen concentration than maximum leaf C:N ratio
      else if (cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) / leafoptn > pfts[indiv.pft_id].cton_leaf_max) then
        cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) / pfts[indiv.pft_id].cton_leaf_max
      else leafoptn

      -- Updating annual optimal leaf C:N ratio
      let indiv = indiv with cton_leaf_aopt = min(cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) / leafoptn, indiv.cton_leaf_aopt)
      -- Leaf nitrogen demand
      let indiv = indiv  with leafndemand = max(leafoptn - indiv.nmass_leaf, 0.0)
      -- Setting daily optimal leaf C:N ratio
      let cton_leaf_opt =
        if (indiv.leafndemand > 0) then
          cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) / leafoptn
        else
          max(pfts[indiv.pft_id].cton_leaf_min, cton_leaf(indiv, true, stand, standpft, pfts, patchpft))
      in (leafoptn, cton_leaf_opt, indiv)
    else
      let indiv = indiv with leafndemand = 0.0
      let cton_leaf_opt = cton_leaf(indiv, true, stand, standpft, pfts, patchpft)
      in (leafoptn, cton_leaf_opt, indiv)

    -- Nitrogen demand
    -- Root nitrogen demand
    let indiv = indiv with rootndemand = max(0.0, cmass_root_today(indiv, pfts[indiv.pft_id], patchpft) / (cton_leaf_opt * pfts[indiv.pft_id].cton_root_avr / pfts[indiv.pft_id].cton_leaf_avr) - indiv.nmass_root)
    -- Sap wood nitrogen demand. Demand is ramped up throughout the year.
    let indiv =
    if (pfts[indiv.pft_id].lifeform == TREE) then
      indiv with sapndemand = max(0.0, indiv.cmass_sap / (cton_leaf_opt * pfts[indiv.pft_id].cton_sap_avr / pfts[indiv.pft_id].cton_leaf_avr) - indiv.nmass_sap) * ((1.0 + (double)date.day)/date.year_length())
    else indiv
    -- Labile nitrogen storage demand
    let indiv = indiv with storendemand = indiv.ndemand_storage(cton_leaf_opt)
    --TODO HO demand
    let indiv = indiv with hondemand = 0.0
    -- Total nitrogen demand
    let ndemand_tot : real = indiv.leafndemand + indiv.rootndemand + indiv.sapndemand + indiv.storendemand + indiv.hondemand
    -- Calculate scalars to possible nitrogen uptake
    -- Current plant mobile nitrogen concentration
    let ntoc : real =  if !negligible(cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) + cmass_root_today(indiv, pfts[indiv.pft_id], patchpft)) then (indiv.nmass_leaf + indiv.nmass_root) / (cmass_leaf_today(indiv, pfts[indiv.pft_id], patchpft) + cmass_root_today(indiv, pfts[indiv.pft_id], patchpft)) else 0.0
    -- Scale to maximum nitrogen conceforest_floorntrations
    let indiv = indiv with cton_status = max(0.0, (ntoc - 1.0 / pfts[indiv.pft_id].cton_leaf_min) / (1.0 / pfts[indiv.pft_id].cton_leaf_avr - 1.0 / pfts[indiv.pft_id].cton_leaf_min))
    -- Nitrogen availablilty scalar due to saturating Michealis-Menten kinetics
    let nmin_scale : real = kNmin + nmin_avail / (nmin_avail + gridcell.pft[pfts[indiv.pft_id].id].Km)
    -- Maximum available soil mineral nitrogen for this individual is base on its root area.
    -- This is considered to be related to FPC which is proportional to crown area which is approx
    -- 4 times smaller than the root area
    let max_indiv_avail : real = min(1.0, indiv.fpc * 4.0) * nmin_avail
    -- Maximum nitrogen uptake due to all scalars (times 2 because considering both NO3- and NH4+ uptake)
    -- and soil available nitrogen within individual projectived coverage
    let maxnup : real = min(2.0 * pfts[indiv.pft_id].nuptoroot * nmin_scale * temp_scale * indiv.cton_status * cmass_root_today(indiv, pft, patchpft), max_indiv_avail)
    -- Nitrogen demand limitation due to maximum nitrogen uptake capacity
    let fractomax : real = if ndemand_tot > 0.0 then min(maxnup/ndemand_tot,1.0) else 0.0
    -- Root and leaf demand from storage pools
    let indiv = indiv with leafndemand_store = indiv.leafndemand * (1.0 - fractomax)
    let indiv = indiv with rootndemand_store = indiv.rootndemand * (1.0 - fractomax)
    -- Nitrogen demand after adjustment to maximum uptake capacity
    let indiv = indiv with leafndemand  = indiv.leafndemand * fractomax
    let indiv = indiv with rootndemand  = indiv.rootndemand * fractomax
    let indiv = indiv with sapndemand   = indiv.sapndemand * fractomax
    let indiv = indiv with storendemand = indiv.storendemand * fractomax
    -- Sum total nitrogen demand individual is capable of taking up
    let indiv = indiv with ndemand = indiv.leafndemand + indiv.rootndemand + indiv.sapndemand + indiv.storendemand
    -- Negative nitrogen demand not allowed
    let indiv =
    if (indiv.ndemand <= 0.0) then
      let indiv = indiv with ndemand = 0.0
      -- Compartments fraction of total nitrogen demand
      let indiv = indiv with leaffndemand  = 0.0
      let indiv = indiv with rootfndemand  = 0.0
      let indiv = indiv with sapfndemand   = 0.0
      let indiv = indiv with storefndemand = 0.0
      in indiv
    else
      -- Compartments fraction of total nitrogen demand
      let indiv = indiv with leaffndemand  = indiv.leafndemand / indiv.ndemand
      let indiv = indiv with rootfndemand  = indiv.rootndemand / indiv.ndemand
      let indiv = indiv with sapfndemand   = indiv.sapndemand  / indiv.ndemand
      let indiv = indiv with storefndemand = max(0.0, 1.0 - (indiv.leaffndemand + indiv.rootfndemand + indiv.sapfndemand))
      in indiv
    in indiv
  ) vegetation
  -- Sum total patch nitrogen demand
  let patch = patch with ndemand = patch.ndemand + (reduce (+) 0.0 <| map (\indiv -> indiv.ndemand) vegetation)
  in (patch, vegetation)
--/ Recalculation of photosynthesis under vmax nitrogen stress
--- If nitrogen supply is not able to meet demand it will lead
--  to down-regulation of vmax resulting in lower photosynthesis
--/
let vmax_nitrogen_stress(patch : Patch, climate : Climate, vegetation : [npft]Individual, soil : Soil, pfts : [npft]Pft, ppfts : [nptf]Patchpft)
  : [npft]Individual =
  -- Supply function for nitrogen and determination of nitrogen stress leading
  -- to down-regulation of vmax.
  -- Nitrogen within projective cover of all individuals
  let tot_nmass_avail : real = soil.nmass_avail(NO) * min(1.0, patch.fpc_total)
  -- Calculate individual uptake fraction of nitrogen demand
  let vegetation =
  if (patch.ndemand > tot_nmass_avail) then
    -- Determine individual nitrogen uptake fractions
    fnuptake(vegetation, tot_nmass_avail)
  else vegetation
  -- Resolve nitrogen stress with longterm stored nitrogen
  let vegetation = nstore_usage(vegetation)

  -- Calculate leaf nitrogen associated with photosynthesis, nitrogen limited photosynthesis,
  -- and annual otimal leaf nitrogen content and nitrogen limitation on vmax
  let vegetation =
  map (\ indiv ->
    let pft : Pft = pfts[indiv.pft_id]
    let ppft : Patchpft = ppfts[pft.pft_id]
    -- Calculate leaf nitrogen associated with photosynthesis (Haxeltine et al. 1996 eqn 27/28)
    -- Added difference between needleleaved and broadleaved mentioned in Friend et al. 1997
    -- Should be done on FPC basis, but is not as it does not matter mathematically
    -- Needs to be calculated for each individual due to possible water stress
    -- Todays leaf nitrogen after uptake
    let nmass_leaf = indiv.nmass_leaf + indiv.leafndemand * indiv.fnuptake
    let indiv = indiv with nactive =
    if (indiv.phen > 0.0) then
      max(0.0, nmass_leaf - N0 * cmass_leaf_today(indiv, pft, patchpft))
    else 0.0
    -- Individuals photosynthesis is nitrogen stressed
    let indiv =
    if (indiv.nstress) then
      let pftco2 : real = get_co2(patch, climate, pft, stand, soil, gridcell)
      let ps_env : PhotosynthesisEnvironment = {co2 = pftco2, temp=climate.temp, par=climate.par, fpar=indiv.fpar, daylength=climate.daylength}
      -- Set stresses
      PhotosynthesisStresses ps_stress
      let ps_stress = {ifnlimvmax=true, moss_ps_limit=get_moss_wtp_limit(pft, stand, soil, gridcell), graminoid_ps_limit=get_graminoid_wtp_limit(patch, pft, stand, soil, gridcell), inund_stress=get_inund_stress(ppft, stand, gridcell)}
      -- Individual photosynthesis
      let ps_result = photosynthesis(ps_env, ps_stress, pft,
                                     pft.lambda_max, indiv.nactive / indiv.nextin, -1,
                                     indiv.photosynthesis_result)
      let indiv = indiv with photosynthesis = ps_result

      let indiv = indiv with gpterm = gpterm(indiv.photosynthesis_result.adtmm, pftco2, pft.lambda_max, climate.daylength)
      let indiv =
      if (date.diurnal()) then
        let (indiv, _) =
        loop (indiv, i)
          = (indiv, 0)
        while (i < date.subdaily) do
          let ps_result = indiv.phots[i]
          let ps_env = {co2=pftco2, temp=climate.temps[i], pars=climate.pars[i], fpar=indiv.fpar, daylength=24}
          let ps_result = photosynthesis(ps_env, ps_stress, pft,
                          pft.lambda_max, indiv.nactive / indiv.nextin, indiv.photosynthesis_result.vm,
                          ps_result)
          let indiv = indiv with gpterms = gpterms with [i] = gpterm(ps_result.adtmm, climate.co2, pft.lambda_max, 24)
          in indiv
        in indiv
        else indiv
      in
      -- Sum annual average nitrogen limitation on vmax
      if (indiv.phen) then
        let indiv = indiv with avmaxnlim = indiv.avmaxnlim + indiv.photosynthesis_result.vmaxnlim
      -- On last day of year determined average nitrogen limitation
      -- based on days with leaf on
        in
        if (date.islastday && date.islastmonth) then
          if (!negligible(indiv.nday_leafon)) then
            indiv with avmaxnlim = indiv.avmaxnlim / realFromInt indiv.nday_leafon
          else
            indiv with avmaxnlim = 0.0
        else indiv
      else indiv
    else indiv
    in indiv
  ) vegetation
  in vegetation
--/ Transpirative demand

--- Two alternative parameterisations of aet_monteith are available:
--      AET_MONTEITH_HYPERBOLIC and AET_MONTEITH_EXPONENTIAL
--  \see canexch.h
--/
let wdemand(patch : Patch, climate : Climate, vegetation : [npft]Individual, day : Day) : (Patch, Vegetation) =
    -- non-water-stressed canopy conductance assuming full leaf cover, patch
    -- vegetated area basis (mm/s)
  let (vegetation, gps, glps) = map (\indiv ->
    let pft : Pft = pfts[indiv.pft_id]
    -- Calculate non-water-stressed canopy conductance assuming full leaf cover
    --        - include canopy-conductance component not linked to
    --          photosynthesis (diffusion through leaf cuticle etc) this is
    --          assumed to be proportional to leaf-on fraction
    -- Call photosynthesis for individual assuming stomates fully open
    -- (lambda = lambda_max)
    in
    if (indiv.growingseason()) then
      let leafon_photosynthesis : PhotosynthesisResult = PhotosynthesisResult()
      -- Call photosynthesis first with fpar_leafon to get gp_leafon below.
      -- Should hopefully not be needed in future, demand_leafon only used
      -- by raingreen phenology.
      let temp : real = if date.diurnal() then climate.temps[day.period] else climate.temp
      let par : real = if date.diurnal() then climate.pars[day.period] else climate.par
      let daylength : real = if date.diurnal() then 24 else climate.daylength
      let pftco2 : real = get_co2(patch, climate, pft, stand, soil, gridcell)
      let ps_env = {co2=pftco2, temp=temp, par=par, fpar=indiv.fpar_leafon, daylight=daylight}
      let ppft : Patchpft = ppfts[patch.patchpft_id]
      let ps_stress = {ifnlimvmax=false, moss_ps_limit=get_moss_wtp_limit(pft, stand, soil, gridcell), graminoid_ps_limit=get_graminoid_wtp_limit(patch, pft, stand, soil, gridcell),inund_stress=get_inund_stress(ppft, stand, gridcell)}

      -- No nitrogen limitation when calculating gp_leafon
      let ps_result = photosynthesis(ps_env, ps_stress, pft,
                                     pft.lambda_max, 1.0, -1,
                                     leafon_photosynthesis)
      let gp_leafon : real = gpterm(leafon_photosynthesis.adtmm, pftco2, pft.lambda_max, daylength) + pft.gmin * indiv.fpc
      -- Increment patch sums of non-water-stressed gp by individual value
      let gp_patch = gp_patch + (if date.diurnal() then indiv.gpterms[day.period] else indiv.gpterm) + pft.gmin * fpc_today(indiv)
      let gp_leafon_patch = gp_leafon_patch + gp_leafon
      in (indiv, gp_patch, gp_leafon_patch)
    else (indiv, 0.0, 0.0)
  ) vegetation
  -- Determination of transpirative demand based on a Monteith parameterisation of
  -- boundary layer dynamics, i.e. demand = f(EET, conductance)
  let gp_patch : real = reduce (+) 0.0 gps
    -- non-water-stressed canopy conductance for patch, patch vegetated area
    -- basis (mm/s)
  let gp_leafon_patch : real = reduce (+) 0.0 glps

  -- Calculate transpirational demand on patch vegetated area basis
  -- Eqn 23, Haxeltine & Prentice 1996
  let (gp_patch, patch) =
    if (!negligible(gp_patch) && !negligible(patch.fpc_total)) then
      let gp_patch = gp_patch / patch.fpc_total
      let patch = patch with demand = aet_monteith(patch.eet_net_veg, gp_patch)
      in (gp_patch, patch)
    else
      let patch = patch with wdemand = 0.0
      let patch = patch with wdemand_day = patch.wdemand_day + patch.wdemand
      in (gp_patch, patch)
  let patch =
    if (day.isend) then
      patch with wdemand_day = patch.wdemand_day / date.subdaily
    else patch
  let patch =
    if (!negligible(gp_leafon_patch) && !negligible(patch.fpc_total)) then
      let gp_leafon_patch = gp_leafon_patch / patch.fpc_total
      let patch = patch with wdemand_leafon = aet_monteith(patch.eet_net_veg, gp_leafon_patch)
      in patch
    else patch with wdemand_leafon = 0.0
  in (patch, vegetation)
--double water_uptake(double wcont[NSOILLAYER], double awc[NSOILLAYER],
--double water_uptake_twolayer(double wcont[NSOILLAYER], double awc[NSOILLAYER],
--double irrigated_water_uptake(Patch& patch, Pft& pft, const Day& day) {
--/ Actual evapotranspiration and water stress
--- Soil water supply at the roots available to meet the transpirational demand
--  Fundamentally, water stress = supply < demand
--/
let aet_water_stress(patch : Patch, stand: Stand, vegetation : [npft]Individual, day : Day, pfts : [npft]Pft, spfts : [npft]Standpft, ppfts : [npft]Patchpft) : (Vegetation, [npft]Patchpft)  =
  -- Supply function for evapotranspiration and determination of water stress leading
  -- to down-regulation of stomatal conductance. Actual evapotranspiration determined
  -- as smaller of supply and transpirative demand (see function demand).
  -- Base value for actual canopy conductance calculated here for water-stressed
  -- individuals and used to derive actual photosynthesis in function npp (below)
  -- Calculate common point supply for each PFT in this patch
  let (p, ppfts) =
  loop (p, ppfts)
    = (0, ppfts)
  while (p<npft) do
    let spft : Standpft = spfts[p] in
    if (!spft.active) then (p+1) else
    -- Retrieve next patch PFT
    let ppft : Patchpft = ppfts[p]
    -- Retrieve PFT
    let pft : Pft = pfts[p]

    let (pft, ppft) =
    if (day.isstart || spft.irrigated && pft.pft_id == stand.pft_id) then
      -- Calculate effective water supply from plant roots
      -- Rescale available water by patch FPC if exceeds 1
      -- (this then represents the average amount of water available over an
      -- individual's FPC, assuming individuals are equal in competition for water)
      let wr : real =
        if (spft.irrigated && pft.id == patch.stand.pftid) then
          irrigated_water_uptake(patch, pft, day)
        else
          let wcont_local : [NSOILLAYER]real = copy soil.layer_soil_water
          in
          if (iftwolayersoil) then
            water_uptake_twolayer(wcont_local, patch.soil.soiltype.awc,
                  pft.rootdist, pft.emax, patch.fpc_rescale, ppft.fwuptake,
                  pft.lifeform == TREE, pft.drought_tolerance)
          else
            if (patch.stand.is_highlatitude_peatland_stand()) then -- Use awc_peat
              water_uptake(wcont_local, patch.soil.soiltype.awc_peat,
                  pft.rootdist, pft.emax, patch.fpc_rescale, ppft.fwuptake,
                  pft.lifeform == TREE, pft.drought_tolerance)
            else
              water_uptake(wcont_local, patch.soil.soiltype.awc,
                  pft.rootdist, pft.emax, patch.fpc_rescale, ppft.fwuptake,
                  pft.lifeform == TREE, pft.drought_tolerance)
      -- Calculate supply (Eqn 24, Haxeltine & Prentice 1996)
      let ppft = ppft with wsupply_leafon =
        if (patch.stand.landcover!=CROPLAND || ppft.cropphen.growingseason) then
          pft.emax * wr
        else 0.0
      let ppft = ppft with wsupply = ppft.wsupply_leafon * ppft.phen
      in (pft, ppft)
      else (pft, ppft)
    let ppft = ppft with wstress = ppft.wsupply < patch.wdemand && !negligible(ppft.phen) && !(pft.phenology==CROPGREEN && !largerthanzero(patch.wdemand-ppft.wsupply, -10))
    -- Calculate water-stressed canopy conductance on FPC basis assuming
    -- FPAR=1 and deducting canopy conductance component not associated
    -- with CO2 uptake valid for all individuals of this PFT in this patch
    -- today.
    -- Eqn 25, Haxeltine & Prentice 1996
    -- Fix, valid for monocultures, for faulty equation, manifesting itself in problems with crops in high scenario CO2-levels.
    -- No fix for natural vegetation yet.
    let gmin : real = if pft.phenology==CROPGREEN then ppft.phen * pft.gmin else pft.gmin
    let ppft = ppft with gcbase = if ppft.wstress then max(gc_monteith(ppft.wsupply, patch.eet_net_veg)-gmin * ppft.wsupply / patch.wdemand, 0.0) else 0
    let ppft =
    if (!date.diurnal()) then
      let ppft = ppft with wstress_day = ppft.wstress
      let ppft= ppft with gcbase_day = ppft.gcbase
      in ppft
    else if (day.isend) then
      let ppft = ppft with wstress_day = ppft.wsupply < patch.wdemand_day && !negligible(ppft.phen) && !(pft.phenology==CROPGREEN && !largerthanzero(patch.wdemand-ppft.wsupply, -10))
      let ppft = ppft with gcbase_day = if ppft.wstress_day then max(gc_monteith(ppft.wsupply,
          patch.eet_net_veg) - gmin * ppft.wsupply / patch.wdemand_day, 0.0) else 0
      in ppft
      else ppft
    in (p+1, ppfts with [p] = ppft)

  -- Calculate / transfer supply to individuals
  let vegetation = map(\indiv ->
    let ppft : Patchpft = patchpfts[indiv.pft_id]
    let indiv =
      if (day.isstart) then
        let indiv = indiv with aet = 0 in
        if (date.day == 0) then
          indiv with aaet = 0.0
        else indiv
      else indiv
    let indiv = indiv with wstress = ppft.wstress
    let indiv =
      if (indiv.alive) then
        if (indiv.wstress) then
          indiv with aet = indiv.aet + ppft.wsupply
        else indiv
      else
        indiv with aet = indiv.aet + if negligible(indiv.phen) then 0.0 else patch.wdemand
    let indiv = indiv with aet =
      if (day.isend) then
        indiv.aet * indiv.fpc / date.subdaily
      else aet
    let indiv = indiv with aaet =
      if (day.isend) then
        indiv.aaet + indiv.aet
      else indiv.aaet
    in indiv
  ) vegetation
  in (vegetation, ppfts)


--/ Water scalar
let water_scalar(stand : Stand, patch : Patch, vegetation : [npft]Individual, day : Day, ppfts : [npft]Patchpft, pfts : [nptf]Pft)
  : [nptf]Patchpft =
  -- Derivation of daily and annual versions of water scalar (wscal, omega)
  -- Daily version is used to determine leaf onset and abscission for raingreen PFTs.
  -- Annual version determines relative allocation to roots versus leaves for
  -- subsequent year
  let (_, ppfts) = loop (p, ppfts) = (0, ppfts) while (p<npft) do
    -- Retrieve next patch PFT
    let ppft : Patchpft = ppfts[p] in
      if (!patch.stand.pft[p].active) then (p+1, ppfts) else
    let pft = pfts[p]
    let ppft =
      if (day.isstart) then
      let ppft = ppft with wscal = 0 in
        if (date.day == 0) then
        let ppft = ppft with wscal_mean = 0 in
          if (pft.phenology==CROPGREEN || pft.isintercropgrass) then
            ppft with cropphen = cropphen with growingdays_y = 0
          else ppft
        else ppft
      else ppft
    let ppft = ppft with wscal =
      -- Calculate patch PFT water scalar value
      if (!negligible(patch.wdemand_leafon)) then
        ppft.wscal + min(1.0, ppft.wsupply_leafon/patch.wdemand_leafon)
      else
        ppft.wscal + 1.0
    let ppft =
      if (day.isend) then
        let ppft = ppft with wscal = ppft.wscal / realFromInt date.subdaily in
        if (stand.landcover!=CROPLAND || pft.phenology==ANY && ppft_pft_id==stand.pft_id) --natural, urban, pasture, forest and peatland stands
          then --normal grass growth
          let ppft = ppft with wscal_mean = ppft.wscal_mean + ppft.wscal in
            -- Convert from sum to mean on last day of year
            if (date.islastday && date.islastmonth) then
              ppft with wscal_mean = ppft.wscal_mean / date.year_length()
            else ppft
        else
          if (ppft.cropphen.growingseason || ppft.pft.phenology == CROPGREEN && date.day == ppft.cropphen.hdate || ppft.pft.isintercropgrass && date.day == patch.pft[stand.pft_id].cropphen.eicdate) then
            let ppft = ppft with cropphen = ppft.cropphen with growingdays_y = ppft.cropphen.growingdays_y+1
            let ppft = ppft with wscal_mean = ppft.wscal_mean + (ppft.wscal - ppft.wscal_mean) / ppft.cropphen.growingdays_y
            in ppft
          else ppft
      else ppft
    let ppfts = ppfts with [p] = ppft
    in (p+1, ppfts)
  in ppfts
-- ASSIMILATION_WSTRESS
-- Internal function (do not call directly from framework)
let assimilation_wstress
    (pft: Pft,
    co2: real,
    temp: real,
    par: real,
    daylength: real,
    fpar: real,
    fpc: real,
    gcbase: real,
    vmax: real,
    nactive: real,
    ifnlimvmax: bool,
    moss_wtp_limit: real,
    graminoid_wtp_limit: real,
    inund_stress: real)
    : (PhotosynthesisResult, real) =
    -- DESCRIPTION
    -- Calculation of net C-assimilation under water-stressed conditions
    -- (demand>supply see function canopy_exchange). Utilises a numerical
    -- iteration procedure to find the level of stomatal aperture (characterised by
    -- lambda, the ratio of leaf intercellular to ambient CO2 concentration) which
    -- satisfies simulataneously a canopy-conductance based and light-based
    -- formulation of photosynthesis (Eqns 2, 18 and 19, Haxeltine & Prentice (1996)).
    -- Numerical method is a tailored implementation of the bisection method,
    -- assuming root (f(lambda)=0) bracketed by f(0.02)<0 and
    -- f(lambda_max)>0 (Press et al 1986)
    -- The bisection method terminates when we're close enough to a root
    -- (absolute value of f(lambda) < EPS), or after a maximum number of
    -- iterations.
    -- Note that the function sometimes doesn't search for a lambda,
    -- and returns zero assimilation (for instance if there is no
    -- root within the valid interval, or if daylength is zero).
    -- So if zero assimilation is returned, the returned lambda should
    -- not be used!
    -- OUTPUT PARAMETER
    -- phot_result = result of photosynthesis for the found lambda
    -- lambda      = the lambda found by the bisection method (see above)
    -- Set lambda to something for cases where we don't actually search for
    -- a proper lambda. This value shouldn't be used (see documentation
    -- above), but we'll set it to something anyway so we don't return
    -- random garbage.
  let lambda = -1
  in if (negligible(fpc) || negligible(fpar) || negligible(gcbase * daylength * 3600)) -- Return zero assimilation
  then (PhotosynthesisResult(), lambda) else
  -- Canopy conductance component associated with photosynthesis on a
  -- daily basis (mm / m2 / day)
  let gcphot : real = gcbase * daylength * 3600 / 1.6 * co2 * CO2_CONV
  -- At this point the function f(x) = g(x) - h(x) can be calculated as:
  --
  -- g(x) = phot_result.adtmm / fpc (after a call to photosynthesis with lambda x)
  -- h(x) = gcphot * (1 - x)
  -- Evaluate f(lambda_max) to see if there's a root
  -- in the interval we're searching
  let ps_env = {co2=co2, temp=temp, par=par, fpar=fpar, daylength=daylength}
  let ps_stress = {ifnlimvmax=ifnlimvmax, moss_ps_limit=moss_wtp_limit, graminoid_ps_limit=graminoid_wtp_limit, inund_stress=inund_stress}
  let phot_result = photosynthesis(ps_env, ps_stress, pft, pft.lambda_max, nactive, vmax)
  let f_lambda_max : real = phot_result.adtmm / fpc - gcphot * (1 - pft.lambda_max)
  in if (f_lambda_max <= 0) -- Return zero assimilation
  then (PhotosynthesisResult(), lambda) else
  let EPS : real = 0.1 -- minimum precision of solution in bisection method
  -- Implement numerical solution
  let x1 : real = 0.02                      -- minimum bracket of root
  let x2 : real = pft.lambda_max            -- maximum bracket of root
  let rtbis : real = x1                     -- root of the bisection
  let dx : real = x2 - x1
  let MAXTRIES : int = 6 -- maximum number of iterations towards a solution
  let b : int = 0        -- number of tries so far towards solution
  let fmid : real = EPS + 1.0
  let (_, _, xmid, _, _, _, _, phot_result) =
  loop (b, dx, _, rtbis, fmid, ps_env, ps_stress, _)
      = (b, dx, 0.0, rtbis, fmid, ps_env, ps_stress, phot_result)
    while (abs(fmid) > EPS && b <= MAXTRIES) do
      let dx = dx * 0.5
      let xmid = rtbis + dx -- current guess for lambda
      -- Call function photosynthesis to calculate alternative value
      -- for total daytime photosynthesis according to Eqns 2 & 19,
      -- Haxeltine & Prentice (1996), and current guess for lambda
      let phot_result = photosynthesis(ps_env, ps_stress, pft, xmid, nactive, vmax)
     -- Evaluate fmid at the point lambda=xmid
     -- fmid will be an increasing function of xmid, with a solution
     -- (fmid=0) between x1 and x2
     -- Second term is total daytime photosynthesis (mm/m2/day) implied by
     -- canopy conductance and current guess for lambda (xmid)
     -- Eqn 18, Haxeltine & Prentice 1996
      let fmid = phot_result.adtmm / fpc - gcphot * (1 - xmid)
      let rtbis = if (fmid < 0) then xmid else rtbis
      in (b + 1, dx, xmid, rtbis, fmid, ps_env, ps_stress, phot_result)
  -- bvoc
  let lambda = xmid
  in (phot_result, lambda)
---------------------------------------------------------------------------------------
-- AUTOTROPHIC RESPIRATION
-- Internal function (do not call directly from framework)
let respiration( gtemp_air : real,  gtemp_soil : real,  lifeform : lifeformtype,
   respcoeff : real,  cton_sap : real,  cton_root : real,
   cmass_sap : real,  cmass_root_today : real,  assim : real) : real =
  -- DESCRIPTION
  -- Calculation of daily maintenance and growth respiration for individual with
  -- specified life form, phenological state, tissue C:N ratios and daily net
  -- assimilation, given current air and soil temperatures.
  -- Sitch et al. (2000), Lloyd & Taylor (1994), Sprugel et al (1996).
  -- NOTE: leaf respiration is not calculated here, but included in the calculation
  -- of net assimilation (function production above) as a proportion of rubisco
  -- capacity (Vmax).
  -- INPUT PARAMETERS
  -- gtemp_air  = respiration temperature response incorporating damping of Q10
  --              response due to temperature acclimation (Eqn 11, Lloyd & Taylor
  --              1994) Eqn B2 below
  -- gtemp_soil = as gtemp_air given soil temperature
  -- lifeform   = PFT life form class (TREE or GRASS)
  -- respcoeff  = PFT respiration coefficient
  -- cton_sap   = PFT sapwood C:N ratio
  -- cton_root  = PFT root C:N ratio
  -- phen       = vegetation phenological state (fraction of potential leaf cover)
  -- cmass_sap  = sapwood C biomass on grid cell area basis (kgC/m2)
  -- cmass_root = fine root C biomass on grid cell area basis (kgC/m2)
  -- assim      = net assimilation on grid cell area basis (kgC/m2/day)
  -- OUTPUT PARAMETER
  -- resp       = sum of maintenance and growth respiration on grid cell area basis
  --              (kgC/m2/day)
  -- guess2008 - following a comment by Annett Wolf, the following parameter value was changed:
  -- const double K=0.0548 -- OLD value
  let K : real = 0.095218 in  -- NEW parameter value in respiration equations
  -- See the comment after Eqn (4) below.
  --double resp_sap    -- sapwood respiration (kg/m2/day)
  --double resp_root   -- root respiration (kg/m2/day)
  --double resp_growth -- growth respiration (kg/m2/day)
  -- Calculation of maintenance respiration components for each living tissue:
  --
  -- Based on the relations
  --
  -- (A) Tissue respiration response to temperature
  --     (Sprugel et al. 1996, Eqn 7)forest_floor
  --
  --     (A1) Rm = 7.4e-7 * N * f(T)
  --     (A2) f(T) = EXP (beta * T)
  --
  --       where Rm   = tissue maintenance respiration rate in mol C/sec
  --             N    = tissue nitrogen in mol N
  --             f(T) = temperature response function
  --             beta = ln Q10 / 10
  --             Q10  = change in respiration rate with a 10 K change
  --                    in temperature
  --             T    = tissue absolute temperature in K
  --
  -- (B) Temperature response of soil respiration across ecosystems
  --     incorporating damping of Q10 response due to temperature acclimation
  --     (Lloyd & Taylor 1994, Eqn 11)
  --
  --     (B1) R = R10 * g(T)
  --     (B2) g(T) = EXP [308.56 * (1 / 56.02 - 1 / (T - 227.13))]
  --
  --       where R    = respiration rate
  --             R10  = respiration rate at 10 K
  --             g(T) = temperature response function at 10 deg C
  --             T    = soil absolute temperature in K
  --
  -- Mathematical derivation:
  --
  -- For a tissue with C:N mass ratio cton, and C mass, c_mass, N concentration
  -- in mol given by
  --  (1) N = c_mass / cton / atomic_mass_N
  -- Tissue respiration in gC/day given by
  --  (2) R = Rm * atomic_mass_C * seconds_per_day
  -- From (A1), (1) and (2),
  --  (3) R = 7.4e-7 * c_mass / cton / atomic_mass_N * atomic_mass_C
  --         ---seconds_per_day * f(T)
  -- Let
  --  (4) k = 7.4e-7 * atomic_mass_C / atomic_mass_N * seconds_per_day
  --        = 0.0548
  -- guess2008 - there is an ERROR here, spotted by Annett Wolf
  -- If we calculate the respiration at 20 degC using g(T) and compare it to
  -- Sprugel's eqn 3, for 1 mole tissue N, say, we do NOT get the same result with this
  -- k value. This is because g(T) = 1 at 10 degC, not 20 degC. Changing k from 0.0548
  -- to 0.095218 gives exactly the same results as Sprugel at 20 degC. The scaling factor
  -- 7.4e-7 used here is taken from Sprugel's eqn. (7), but they used f(T), not g(T), and
  -- these are defined on different bases.
  -- from (3), (4)
  --  (5) R = k * c_mass / cton * f(T)
  -- substituting ecosystem temperature response function g(T) for f(T) (Eqn B2),
  --  (6) R = k * c_mass / cton * g(T)
  -- incorporate PFT-specific respiration coefficient to model acclimation
  -- of respiration rates to average (temperature) conditions for PFT (Ryan 1991)
  --  (7) R_pft = respcoeff_pft * k * c_mass / cton * g(T)
  if (lifeform == TREE) then
    -- Sapwood respiration (Eqn 7)
    let resp_sap = respcoeff * K * cmass_sap / cton_sap * gtemp_air
    -- Root respiration (Eqn 7)
    -- Assumed that root phenology follows leaf phenology
    let resp_root = respcoeff * K * cmass_root_today / cton_root * gtemp_soil
    -- Growth respiration = 0.25 ( GPP - maintenance respiration)
    let resp_growth = (assim - resp_sap - resp_root) * 0.25
    -- guess2008 - disallow negative growth respiration
    -- (following a comment (060823) from Annett Wolf)
    let resp_growth = if (resp_growth < 0.0) then 0.0 else resp_growth
    -- Total respiration is sum of maintenance and growth respiration
    in resp_sap + resp_root + resp_growth
  else if (lifeform == GRASS || lifeform == MOSS) then
    -- Root respiration
    let resp_root = respcoeff * K * cmass_root_today / cton_root * gtemp_soil
    -- Growth respiration (see above)
    let resp_growth = (assim - resp_root) * 0.25
    -- guess2008 - disallow negative growth respiration
    -- (following a comment (060823) from Annett Wolf)
    let resp_growth = if (resp_growth < 0.0) then 0.0 else resp_growth
    -- Total respiration (see above)
    in resp_root + resp_growth
  else nan --else fail ("Canopy exchange function respiration: unknown life form")
--
--/ Net Primary Productivity
--- Includes BVOC calculations \see bvoc.cpp
--/
let npp( patch : Patch
       , climate : Climate
       , vegetation : [npft]Individual
       , day : Day
       , pfts : [npft]Pft
       , spfts : [npft]Standpft)
       : (Patch, Vegetation) =
  -- Determination of daily NPP. Leaf level net assimilation calculated for non-
  -- water-stressed individuals (i.e. with fully-open stomata) using base value
  -- from function demand (above) for water-stressed individuals using base value
  -- for canopy conductance by a simultaneous solution of light-based and canopy
  -- conductance-based equations for net daily photosynthesis (see function
  -- assimilation wstress above). The latter uses the PFT-specific base value for
  -- conductance from function aet_water_stress (above).
  -- Plant respiration obtained by a call to function respiration (above).
  let (par, temp, hours, rad, gtemp) =
    if (date.diurnal()) then
      (climate.pars[day.period]
      ,climate.temps[day.period]
      ,24-- diurnal "daylength" to convert to daily units
      ,climate.rads[day.period]
      ,climate.gtemps[day.period])
    else
      (climate.par
      ,climate.temp
      ,climate.daylength
      ,climate.rad
      ,climate.gtemp)

  let (vegetation, argpile) = unzip <| map (\indiv ->
    -- For this individual ...
    let pft : Pft = pfts[indiv.pft_id]
    let ppft : Patchpft = ppfts[pft.pft_id]
    --Don't do calculations for crops outside their growingseason
    in if (!indiv.growingseason()) then (indiv with dnpp = 0.0, (false, false, 0, 0, 0, 0.0, 0.0, 0.0, date, 0))
    else
      let pftco2 : real = get_co2(patch, climate, pft, stand, soil, gridcell)
      let inund_stress : real = get_inund_stress(ppft, stand, gridcell)
      let graminoid_wtp_limit : real = get_graminoid_wtp_limit(patch, pft, stand, soil, gridcell)
      let moss_wtp_limit : real = get_moss_wtp_limit(pft, stand, soil, gridcell)
      let phot : PhotosynthesisResult = if date.diurnal() then indiv.phots[day.period] else indiv.photosynthesis_result
      let (phot, lambda) =
        if (indiv.wstress) then
          -- Water stress - derive assimilation by simultaneous solution
          -- of light- and conductance-based equations of photosynthesis
          assimilation_wstress(pft, pftco2, temp, par, hours, indiv.fpar, indiv.fpc,
            ppft.gcbase, phot.vm, phot, lambda,
            indiv.nactive / indiv.nextin, ifnlim, moss_wtp_limit, graminoid_wtp_limit, inund_stress)
        else (phot, lambda)
      let assim = phot.net_assimilation()
      --if (ifbvoc) { --TODO FIXME: this is false in global.ins, so it is omitted here for now
      --  PhotosynthesisResult phot_nostress = date.diurnal() ? indiv.phots[day.period] : indiv.photosynthesis_result
      --  bvoc(temp, hours, rad, climate, patch, indiv, pft, phot_nostress, phot.adtmm, day)
      --}
      -- Calculate autotrophic respiration
      let cmass_root =
        if (indiv.cropindiv && indiv.cropindiv.isintercropgrass && indiv.phen == 0.0) then
          0.0
        else
          cmass_root_today(indiv, pft, patchpft)
      -- Static root and sap wood C:N ratio if no N limitation
      -- to not let N affect respiration for C only version of model
      let (cton_sap, cton_root) =
        if (ifnlim) then
          (indiv.cton_sap(), indiv.cton_root())
        else
          (pft.cton_sap_avr, pft.cton_root_avr)

      let (X) = respiration(gtemp, soil.gtemp, pfts[indiv.pft_id].lifeform,
                            pfts[indiv.pft_id].respcoeff, cton_sap, cton_root,
                            indiv.cmass_sap, cmass_root, assim, resp)
      -- Convert to averages for this period for accounting purposes
      let assim = assim / date.subdaily
      let resp = resp / date.subdaily
      -- Update accumulated annual NPP and daily vegetation-atmosphere flux
      let indiv = indiv with dnpp = assim - resp
      let indiv = indiv with anpp = indiv.anpp + indiv.dnpp

      let alive : bool = indiv.alive
      let itoicg : bool = istruecrop_or_intercropgrass indiv pft
      let pft_id = indiv.pft_id

      -- report fluxes reduction extracted from map, see below

      let indiv =
        if (lai_today(indiv, pfts[indiv.pft_id], ppfts[indiv.pft_id]) > indiv.mlai_max[date.month]) then
          indiv with mlai_max = indiv.mlai_max with [date.month] = lai_today(indiv, pfts[indiv.pft_id], ppfts[indiv.pft_id])
        else indiv

      let indiv =
        if (day.isend) then
          indiv with mlai = mlai with [date.month] = indiv.mlai[date.month] + lai_today(indiv, pfts[indiv.pft_id], ppfts[indiv.pft_id]) / realFromInt date.ndaymonth[date.month]
        else indiv
      in (indiv, (alive, itoicg, NPP, GPP, RA, indiv.dnpp, assim, resp, date, pft_id))
    ) vegetation

  let patch =
    reduce patch (\patch (alive, itoicg, NPP, GPP, RA, dnpp, assim, resp, date, pft_id) ->
      let fluxes = patch.fluxes
      let fluxes = Individual_report_flux_PerPFTFluxType(fluxes, alive, itoicg, NPP, indiv.dnpp, date, pft_id)
      let fluxes = Individual_report_flux_PerPFTFluxType(fluxes, alive, itoicg, GPP, assim, date, pft_id)
      let fluxes = Individual_report_flux_PerPFTFluxType(fluxes, alive, itoicg, RA, resp, date, pft_id)
      let patch = patch with fluxes = fluxes
      in patch
    )  argpile

  in (patch, vegetation)
--/ Forest-floor conditions
--- Called in cohort/individual mode (not population mode) to quantify growth
--  conditions at the forest floor for each PFT
--  Calculates net assimilation at top of grass canopy (or at soil surface if
--  there is none).
--/
let forest_floor_conditions(stand: Stand, patch : Patch, climate : Climate, spfts : [npft]Standpft, ppfts : [npft]Patchpft, pfts : [npft]Pft)
  : ([npft]Patchpft, [npft]Standpft) =
  --double lambda      -- not used here
  let phot : PhotosynthesisResult = PhotosynthesisResult()
  let (ppfts, spfts, p) =
  loop (ppfts, spfts, p)
    = (ppfts, spfts, 0)
  while (p<npft) do
    let ppft : Patchpft = ppfts[p]
    let spft : Standpft = spfts[p]
    in if (!spft.active) then (ppfts, spfts, p+1) else
      let pft : Pft = pfts[spft.pft_id]
      -- peatland limits on photosynthesis
      let pftco2 : real = get_co2(patch, climate, pft, stand, soil, gridcell) -- was patch.stand.get_gridcell().climate
      let inund_stress : real = get_inund_stress(ppft, stand, gridcell)
      let graminoid_wtp_limit : real = get_graminoid_wtp_limit(patch, pft, stand, soil, gridcell)
      let moss_wtp_limit : real = get_moss_wtp_limit(pft, stand, soil, gridcell)
      -- Initialise net photosynthesis sum on first day of year
      let ppft =
        if (date.day == 0) then
          ppft with anetps_ff = 0.0
        else ppft

      let ppft =
        if (patch.stand.landcover != CROPLAND || pft.phenology != CROPGREEN && ppft.cropphen.growingseason) then
          let assim =
            if (ppft.wstress_day) then
              --PhotosynthesisStresses ps_stress
              --ps_stress.set(false, get_moss_wtp_limit(pft, stand, soil, gridcell), get_graminoid_wtp_limit(patch, pft, stand, soil, gridcell), get_inund_stress(ppft, stand, gridcell))
              let (phot_result, lambda) =
                                assimilation_wstress(pft, pftco2, climate.temp, climate.par,
                                  climate.daylength, patch.fpar_grass * ppft.phen, 1.0, ppft.gcbase_day,
                                  spft.photosynthesis.vm, phot, lambda, 1.0, false,
                                  moss_wtp_limit, graminoid_wtp_limit, inund_stress)
              in phot.net_assimilation()
            else spft.photosynthesis.net_assimilation() * ppft.phen * patch.fpar_grass
          -- Accumulate annual value
          in ppft with anetps_ff = ppft.anetps_ff + assim
        else ppft

      let (ppft, spft) =
        if (ppfts, spfts)(date.islastmonth && date.islastday) then
          -- Avoid negative ppft.anetps_ff
          let ppft = ppft with anetps_ff = max(0.0, ppft.anetps_ff)
          let spft =
            if (ppft.anetps_ff > spft.anetps_ff_max) then
              spft with anetps_ff_max = ppft.anetps_ff
            else spft
          in (ppft, spft)
        else (ppft, spft)
    in (ppfts with [p] = ppft, spfts with [p] = spft, p+1)
  in (ppfts, spfts)

--/ Leaf senescence for crops Eqs. 8,9,13 and 14 in Olin 2015
let leaf_senescence(stand: Stand, patch: Patch, vegetation: [npft]Individual, ppfts: [npft]Patchpft) : [npft]Individual =
  -- FIXME when you know returntype if (!(stand.is_true_crop_stand() && ifnlim)) { return }

  let vegetation = map (\indiv ->
    -- Age dependent N retranslocation, Sec. 2.1.3 Olin 2015
    let indiv =
      if (ppfts[indiv.pft_id].cropphen.dev_stage > 1.0) then
        let senNr = 0.1
        let senN = senNr * (indiv.nmass_leaf-cmass_leaf_today(indiv, pft, patchpft) / (pfts[indiv.pft_id].cton_leaf_max))
        in
        -- Senescence is not done during spinup
        if (date.year > nyear_spinup && senN > 0) then
          let indiv = indiv with nmass_leaf = indiv.nmass_leaf - senN
          let indiv = indiv with cropindiv = indiv.cropindiv with nmass_agpool = indiv.cropindiv.nmass_agpool + senN
          in indiv
        else indiv
      else indiv

    let r =
    -- N dependant C mass loss, with an inertia of 1/10, Eq. 13 Olin 2015
      if (cmass_leaf_today(indiv, pft, patchpft) > 0.0) then
        let Ln : real = indiv.lai_nitrogen_today()
        let Lnld : real = lai_today(indiv, pfts[indiv.pft_id], ppfts[indiv.pft_id])
        in (Lnld - min(Lnld, Ln))/pfts[indiv.pft_id].sla/10.0
      else 0.0
    -- No senescence during the initial growing period
    let indiv = indiv with daily_cmass_leafloss =
      if (ppfts[indiv.pft_id].cropphen.fphu < 0.05) then
        0.0
      else
        max(0.0, r)
    let indiv = indiv with daily_nmass_leafloss = 0.0
    in indiv
    ) vegetation
  in vegetation

--/ Initiate required variables for the module
let init_canexch(patch : Patch, climate : Climate, vegetation : [npft]Individual, date : Date) : (Patch, Vegetation) =
  let vegetation =
    if (date.day == 0) then
      map (\indiv ->
        let indiv = indiv with anpp = 0.0
        let indiv = indiv with leafndemand = 0.0
        let indiv = indiv with rootndemand = 0.0
        let indiv = indiv with sapndemand = 0.0
        let indiv = indiv with storendemand = 0.0
        let indiv = indiv with hondemand = 0.0
        let indiv = indiv with nday_leafon = 0
        let indiv = indiv with avmaxnlim = 1.0
        let indiv = indiv with cton_leaf_aavr = 0.0
        in
        if (!negligible(indiv.cmass_leaf) && !negligible(indiv.nmass_leaf)) then
          indiv with cton_leaf_aopt = indiv.cmass_leaf / indiv.nmass_leaf
        else
          let indiv = indiv with cton_leaf_aopt = pfts[indiv.pft_id].cton_leaf_max
          let indiv = indiv with mlai = replicate 12 realzero
          let indiv = indiv with mlai_max = replicate 12 realzero
          in indiv
      ) vegetation
    else vegetation
  in (patch with wdemand_day = 0, vegetation)

--let canopy_exchange(patch : Patch, climate : Climate, individuals : [][][]Individuals) : Patch =
--  let vegetation : []Individual = individuals[patch.gridcell_id, patch.stand_id, patch.patch_id]
--  in patch
--- Canopy exchange
--- Vegetation-atmosphere exchange of CO2 and water including calculations
--  of actual evapotranspiration (AET), canopy conductance, carbon assimilation
--  and autotrophic respiration.
--  Should be called each simulation day for each modelled area or patch,
--  following update of leaf phenology and soil temperature and prior to update
--  of soil water.
---
let canopy_exchange(patch : Patch, vegetation : [npft]Individual, climate : Climate, pfts : [npft]Pft, date : Date, spfts: [npft]Standpft, npfts: [npft]Patchpft) : (Patch, Vegetation, [npft]Standpft, [npft]Patchpft) =
  -- NEW ASSUMPTIONS CONCERNING FPC AND FPAR (Ben Smith 2002-02-20)
  -- FPAR = average individual fraction of PAR absorbed on patch basis today,
  --        including effect of current leaf phenology (this differs from previous
  --        versions of LPJ-GUESS in which FPAR was on an FPC basis)
  -- FPC =  PFT population (population mode), cohort (cohort mode) or individual
  --        (individual mode) fractional projective cover as a fraction of patch area
  --        (in population mode, corresponds to LPJF variable fpc_grid). Updated
  --        annually based on leaf-out LAI (see function allometry in growth module).
  --        (FPC was previously equal to summed crown area as a fraction of patch
  --        area in cohort/individual mode)
  -- Retrieve Vegetation and Climate objects for this patch
  -- Initial no-stress canopy exchange processes
  let (patch, vegetation) = init_canexch(patch, climate, vegetation)
  -- Canopy exchange processes
  let (patch, vegetation) = fpar(patch, vegetation, climate, pfts, date)
  -- Calculates no-stress daily values of photosynthesis and gpterm
  let (vegetation, spfts, ppfts) = photosynthesis_nostress(patch, climate, spfts ,ppfts, pfts)
  -- Nitrogen demand
  let (vegetation, patch) = ndemand(patch, vegetation, gridcell , soil , pfts)
  -- Nitrogen stress
  let vegetation = vmax_nitrogen_stress(patch, climate, vegetation)
  -- Only these processes are affected in diurnal mode
  let day = Day(date)

  let (patch, vegetation, ppfts, pfts, _) =
  loop (patch, vegetation, ppfts, pfts, day) =
       (patch, vegetation, ppfts, pfts, Day(date))
  while(day.period != date.subdaily) do
    let (patch, vegetation) = wdemand(patch, climate, vegetation, day)
    let (vegetation, ppfts) = aet_water_stress(patch, stand, vegetation, day, pfts, spfts, ppfts)

    let pfts = water_scalar(stand, patch, vegetation, day, ppfts, pfts)
    let (patch, vegetation) = npp(patch, climate, vegetation, day, pfts, spfts)
    in (patch, vegetation, ppfts, pfts, day.next())

  let vegetation = leaf_senescence(stand, patch, vegetation, ppfts)
  -- Forest-floor conditions
  let (ppfts, spfts) = forest_floor_conditions(stand, patch, climate, spfts, ppfts, pfts)
  -- Total potential evapotranspiration for patch (mm, patch basis)
  -- is a sum of: (1) potential transpirative demand of the vegetation
  -- (2) evaporation of canopy-intercepted precipitation and
  -- (3) evaporation from the soil surface of non-vegetated parts of patch -
  --     currently with boleht at 0, a patch has only two surfaces - vegetated
  --     and non-vegetated.
  -- This value is only diagnostic, it is not to be used in further calculations.
  -- Correct value should use daily value of patch.demand_leafon.
  let pet_patch : real = patch.wdemand_day * patch.fpc_total + patch.intercep +
        climate.eet * PRIESTLEY_TAYLOR * max(1.0-patch.fpc_total, 0.0)
  let patch = patch with apet = patch.apet + pet_patch
  let patch = patch with mpet = patch.mpet with [date.month] = patch.mpet[date.month] + pet_patch
  in (patch, vegetation, spfts, ppfts)

--type GlobalData = {
--    pfts : [npft]Pft,
--    standtypes: []StandType,
--    managementtypes : []ManagementType,
--    soultypes: []Soiltype,
    --gridcells : []Gridcell
--    climate: []Climate,
--    weathergenstate: []WeatherGenState,
--    landcover: []Landcover,
--    massbalance: []MassBalance,
  --gridcellsts: [][]Gridcellst,
  --    gridcellpfts: [][]Gridcellpft,
--    stands: [][]Stand,
--    standpfts : [][][]Standpft,
--    patches : [][][]Patch,
--    patchpfts : [][][][]Patchpft,
--    individuals : [][][][]Individual
--}
