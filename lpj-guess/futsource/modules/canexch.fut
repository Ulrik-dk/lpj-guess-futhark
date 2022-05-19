-------------------------------------Combined imports-----------------------------------
open import "../framework/guess"

--#include "config.h"
--#include "canexch.h"
--#include "driver.h"
open import "../modules/q10"
open import "../futhark-extras"
--#include "bvoc.h"
--#include "ncompete.h"
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
let lambertbeer (lai : real) = exp(-0.5 * lai)



-------------------------------------PART FROM .cpp------------------------------------



-- leaf nitrogen (kgN/kgC) not associated with photosynthesis
-- (value given by Haxeltine & Prentice 1996b) */
let N0 : real = 7.15 * G_PER_MG

-- Lookup tables for parameters with Q10 temperature responses

-- lookup table for Q10 temperature response of CO2/O2 specificity ratio
let lookup_tau = LookupQ10(0.57, 2600.0)

-- lookup table for Q10 temperature response of Michaelis constant for O2
let lookup_ko = LookupQ10(1.2, 3.0e4)

-- lookup table for Q10 temperature response of Michaelis constant for CO2
let lookup_kc = LookupQ10(2.1, 30.0)



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
--/

let photosynthesis_ugly_naive_dont_use_this_one(ps_env : PhotosynthesisEnvironment,
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

  -- No photosynthesis during polar night, outside of temperature range or no RuBisCO activity
  in if (negligible(daylength) || negligible(fpar) || temp > pft.pstemp_max || temp < pft.pstemp_min || !(boolFromReal vm)) then PhotosynthesisResult() else

  -- Scale fractional PAR absorption at plant projective area level (FPAR) to
  -- fractional absorption at leaf level (APAR)
  -- Eqn 4, Haxeltine & Prentice 1996a
  let apar : real = par * fpar * alphaa(pft)

  -- Calculate temperature-inhibition coefficient
  -- This function (tscal) is mathematically identical to function tstress in LPJF.
  -- In contrast to earlier versions of modular LPJ and LPJ-GUESS, it includes both
  -- high- and low-temperature inhibition.
  let k1 : real = (pft.pstemp_min+pft.pstemp_low) / 2.0
  let tscal : real = (1.0 - 0.01*exp(4.6/(pft.pstemp_max-pft.pstemp_high)*(temp-pft.pstemp_high)))/
                    (1.0+exp((k1-temp)/(k1-pft.pstemp_min)*4.6))

  let (b, c1, c2) : (real, real, real) =
    if (pft.pathway == C3) then      -- C3 photosynthesis
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

      let b : real = if pft.lifeform == MOSS then BC_moss else BC3
      -- see Wania et al. 2009b
      in (b, c1, c2)

    else               -- C4 photosynthesis
      -- Calculation of C1_C4 given actual pi (lambda)
      -- C1_C4 incorporates term accounting for effect of intercellular CO2
      -- concentration on photosynthesis (Eqn 14, 16, Haxeltine & Prentice 1996a)
      let b = BC4
      let c1 = min(lambda/LAMBDA_SC4, 1.0) * ALPHA_C4
      let c2 = 1.0
      in (b, c1, c2)

  -- Calculation of non-water-stressed rubisco capacity (Eqn 11, Haxeltine & Prentice 1996a)
  let ps_result = PhotosynthesisResult()


  let ps_result =
    if (vm < 0) then
      let (vm, vmaxnlim, nactive_opt) = vmax(b, c1, c2, apar, tscal, daylength, temp, nactive, ifnlimvmax)
      let ps_result = ps_result with vm = vm
      let ps_result = ps_result with vmaxnlim = vmaxnlim
      let ps_result = ps_result with nactive_opt = nactive_opt
      in ps_result
    else ps_result

  -- Calculation of daily leaf respiration
  -- Eqn 10, Haxeltine & Prentice 1996a
  let ps_result = ps_result with rd_g = (ps_result.vm * b)

  -- PAR-limited photosynthesis rate (gC/m2/h)
  -- Eqn 3, Haxeltine & Prentice 1996a
  let ps_result = ps_result with je = (c1 * tscal * apar * CMASS * CQ / daylength)

  -- Rubisco-activity limited photosynthesis rate (gC/m2/h)
  -- Eqn 5, Haxeltine & Prentice 1996a
  let jc : real = c2 * ps_result.vm / 24.0

  -- Calculation of daily gross photosynthesis
  -- Eqn 2, Haxeltine & Prentice 1996a
  -- Notes: - there is an error in Eqn 2, Haxeltine & Prentice 1996a (missing
  --       theta in 4*theta*je*jc term) which is fixed here
  let ps_result = ps_result with agd_g = ((ps_result.je + jc - sqrt((ps_result.je + jc) * (ps_result.je + jc) - 4.0 * THETA * ps_result.je * jc)) /(2.0 * THETA) * daylength)

  let ps_result = if (!iftwolayersoil) then
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
    let ps_result = ps_result with agd_g = (ps_result.agd_g * inund_stress)
    let ps_result = ps_result with rd_g = (ps_result.rd_g * inund_stress)

    -- 2a) Moss dessication
    let ps_result =
      if (pft.lifeform == MOSS) then
        -- Reduce agd_g using moss_wtp_limit (lies between [0.3, 1.0])
        let ps_result = ps_result with agd_g = (ps_result.agd_g * moss_ps_limit)
        let ps_result = ps_result with rd_g = (ps_result.rd_g * moss_ps_limit)
        in ps_result
      else ps_result

    -- 2b) Graminoid dessication (NB! all other PFTs have graminoid_ps_limit == 1)
    let ps_result = ps_result with agd_g = ps_result.agd_g * graminoid_ps_limit
    let ps_result = ps_result with rd_g = ps_result.rd_g * graminoid_ps_limit
    in ps_result
  else ps_result

    -- Leaf-level net daytime photosynthesis (gC/m2/day)
  -- Based on Eqn 19, Haxeltine & Prentice 1996a
  let adt : real = ps_result.agd_g - daylength / 24.0 * ps_result.rd_g

  -- Convert to CO2 diffusion units (mm/m2/day) using ideal gas law
  let ps_result = ps_result with adtmm = (adt / CMASS * 8.314 * (temp + K2degC) / PATMOS * 1e3)
  in ps_result

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
  let temp : real = ps_env.temp -- 21.794977
  let co2 : real = ps_env.co2 -- 296.378500
  let fpar : real = ps_env.fpar --1.000000
  let par : real = ps_env.par --7307193.574951
  let daylength : real = ps_env.daylength --10.826868

  -- Get the stresses
  let ifnlimvmax : bool = ps_stresses.ifnlimvmax -- false

  let pstemp_max = pft.pstemp_max --55
  let pstemp_high = pft.pstemp_high--30
  let pstemp_low = pft.pstemp_low--25
  let pstemp_min = pft.pstemp_min--2

  let pathway = pft.pathway
  let lifeform = pft.lifeform

  -- No photosynthesis during polar night, outside of temperature range or no RuBisCO activity
  in if (negligible(daylength) || negligible(fpar) || temp > pstemp_max || temp < pstemp_min || !(boolFromReal vm)) then PhotosynthesisResult() else

  -- Scale fractional PAR absorption at plant projective area level (FPAR) to
  -- fractional absorption at leaf level (APAR)
  -- Eqn 4, Haxeltine & Prentice 1996a

  let (par, fpar) = (7.30719e+06, 1)

  let apar : real = 4.74968e+6--par * fpar * alphaa(pft)

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
    else (ps_result.vm, ps_result.vmaxnlim, ps_result.nactive_opt)


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

    -- 2a) Moss dessication
    -- Reduce agd_g using moss_wtp_limit (lies between [0.3, 1.0])
    let agd_g = if lifeform == MOSS then agd_g * moss_ps_limit else agd_g
    let rd_g = if lifeform == MOSS then rd_g * moss_ps_limit else rd_g


    -- 2b) Graminoid dessication (NB! all other PFTs have graminoid_ps_limit == 1)
    let agd_g = agd_g * graminoid_ps_limit
    let rd_g = rd_g * graminoid_ps_limit
    in (agd_g, rd_g)
  else (agd_g, rd_g)

    -- Leaf-level net daytime photosynthesis (gC/m2/day)
  -- Based on Eqn 19, Haxeltine & Prentice 1996a
  let adt : real = agd_g - daylength / 24.0 * rd_g

  -- Convert to CO2 diffusion units (mm/m2/day) using ideal gas law
  let adtmm = (adt / 12 * 8.314 * (temp + 273.15) / 100000 * 1000) -- expected 20, got 31 -- 0.64
  in {agd_g=agd_g, -- expected 11, got 16 -- 0.67
      adtmm=adtmm, -- expected 20, got 31 -- 0.64
      je=je, -- expected 1.248, got 1.711 -- 0.73
      rd_g=rd_g, -- expected 2.704, got 2.364 -- 1.144
      vm=vm, -- expected 180.282, got 118.218 -- 1.525
      nactive_opt=nactive_opt, -- expected 1.2033e-2, got 7.890385969566833e-3 --1.525
      vmaxnlim=vmaxnlim
  }


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
    : (Pft, PhotosynthesisResult, real) =

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
  then (pft, PhotosynthesisResult(), lambda) else
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
    then (pft, PhotosynthesisResult(), lambda) else
      let EPS : real = 0.1 -- minimum precision of solution in bisection method

      let xmid : real = 0.0 -- TODO

      -- Implement numerical solution
      let x1 : real = 0.02                      -- minimum bracket of root
      let x2 : real = pft.lambda_max            -- maximum bracket of root
      let rtbis : real = x1                     -- root of the bisection
      let dx : real = x2 - x1

      let MAXTRIES : int = 6 -- maximum number of iterations towards a solution
      let b : int = 0        -- number of tries so far towards solution

      let fmid : real = EPS + 1.0

      let ps_env = {co2=co2, temp=temp, par=par, fpar=fpar, daylength=daylength}
      let (_, _, xmid, _, _, _, _, pft, phot_result) =
      loop (b, dx, _, rtbis, fmid, ps_env, ps_stress, pft, _)
          = (b, dx, xmid, rtbis, fmid, ps_env, ps_stress, pft, phot_result)
        while (abs(fmid) > EPS && b <= MAXTRIES) do
          let dx = dx * 0.5
          let xmid = rtbis + dx -- current guess for lambda

          -- Call function photosynthesis to calculate alternative value
          -- for total daytime photosynthesis according to Eqns 2 & 19,
          -- Haxeltine & Prentice (1996), and current guess for lambda

          let phot_result = photosynthesis(ps_env, ps_stress, pft,
                                                      xmid, nactive, vmax)

         -- Evaluate fmid at the point lambda=xmid
         -- fmid will be an increasing function of xmid, with a solution
         -- (fmid=0) between x1 and x2

         -- Second term is total daytime photosynthesis (mm/m2/day) implied by
         -- canopy conductance and current guess for lambda (xmid)
         -- Eqn 18, Haxeltine & Prentice 1996

          let fmid = phot_result.adtmm / fpc - gcphot * (1 - xmid)
          let rtbis = if (fmid < 0) then xmid else rtbis
          in (b + 1, dx, xmid, rtbis, fmid, ps_env, ps_stress, pft, phot_result)
      -- bvoc
      let lambda = xmid
      in (pft, phot_result, lambda)
