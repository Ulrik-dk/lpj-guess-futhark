-------------------------------------Combined imports-----------------------------------
open import "../framework/guess"

--#include "config.h"
--#include "canexch.h"
--#include "driver.h"
--#include "q10.h"
--#include "bvoc.h"
--#include "ncompete.h"
--#include <assert.h>
-------------------------------------PART FROM .h--------------------------------------

-- Constants for photosynthesis calculations
-- conversion factor for solar radiation at 550 nm from J/m2 to mol_quanta/m2 (E=mol quanta) mol J-1
let CQ : f64 = 4.6e-6

-- intrinsic quantum efficiency of CO2 uptake, C3 plants
let ALPHA_C3 : f64 = 0.08

-- intrinsic quantum efficiency of CO2 uptake, C4 plants
let ALPHA_C4 : f64 = 0.053

-- O2 partial pressure (Pa)
let PO2 : f64 = 2.09e4

-- colimitation (shape) parameter
let THETA : f64 = 0.7

-- 'saturation' ratio of intercellular to ambient CO2 partial pressure for C4 plants
let LAMBDA_SC4 : f64 = 0.4

-- leaf respiration as fraction of maximum rubisco, C3 plants
let BC3 : f64 = 0.015

-- leaf respiration as fraction of maximum rubisco, C4 plants
let BC4 : f64 = 0.02

-- leaf respiration as fraction of maximum rubisco, mosses
-- see Wania et al. (2009b)
let BC_moss : f64 = 0.03

let CMASS : f64 = 12.0    -- atomic mass of carbon
let ALPHAA : f64 = 0.5    -- value chosen to give global carbon pool and flux values that
                -- agree with published estimates.
                -- scaling factor for PAR absorption from leaf to plant projective area level
                -- alias "twigloss". Should normally be in the range 0-1

let ALPHAA_NLIM : f64 = 0.65 -- Same as ALPHAA above but chosen to give pools and flux values
                 -- that agree with published estimates when Nitrogen limitation is
                 -- switched on.

let ALPHAA_CROP : f64 = 0.7      -- Value for crops without N limitation.
let ALPHAA_CROP_NLIM : f64 = 0.9  -- Value for crops with N limitation

-- Lambert-Beer extinction law (Prentice et al 1993 Monsi & Saeki 1953)
let lambertbeer (lai : f64) = exp(-0.5 * lai)

-------------------------------------PART FROM .cpp------------------------------------
--------------------------------------------------------------------------------------/
-- ASSIMILATION_WSTRESS
-- Internal function (do not call directly from framework)

let assimilation_wstress
      (pft: Pft, co2: f64, temp: f64, par: f64,
      daylength: f64, fpar: f64, fpc: f64, gcbase: f64,
      vmax: f64, phot_result: PhotosynthesisResult, lambda: f64,
      nactive: f64, ifnlimvmax: bool, moss_wtp_limit: f64, graminoid_wtp_limit: f64, inund_stress: f64)
      : (Pft, PhotosynthesisResult, f32)
      =

    -- DESCRIPTION
    -- Calculation of net C-assimilation under water-stressed conditions
    -- (demand>supply; see function canopy_exchange). Utilises a numerical
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
    let gcphot : f64 = gcbase * daylength * 3600 / 1.6 * co2 * CO2_CONV

    -- At this point the function f(x) = g(x) - h(x) can be calculated as:
    --
    -- g(x) = phot_result.adtmm / fpc (after a call to photosynthesis with lambda x)
    -- h(x) = gcphot * (1 - x)

    -- Evaluate f(lambda_max) to see if there's a root
    -- in the interval we're searching
    let ps_env = {co2=co2, temp=temp, par=par, fpar=fpar, daylength=daylength}
    let ps_stress = {ifnlimvmax=ifnlimvmax, moss_wtp_limit=moss_wtp_limit, graminoid_wtp_limit=graminoid_wtp_limit, inund_stress=inund_stress}
    let phot_result = photosynthesis(ps_env, ps_stress, pft, pft.lambda_max, nactive, vmax, phot_result)

    let f64 f_lambda_max = phot_result.adtmm / fpc - gcphot * (1 - pft.lambda_max) -- Return zero assimilation
    in if (f_lambda_max <= 0)
    then (pft, phot_result.clear(), lambda) else
      let f64 EPS = 0.1 -- minimum precision of solution in bisection method

      let f64 xmid = 0.0 -- TODO

      -- Implement numerical solution
      let f64 x1 = 0.02                      -- minimum bracket of root
      let f64 x2 = pft.lambda_max            -- maximum bracket of root
      let f64 rtbis = x1                     -- root of the bisection
      let f64 dx = x2 - x1

      let int MAXTRIES = 6 -- maximum number of iterations towards a solution
      let int b = 0        -- number of tries so far towards solution

      let f64 fmid = EPS + 1.0

      let ps_env = ps_env.set(co2, temp, par, fpar, daylength)
      in
      loop (b,dx,xmid,phot_result,fmid)
       while (abs(fmid) > EPS && b <= MAXTRIES) do
        let b = b++
        let dx = dx *= 0.5
        let xmid = rtbis + dx -- current guess for lambda

        -- Call function photosynthesis to calculate alternative value
        -- for total daytime photosynthesis according to Eqns 2 & 19,
        -- Haxeltine & Prentice (1996), and current guess for lambda

        let phot_result = photosynthesis(ps_env, ps_stress, pft,
                                          xmid, nactive, vmax,
                                          phot_result)

       -- Evaluate fmid at the point lambda=xmid
       -- fmid will be an increasing function of xmid, with a solution
       -- (fmid=0) between x1 and x2

       -- Second term is total daytime photosynthesis (mm/m2/day) implied by
       -- canopy conductance and current guess for lambda (xmid)
       -- Eqn 18, Haxeltine & Prentice 1996

        let fmid = phot_result.adtmm / fpc - gcphot * (1 - xmid)
        let rtbis = if (fmid < 0) then xmid else rtbis
        in (b,dx,xmid,phot_result,fmid)

      -- bvoc
      let lambda = xmid
      in (pft, phot_result, lambda)
