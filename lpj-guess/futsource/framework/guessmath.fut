----------section from guessmath.h----------
--#include <assert.h>
--#include "archive.h"
--#include <cmath>
open import "../futhark-extras"

let PI : real = 4 * atan 1.0

let DEGTORAD : real = PI / 180.0

let SECS_PER_DAY : real = 24.0*60.0*60.0
let M3_PER_MM3: real = 1E-9
let MM3_PER_M3: real = 1E9
let HA_PER_M2: real = 1E-4
let M2_PER_HA: real = 1E4
let CM2_PER_M2: real = 1E4
let MM2_PER_M2: real = 1E6
let SQ_M : real = 1.0 -- 1 m2
let M_PER_MM: real = 1E-3
let CM_PER_M: real = 1E2
let CM_PER_MM: real = 1E-1
let MM_PER_M: real = 1E3
let MM_PER_CM: real = 1E1
let KG_PER_G: real = 1E-3
let G_PER_KG: real = 1E3
let G_PER_MG: real = 1E-3
let MG_PER_G: real = 1E3 -- Milligrams
let KG_PER_MT: real = 1E3 -- Metric tonnes
let MMOL_PER_MOL : real = 1E3
let J_PER_KJ: real = 1E3
let KMH_PER_MS : real = 3.6
let FRACT_TO_PERCENT : real = 100.0
let PERCENT_TO_FRACT : real = 0.01
let R_EARTH : real = 6371.2213 -- mean earth-radius[km]


-- TODO default limit should be 0.0
let negligible (dval : real) : bool = abs(dval) < 1.0e-30
-- Returns true if |dval| < EPSILON, otherwise false

  -- TODO default limit should be 0.0
let negligible_2 (dval : real, limit : real) : bool =
	-- Returns true if |dval| < EPSILON, otherwise false
	if (limit == 0.0) then (abs dval) < (pow (10.0, limit)) else abs(dval) < 1.0e-30


-- TODO: more to be translated here
