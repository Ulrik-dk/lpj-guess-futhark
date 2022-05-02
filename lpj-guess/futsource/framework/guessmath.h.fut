--#include <assert.h>
--#include "archive.h"

--#include <cmath>

let PI : f64 = 4 * f64.atan 1.0

let DEGTORAD : f64 = PI / 180.0

let SECS_PER_DAY : f64 = 24.0*60.0*60.0
let M3_PER_MM3: f64 = 1E-9
let MM3_PER_M3: f64 = 1E9
let HA_PER_M2: f64 = 1E-4
let M2_PER_HA: f64 = 1E4
let CM2_PER_M2: f64 = 1E4
let MM2_PER_M2: f64 = 1E6
let SQ_M : f64 = 1.0 -- 1 m2
let M_PER_MM: f64 = 1E-3
let CM_PER_M: f64 = 1E2
let CM_PER_MM: f64 = 1E-1
let MM_PER_M: f64 = 1E3
let MM_PER_CM: f64 = 1E1
let KG_PER_G: f64 = 1E-3
let G_PER_KG: f64 = 1E3
let G_PER_MG: f64 = 1E-3
let MG_PER_G: f64 = 1E3 -- Milligrams
let KG_PER_MT: f64 = 1E3 -- Metric tonnes
let MMOL_PER_MOL : f64 = 1E3
let J_PER_KJ: f64 = 1E3
let KMH_PER_MS : f64 = 3.6
let FRACT_TO_PERCENT : f64 = 100.0
let PERCENT_TO_FRACT : f64 = 0.01
let R_EARTH : f64 = 6371.2213 -- mean earth-radius[km]
