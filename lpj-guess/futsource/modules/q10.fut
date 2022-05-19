open import "../framework/guessmath"
open import "../futhark-extras"

---------------------------------------------------------------------------------------
--- \file q10.h
--- \brief Q10 calculations for photosynthesis
---
--- Calculations of Q10 values for photosynthesis, formerly
--- placed in canexch.cpp now moved to separate header file
--- because of their application for BVOC calculations as well.
---
--- \author Guy Schurgers (based on LPJ-GUESS 2.1 / Ben Smith)
--- $Date: 2012-07-11 18:10:03 +0200 (Wed, 11 Jul 2012) $
---
---------------------------------------------------------------------------------------

-- WHAT SHOULD THIS FILE CONTAIN?
-- Module header files need normally contain only declarations of functions
-- defined in the module that are to be accessible to the calling framework or
-- to other modules.

-- Constants required for Q10 lookup tables used by photosynthesis
let Q10_MINTEMP : real = -70.0  -- minimum temperature ever (deg C)
let Q10_MAXTEMP : real = 70.0  -- maximum temperature ever (deg C)
let Q10_PRECISION : real = 0.01  -- rounding precision for temperature
let Q10_NDATA : int = (intFromReal) ((Q10_MAXTEMP-Q10_MINTEMP)/Q10_PRECISION + 1.5)
  -- maximum number of values to store in each lookup table

--- Q10 lookup table class
-- Stores pre-calculated temperature-adjusted values based on Q10 and
-- a 25-degree base value.
type LookupQ10 = {
  --- The temperature-adjusted values
  data : [Q10_NDATA]real
}
--- Creates a lookup table
-- \param q10    The Q10 to be used for the table
-- \param base25 Base value for 25 degrees C
let LookupQ10(q10 : real, base25 : real) = map (\i ->
   base25 * pow(q10, (Q10_MINTEMP + (realFromInt i)*Q10_PRECISION - 25.0) / 10.0))
   (iota Q10_NDATA)

-- "Array element" operator
-- \param temp  Temperature (deg C)
-- \returns     Temperature-adjusted value based on Q10 and 25-degree base value
let readQ10 (temp: real, table: [Q10_NDATA]real) =
  -- Element number corresponding to a particular temperature
  let temp = if temp < Q10_MINTEMP then Q10_MINTEMP else if temp > Q10_MAXTEMP then Q10_MAXTEMP else temp

  let i = (intFromReal) ((temp-Q10_MINTEMP)/Q10_PRECISION+0.5)

  in table[i]
