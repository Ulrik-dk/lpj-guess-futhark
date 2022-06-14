type real = f64
type int = i64
type uint = u64

let ifbvoc : int = 0

let realzero : real = 0.0
let Date_subdaily : int = 1 -- get this somewhere else
let Date_MAX_YEAR_LENGTH : int = 365
let Date_MAX_YEAR_LENGTH_plusone : int = 366
type xtring = int -- we wont care about strings
let euler : real = 2.71828
let pow(a: real, b: real) = a ** b
let exp(a: real) = pow(euler, a)
let log10 = f64.log10
let log = log10
let sqrt = f64.sqrt
let abs = f64.abs
let sin = f64.sin
let cos = f64.cos

let nan = f64.nan
let intnan = i64.lowest --this may serve the purpose
let atan = f64.atan
let min (a: f64, b: f64) = if a > b then b else a
let max (a: f64, b: f64) = if a > b then a else b

let boolFromReal = bool.f64
let realFromInt = f64.i64
let intFromReal = i64.f64

type enum_type = i64


let npft : int = 6 -- TODO: this is a dynamic global variable that should be aquired otherwise than this!



-- from global.ins
let run_landcover : bool = false
let lambda_max : real = 0.8
--parameters.cpp?
let disturb_pasture : bool = false
