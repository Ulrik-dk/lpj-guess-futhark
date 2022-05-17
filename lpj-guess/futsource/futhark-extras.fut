type real = f64
type int = i64

let Date_MAX_YEAR_LENGTH_plusone : int = 367
type xtring = int -- we wont care about strings
let euler : real = 2.71828
let pow(a: real, b: real) = a ** b
let exp(a: real) = pow(euler, a)
let log10 = f64.log10
let sqrt = f64.sqrt
let abs = f64.abs

let nan = f64.nan
let atan = f64.atan
let min (a: f64, b: f64) = if a > b then b else a
let max (a: f64, b: f64) = if a > b then a else b

let boolFromReal = bool.f64
let realFromInt = f64.i64
let intFromReal = i64.f64

type enum_type = i64
