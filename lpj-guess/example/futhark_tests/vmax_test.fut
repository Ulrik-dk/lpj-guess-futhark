open import "../../futsource/everything"

let input =
  let b : f64 = 0.015
  let c1 : f64 = 0.0535267
  let c2 : f64 = 0.314037
  let apar : f64 = 4.74968e+06
  let tscal : f64 = 0.962908
  let daylength : f64 = 10.8269
  let temp : f64 = 21.795
  let nactive : f64 = 1
  let ifnlimvmax : bool = false
  in (b, c1, c2, apar, tscal, daylength, temp, nactive, ifnlimvmax)

-- Autogenerated test of vmax output (vm, vmaxnlim, nactive_opt) field: vm
-- ==
-- entry: vm_test
-- input {}
-- output { 180.282661 }

entry vm_test =
  let (vm, vmaxnlim, nactive_opt) = vmax input
  in vm

-- Autogenerated test of vmax output (vm, vmaxnlim, nactive_opt) field: vmaxnlim
-- ==
-- entry: vmaxnlim_test
-- input {}
-- output { 1.000000 }

entry vmaxnlim_test =
  let (vm, vmaxnlim, nactive_opt) = vmax input
  in vmaxnlim

-- Autogenerated test of vmax output (vm, vmaxnlim, nactive_opt) field: nactive_opt
-- ==
-- entry: nactive_opt_test
-- input {}
-- output { 0.012033 }

entry nactive_opt_test =
  let (vm, vmaxnlim, nactive_opt) = vmax input
  in nactive_opt

