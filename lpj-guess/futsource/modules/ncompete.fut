open import "../framework/guessmath"
open import "../framework/guess"

type NCompetingIndividual = {
  ndemand : real,
  strength : real,
  fnuptake : real
}

let NCompetingIndividual(individual : Individual) : NCompetingIndividual = {
  ndemand = individual.ndemand,
  strength = 0.0,
  fnuptake = individual.fnuptake
} -- FIXME: a guess at a constructor

let ncompete [n] (individuals : [n]NCompetingIndividual, nmass_avail : real) =
  let nsupply : real = nmass_avail    -- Nitrogen available for uptake
  -- calculate total nitrogen uptake strength, and set all uptake to
  -- zero initially
  let total_ups = reduce (+) realzero <| map (\i -> i.strength) individuals -- Total uptake strength

  let individuals = map (\i -> {ndemand=i.ndemand, strength=i.strength, fnuptake=0}) individuals

  let full_uptake : bool = true
                    -- If an individual could get more than its demand, then
                    -- that indiv gets fnuptake = 1 and everything has to be
                    -- redone for all indiv with fnuptake < 1 as more nitrogen
                    -- could be taken up per unit strength (starts with true to
                    -- get into while loop)

  let (_, individuals, _, _) =
  loop (full_uptake, individuals, nsupply, total_ups)
  while (full_uptake) do
    let full_uptake = false

    let ratio_uptake : real =        -- Nitrogen per uptake strength
    -- decide how much nitrogen that will be taken up by each uptake strength
      if (total_ups > 0.0) then nsupply / total_ups
      else 0.0
    -- Go through individuals
    let (full_uptake, individuals, nsupply, total_ups, ratio_uptake, i) =
    loop (full_uptake, individuals, nsupply, total_ups, ratio_uptake, i) =
         (full_uptake, individuals, nsupply, total_ups, ratio_uptake, 0)
    while (i < n && !full_uptake) do
      let indiv : NCompetingIndividual = individuals[i]
      -- if fnuptake doesn't meet indiv nitrogen demand, then calculate a new value for fnuptake

      let (indiv, nsupply, total_ups, full_uptake) =
        if (indiv.fnuptake != 1.0) then
          let entitlement : real = ratio_uptake * indiv.strength in
          if (entitlement > indiv.ndemand && !negligible(indiv.ndemand)) then
            let indiv = indiv with fnuptake = 1.0
            let nsupply = nsupply - indiv.ndemand
            let total_ups = total_ups - indiv.strength
            let full_uptake = true
            in (indiv, nsupply, total_ups, full_uptake)
          else if (indiv.ndemand > 0.0) then
            let indiv = indiv with fnuptake = min(1.0, entitlement / indiv.ndemand)
            in (indiv, nsupply, total_ups, full_uptake)
          else
            let indiv = indiv with fnuptake = 0.0
            in (indiv, nsupply, total_ups, full_uptake)
        else (indiv, nsupply, total_ups, full_uptake)

      let individuals = individuals with [i] = indiv
      in (full_uptake, individuals, nsupply, total_ups, ratio_uptake, i)
    in (full_uptake, individuals, nsupply, total_ups)
  in individuals
