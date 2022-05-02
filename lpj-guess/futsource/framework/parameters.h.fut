--------------------------------------------------------------------------------------/
-- Enums needed by some of the global instruction file parameters defined below

--/ Vegetation 'mode', i.e. what each Individual object represents
-- Can be one of:
--  1. The average characteristics of all individuals comprising a PFT
--     population over the modelled area (standard LPJ mode)
--  2. A cohort of individuals of a PFT that are roughly the same age
--  3. An individual plant
--/
type vegmodetype = i8
let NOVEGMODE : vegmodetype = 0
let INDIVIDUAL : vegmodetype = 1
let COHORT : vegmodetype = 2
let POPULATION : vegmodetype = 3

--/ Land cover type of a stand. NLANDCOVERTYPES keeps count of number of items.
--  NB. set_lc_change_array() must be modified when adding new land cover types
--/
type landcovertype = i8
let URBAN : landcovertype = 0
let CROPLAND : landcovertype = 1
let PASTURE : landcovertype = 2
let FOREST : landcovertype = 3
let NATURAL : landcovertype = 4
let PEATLAND : landcovertype = 5
let BARREN : landcovertype = 6
let NLANDCOVERTYPES : landcovertype = 7

--/ Water uptake parameterisations
-- \see water_uptake in canexch.cpp
type wateruptaketype = i8
let WR_WCONT : wateruptaketype = 0
let WR_ROOTDIST : wateruptaketype = 1
let WR_SMART : wateruptaketype = 2
let WR_SPECIESSPECIFIC : wateruptaketype = 3

--/bvoc: define monoterpene species used
type monoterpenecompoundtype = i64
let APIN : monoterpenecompoundtype = 0
let BPIN : monoterpenecompoundtype = 1
let LIMO : monoterpenecompoundtype = 2
let MYRC : monoterpenecompoundtype = 3
let SABI : monoterpenecompoundtype = 4
let CAMP : monoterpenecompoundtype = 5
let TRIC : monoterpenecompoundtype = 6
let TBOC : monoterpenecompoundtype = 7
let OTHR : monoterpenecompoundtype = 8
let NMTCOMPOUNDTYPES : monoterpenecompoundtype = 9

--/ Fire model setting. Either use
--	One of
--	BLAZE 		Use the BLAZE model to generate fire fluxes
--                      (must be accompanied by ignitionmode; DEFAULT)
--	GLOBFIRM	fire parameterization following Thonicke et al. 2001
--	NOFIRE		no fire model
--/
type firemodeltype = i8
let BLAZE : firemodeltype = 0
let GLOBFIRM : firemodeltype = 1
let NOFIRE : firemodeltype = 2

--/ Type of weathergenerator used
--     One of:
--      GWGEN           Global Weather GENerator (needed by BLAZE, due to
--                      additional rel. humidity and wind; DEFAULT)
--      INTERP          use standard interpolation scheme
--      NONE            Should be set if daily input is used (e.g. in cfinput)
--/
type weathergeneratortype = i8
let GWGEN : weathergeneratortype = 0
let INTERP : weathergeneratortype = 1
let NONE : weathergeneratortype = 2

--/How to determine root distribution in soil layers
type rootdisttype = i8
let ROOTDIST_FIXED : rootdisttype = 0
let ROOTDIST_JACKSON : rootdisttype = 1
