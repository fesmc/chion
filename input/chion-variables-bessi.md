# bessi

Output variable table for `model = "bessi"`, read by `chion_io.f90`
(`load_var_io_table`, fesm-utils `variable_io`). Adding a variable here is not
enough on its own: it also needs a `case` in `chion_write_var_bessi`.

Names, units and long names are taken verbatim from `Chion.jl/src/io.jl`
`NETCDF_METADATA` so that a chion output file and a Chion.jl output file are
directly comparable variable by variable (WP16). Four entries (`vapor_mass`,
`smb`, `albedo_snow`, `ice_temperature`) have no Chion.jl counterpart and are
flagged below. Chion.jl's monthly-only `surface_smb` and `latent_heat_flux` are
not written: chion has no monthly writer, and a host forms its surface SMB
from the cumulative fields itself (precipitation - runoff - sublimation).

The `dimensions` column is LOGICAL, not literal. `column` is expanded by the
writer to `xc, yc` when a spatial mapping has been attached with
`chion_set_grid`, and to a bare `column` dimension otherwise. `time` is
appended by the writer and is always the unlimited dimension.

| id | variable             | dimensions    | units       | long_name                                            |
|----|----------------------|---------------|-------------|------------------------------------------------------|
|  1 | thickness            | column        | m           | Snow thickness                                       |
|  2 | wet_mass             | column        | mmWE        | Snow wet mass                                        |
|  3 | bulk_density         | column        | kg m-3      | Bulk snow density                                    |
|  4 | liquid_water         | column        | kg m-2      | Liquid water mass                                    |
|  5 | mass_base            | column        | mmWE        | Firn mass exported to the ice model                  |
|  6 | smb_ice              | column        | mmWE        | Net mass forcing to the ice sheet                    |
|  7 | runoff               | column        | mmWE        | Cumulative runoff                                    |
|  8 | melt                 | column        | mmWE        | Cumulative melt                                      |
|  9 | refreezing           | column        | mmWE        | Cumulative refreezing                                |
| 10 | sublimation          | column        | mmWE        | Cumulative sublimation                               |
| 11 | vapor_mass           | column        | mmWE        | Cumulative surface vapour mass flux                  |
| 12 | latent_heat_flux_sum | column        | W m-2       | Integrated turbulent latent heat flux                |
| 13 | Tsrf                 | column        | K           | Surface temperature                                  |
| 14 | albedo               | column        | 1           | Surface albedo                                       |
| 15 | N                    | column        | 1           | Number of active snow layers                         |
| 16 | smb                  | column        | kg m-2 s-1  | Net mass flux to the ice sheet                       |
| 17 | mass                 | layer, column | kg m-2      | Layer snow mass                                      |
| 18 | mass_w               | layer, column | kg m-2      | Layer liquid-water mass                              |
| 19 | density              | layer, column | kg m-3      | Layer density                                        |
| 20 | temperature          | layer, column | K           | Layer temperature                                    |
| 21 | snow_age_days        | column        | day         | Time since the latest snowfall event                 |
| 22 | albedo_snow          | column        | 1           | Snow albedo before the thin-snow blend               |
| 23 | ice_temperature      | ice_layer, column | K       | Ice substrate layer temperature                      |

Notes.

* `vapor_mass` (id 11) is a chion addition. Chion.jl carries the accumulator in
  `BESSIState` but does not list it in `NETCDF_METADATA`, so only
  `sublimation` is comparable directly. Sign convention: positive = deposition.
* `smb` (id 16) is a chion addition: the model-agnostic ice-facing flux returned
  by `chion_get_smb`, averaged over the step just completed. It is 0 before the
  first `chion_update`. Chion.jl has no equivalent.
* `albedo_snow` (id 22) is a chion addition (D40): the snow's own albedo, which
  the albedo schemes age and refresh. `albedo` (id 14) is what the surface
  energy balance sees, `albedo_snow` blended with the background by the
  snow-cover fraction; the two agree under full snow cover and in Chion.jl,
  which has no blend.
* `ice_temperature` (id 23) is a chion addition: the thermal ice substrate
  (Chion.jl 03bb445 keeps it in the state but does not write it). It is
  written only when `ice_substrate_layers > 0`, on its own `ice_layer`
  dimension (layer 1 at the snow/ice interface, thicknesses doubling down).
* `smb` (id 16) is the ice-facing flux in kg m-2 s-1 here and in PDD; ITM's
  `smb` is Chion.jl's total surface mass balance rate in mmWE day-1
  (`input/chion-variables-itm.md`).
* `N` (id 15) is written as a float in output files so that unmapped grid cells
  can carry the missing value. The restart file writes it as an integer.
* `latent_heat_flux_sum` is in W m-2 days, not W m-2 — the accumulator is
  incremented by `flux*dt_days`. The unit string reproduces Chion.jl's, which
  is wrong upstream (`units="W m-2"`); see the report for WP14.
