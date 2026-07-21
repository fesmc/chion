# bessi

Output variable table for `model = "bessi"`, read by `chion_io.f90`
(`load_var_io_table`, fesm-utils `variable_io`). Adding a variable here is not
enough on its own: it also needs a `case` in `chion_write_var_bessi`.

Names, units and long names are taken verbatim from `Chion.jl/src/io.jl`
`NETCDF_METADATA` so that a chion output file and a Chion.jl output file are
directly comparable variable by variable (WP16). Two entries (`vapor_mass`,
`smb`) have no Chion.jl counterpart and are flagged below.

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

Notes.

* `vapor_mass` (id 11) is a chion addition. Chion.jl carries the accumulator in
  `BESSIState` but does not list it in `NETCDF_METADATA`, so only
  `sublimation` is comparable directly. Sign convention: positive = deposition.
* `smb` (id 16) is a chion addition: the model-agnostic ice-facing flux returned
  by `chion_get_smb`, averaged over the step just completed. It is 0 before the
  first `chion_update`. Chion.jl has no equivalent.
* `N` (id 15) is written as a float in output files so that unmapped grid cells
  can carry the missing value. The restart file writes it as an integer.
* `latent_heat_flux_sum` is in W m-2 days, not W m-2 — the accumulator is
  incremented by `flux*dt_days`. The unit string reproduces Chion.jl's, which
  is wrong upstream (`units="W m-2"`); see the report for WP14.
