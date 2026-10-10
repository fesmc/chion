# itm

Output variable table for `model = "itm"`, read by `chion_io.f90`.

Names follow `Chion.jl/src/io.jl` `ITM_OUTPUT_VARS` one-to-one, so a chion and a
Chion.jl ITM file are directly comparable. Units and long names are Chion.jl's
(`NETCDF_METADATA`) except for `melt`, `runoff` and `refreezing`: Chion.jl
looks these up under BESSI's cumulative entries ("Cumulative melt", mmWE),
but for ITM they are the step's rates, so chion writes the rate unit and long
name (reported upstream).

The rates [mmWE day-1] are those of the step just completed; the `*_cum`
fields and `smb_ice` are cumulative [mmWE]. mm w.e. is kg m-2 by definition.

**`smb` differs from BESSI's and PDD's `smb`.** Here it is ITM's TOTAL surface
mass balance rate (`sf + rf - runoff`) in mmWE day-1, as in Chion.jl. BESSI and
PDD write chion's ice-facing flux from `chion_get_smb` in kg m-2 s-1 under that
name. ITM's ice-facing rate is `smbi` (mmWE day-1); `chion_get_smb` returns
`smbi/86400` in kg m-2 s-1.

The `dimensions` column is LOGICAL: `column` becomes `xc, yc` when a spatial
mapping is attached with `chion_set_grid`, and a bare `column` dimension
otherwise. `time` is appended by the writer and is always unlimited.

| id | variable             | dimensions    | units       | long_name                                            |
|----|----------------------|---------------|-------------|------------------------------------------------------|
|  1 | H_snow               | column        | mmWE        | ITM snowpack water equivalent                        |
|  2 | alb_s                | column        | 1           | ITM surface albedo                                   |
|  3 | smb                  | column        | mmWE day-1  | ITM total surface mass balance rate                  |
|  4 | smbi                 | column        | mmWE day-1  | ITM ice-facing mass balance rate                     |
|  5 | melt                 | column        | mmWE day-1  | ITM melt rate                                        |
|  6 | runoff               | column        | mmWE day-1  | ITM runoff rate                                      |
|  7 | refreezing           | column        | mmWE day-1  | ITM refreezing rate                                  |
|  8 | Tsrf                 | column        | K           | Surface temperature                                  |
|  9 | melt_net             | column        | mmWE day-1  | ITM net melt rate                                    |
| 10 | smb_cum              | column        | mmWE        | ITM cumulative surface mass balance                  |
| 11 | smb_ice              | column        | mmWE        | Net mass forcing to the ice sheet                    |
| 12 | melt_cum             | column        | mmWE        | ITM cumulative melt                                  |
| 13 | runoff_cum           | column        | mmWE        | ITM cumulative runoff                                |
| 14 | refreezing_cum       | column        | mmWE        | ITM cumulative refreezing                            |

Notes.

* Fortran state names (`chn%itm%now`, and the restart file) keep smbpal's:
  `refrz` = `refreezing`, `tsrf` = `Tsrf`, `smbi_cum` = `smb_ice`,
  `refrz_cum` = `refreezing_cum`.
* `melt_net` is `refrz - melt` (on ice) or `refrz - snow melt` (land), the
  input of the firn-warming surface temperature.
* `smb_ice` (`smbi_cum`, `snow_to_ice + refrz - melted_ice` integrated) is the
  ice-facing balance; `smb_cum` the whole-column one (`sf + rf - runoff`).
