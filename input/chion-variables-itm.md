# itm

Output variable table for `model = "itm"`, read by `chion_io.f90`.

ITM has NO Chion.jl counterpart (`build_model(:itm,...)` errors upstream), so
there is no `NETCDF_METADATA` entry to match. Where a quantity has the same
meaning as a BESSI/Chion.jl one, the Chion.jl NAME and UNITS are reused
deliberately, so that a host can read `smb_ice`, `runoff`, `melt`,
`refreezing`, `albedo` and `Tsrf` from any chion output file without knowing
which model produced it. The remaining entries carry smbpal's own names.

Note ITM works internally in [mm w.e.] and [mm w.e. d-1]; mm w.e. is kg m-2 by
definition, so `mmWE` is used throughout for consistency with the other tables.

The `dimensions` column is LOGICAL: `column` becomes `xc, yc` when a spatial
mapping is attached with `chion_set_grid`, and a bare `column` dimension
otherwise. `time` is appended by the writer and is always unlimited.

| id | variable             | dimensions    | units       | long_name                                            |
|----|----------------------|---------------|-------------|------------------------------------------------------|
|  1 | H_snow               | column        | mmWE        | Snowpack thickness                                   |
|  2 | albedo               | column        | 1           | Surface albedo                                       |
|  3 | Tsrf                 | column        | K           | Surface temperature                                  |
|  4 | smb_ice              | column        | mmWE        | Net mass forcing to the ice sheet                    |
|  5 | runoff               | column        | mmWE        | Cumulative runoff                                    |
|  6 | melt                 | column        | mmWE        | Cumulative melt                                      |
|  7 | refreezing           | column        | mmWE        | Cumulative refreezing                                |
|  8 | smb_total            | column        | mmWE        | Cumulative whole-column surface mass balance         |
|  9 | smb                  | column        | kg m-2 s-1  | Net mass flux to the ice sheet                       |

Notes.

* `albedo` is smbpal's `alb_s`, `Tsrf` is its `tsrf`.
* `smb_ice` is `smbi_cum`, the cumulative form of smbpal's `smbi`
  (`snow_to_ice + refrz - melted_ice`). `smb_total` is `smb_cum`, the cumulative
  form of smbpal's whole-column `smb` (`sf + rf - runoff`). The two are
  different quantities and both are reported, because smbpal reports both.
* `runoff`, `melt` and `refreezing` are the cumulative accumulators
  (`runoff_cum`, `melt_cum`, `refrz_cum`), not smbpal's per-step rates, so they
  have the same meaning as BESSI's fields of the same name.
* `smb` (id 9) is `chion_get_smb`, i.e. `smbi_cum` differenced over the step and
  converted to kg m-2 s-1. It reproduces `itm%now%smbi/86400` exactly.
