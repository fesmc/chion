# pdd

Output variable table for `model = "pdd"`, read by `chion_io.f90`.

Names and units follow `Chion.jl/src/io.jl` `PDD_OUTPUT_VARS` and
`NETCDF_METADATA`. `smb` is a chion addition (see below).

The `dimensions` column is LOGICAL: `column` becomes `xc, yc` when a spatial
mapping is attached with `chion_set_grid`, and a bare `column` dimension
otherwise. `time` is appended by the writer and is always unlimited.

| id | variable             | dimensions    | units       | long_name                                            |
|----|----------------------|---------------|-------------|------------------------------------------------------|
|  1 | snowpack_swe         | column        | mmWE        | Snowpack water equivalent                            |
|  2 | smb_ice              | column        | mmWE        | Net mass forcing to the ice sheet                    |
|  3 | runoff               | column        | mmWE        | Cumulative runoff                                    |
|  4 | pdd_sum              | column        | degC day    | Cumulative positive degree days                      |
|  5 | smb                  | column        | kg m-2 s-1  | Net mass flux to the ice sheet                       |

Notes.

* `smb_ice` here is Chion.jl's, i.e. a WHOLE-COLUMN mass change rather than the
  ice-facing flux its long name claims (docs/porting_notes.md C2). It is written
  unmodified so that the file is comparable with Chion.jl's.
* `smb` (id 5) is the reconciled ice-facing flux from `chion_get_smb`
  (`smb_ice - snowpack_swe`, differenced over the step). For PDD this is
  identically `-ice_melt` and therefore never positive — see
  docs/porting_notes.md D13.
