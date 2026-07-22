"""
Run the two models on one forcing file, configured to agree.

Every place the two readers could drift apart is set EXPLICITLY on both sides
rather than left to a default, because most of the defaults disagree:

  * `load_forcing_file` defaults `air_temperature_in_celsius=true` and
    `precipitation_in_mmwe_day=true`; the chion namelist defaults both to
    .FALSE. The forcing file is written in K and kg m-2 s-1, so both are set
    false here.
  * `load_forcing_file` reads LWD / SHF / LHF / RHZ by default if present.
    `chion_grid.x` reads none of them. They are passed as `nothing` so that a
    stray variable in a future forcing file cannot make the two runs differ.
  * `wind_default` must match chion's ctrl:wind_default.

DIURNAL SUBSTEPPING IS OFF, deliberately. Enabling it changes the albedo
scheme rather than only the shortwave resolution -- the aging law carries no dt
(trap 5), so it fires once per substep (upstream defect 19), and snowfall
brightening is non-linear under substepping (defect 21). It is also the only
consumer of day_of_year / solar_longitude_deg, which the two drivers derive
differently (chion_grid.f90 from modulo(time,365)+1, Chion.jl from a calendar
axis). Comparing it would be comparing two known-divergent schemes.
"""

using NCDatasets
using Chion

const CHION_ROOT = normpath(joinpath(@__DIR__, ".."))

"""
Machine epsilon of a chion build, used to express differences in ulp.

Note this is NOT the resolution of the comparison -- see EPS_FILE in
validate.jl. Chion.jl computes in Float64 but WRITES Float32, so no comparison
through its NetCDF output can resolve below Float32 storage quantization
whatever precision either model runs at.
"""
eps_of(precision::Symbol) = Float64(precision === :dp ? eps(Float64) : eps(Float32))

bindir(precision::Symbol) = precision === :dp ? "libchion/bin-dp" : "libchion/bin"

"""
    run_chion(; precision, forcing, outfile, workdir, model, dt_out, nml_extra)

Write a namelist and run `chion_grid.x` in `workdir`. Returns the output path.

`chion_grid.x` writes to its current working directory and requires
`input/chion_defaults.nml` to be reachable from there, so the run directory gets
a symlink to the repository's `input/`.
"""
function run_chion(; precision::Symbol, forcing::AbstractString,
                   outfile::AbstractString, workdir::AbstractString,
                   model::AbstractString="bessi", dt_out::Float64=1.0,
                   dt::Float64=-1.0, nml_extra::AbstractString="")
    mkpath(workdir)
    link = joinpath(workdir, "input")
    islink(link) || ispath(link) || symlink(joinpath(CHION_ROOT, "input"), link)

    nml = joinpath(workdir, "chion_$(model)_$(precision).nml")
    open(nml, "w") do io
        print(io, """
&ctrl
    file_forcing        = "$(abspath(forcing))"
    file_out            = "$(outfile)"
    name_x              = "x"
    name_y              = "y"
    name_time           = "time"
    name_t2m            = "TT"
    name_sf             = "SF"
    name_rf             = "RF"
    name_swd            = "SWD"
    name_mask           = "mask"
    mask_threshold      = 0.0
    name_lat            = "LAT"
    name_zs             = "SH"
    t2m_in_celsius      = .FALSE.
    precip_in_mmwe_day  = .FALSE.
    wind_default        = 5.0
    dt                  = $(dt)
    dt_out              = $(dt_out)
/

&chion
    model               = "$(model)"
    phys_const          = "Earth"
    phys_const_file     = "input/chion_phys_const.nml"
    restart             = "none"
/

&bessi
    Ntot                = 15
    mass_max            = 500.0
    mass_split          = 300.0
    mass_min            = 100.0
    density_init        = 300.0
    temperature_init    = 273.0
    diurnal_shortwave_substeps            = .FALSE.
    diurnal_shortwave_threshold           = 0.0
    diurnal_shortwave_max_substeps        = 3
    diurnal_shortwave_min_air_temperature = 265.15
    diurnal_temperature_cycle             = .FALSE.
    diurnal_temperature_amplitude         = 5.0
/

&pdd
    pdd_method          = "pism"
    ddf_snow            = 3.0
    ddf_ice             = 8.0
    refreezing_fraction = 0.6
    temperature_sigma   = 5.0
/
$(nml_extra)
""")
    end

    exe = joinpath(CHION_ROOT, bindir(precision), "chion_grid.x")
    isfile(exe) || error("$exe not built. Run: make grid" *
                         (precision === :dp ? " precision=dp" : ""))
    logfile = joinpath(workdir, "chion_$(model)_$(precision).log")
    open(logfile, "w") do log
        run(pipeline(Cmd(`$exe $(basename(nml))`; dir=workdir); stdout=log, stderr=log))
    end
    return joinpath(workdir, outfile)
end

"""
    run_julia_bessi(; forcing, outfile, workdir, ntot)

Run Chion.jl's BESSI on the same file. `netcdf_variables` is the explicit 18-var
list: requesting `latent_heat_flux` would flip the run into monthly-aggregation
mode (`_uses_monthly_output`), which writes one record per month instead of one
per step and would not be comparable.
"""
function run_julia_bessi(; forcing::AbstractString, outfile::AbstractString,
                         workdir::AbstractString, ntot::Int=15, years::Int=1)
    mkpath(workdir)
    out = joinpath(workdir, outfile)
    isfile(out) && rm(out)

    loaded = load_forcing_file(
        abspath(forcing);
        x_name="x", y_name="y", time_name="time",
        air_temperature_name="TT", snowfall_name="SF",
        rainfall_name="RF", shortwave_name="SWD",
        wind_speed_name=nothing,
        q_lw_down_name=nothing, q_sh_name=nothing, q_lh_name=nothing,
        relative_humidity_name=nothing, air_pressure_name=nothing,
        prescribed_albedo_name=nothing,
        surface_height_name="SH", latitude_name="LAT",
        mask_name="mask", mask_threshold=0.0,
        air_temperature_in_celsius=false,
        precipitation_in_mmwe_day=false,
        wind_default=5.0,
    )

    model = BESSIModel(loaded.grid; Ntot=ntot, albedo=:dynamic,
                       densification=:bessi, fresh_snow_density=:constant,
                       mass_max=500.0, mass_split=300.0, mass_min=100.0,
                       density_init=300.0, temperature_init=273.0,
                       diurnal_shortwave_substeps=false,
                       diurnal_temperature_cycle=false)

    sim = Simulation(model; forcing=loaded.forcing, years=years,
                     backend=:threads, write_netcdf=true,
                     netcdf_variables=BESSI_VARS, netcdf_path=out,
                     name="wp16_bessi")
    run!(sim)
    return out
end

"""Run Chion.jl's PDD. The explicit var list keeps it in per-step output mode."""
function run_julia_pdd(; forcing::AbstractString, outfile::AbstractString,
                       workdir::AbstractString, years::Int=1)
    mkpath(workdir)
    out = joinpath(workdir, outfile)
    isfile(out) && rm(out)

    loaded = load_forcing_file(
        abspath(forcing);
        x_name="x", y_name="y", time_name="time",
        air_temperature_name="TT", snowfall_name="SF",
        rainfall_name="RF", shortwave_name="SWD",
        wind_speed_name=nothing,
        q_lw_down_name=nothing, q_sh_name=nothing, q_lh_name=nothing,
        relative_humidity_name=nothing, air_pressure_name=nothing,
        prescribed_albedo_name=nothing,
        surface_height_name="SH", latitude_name="LAT",
        mask_name="mask", mask_threshold=0.0,
        air_temperature_in_celsius=false,
        precipitation_in_mmwe_day=false,
        wind_default=5.0,
    )

    model = PDDModel(loaded.grid; ddf_snow=3.0, ddf_ice=8.0,
                     refreezing_fraction=0.6, temperature_sigma=5.0)

    sim = Simulation(model; forcing=loaded.forcing, years=years,
                     backend=:threads, write_netcdf=true,
                     netcdf_variables=PDD_VARS, netcdf_path=out,
                     name="wp16_pdd")
    run!(sim)
    return out
end
