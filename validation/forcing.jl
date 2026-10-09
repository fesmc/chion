"""
Scenario definitions and forcing-file generation for the WP16 harness.

ONE forcing file drives both models. That is the whole point: any difference in
the output is then attributable to the models, not to two generators that drift
apart. chion reads it with `chion_grid.x` and Chion.jl with `load_forcing_file`,
and the two readers are configured to agree field for field (see runners.jl).

WHAT IS DELIBERATELY ABSENT FROM THE FILE
-----------------------------------------
Only TT, SF, RF, SWD, LAT, SH and mask are written (plus ITM's HI and PDDA in
the ITM file, and RHZ and PS in the humidity-on BESSI file, see below).
`load_forcing_file` picks up
LWD / SHF / LHF / RHZ by their default names whenever they are present, while
`chion_grid.x` reads none of them, so writing any of those would silently hand
Chion.jl forcing that chion never sees. Omitting RHZ additionally keeps humidity
forcing off in the default configurations (before Chion.jl 03bb445 this was
required for the BESSI mass-closure identity: upstream defect 1, the unclipped
vapour diagnostic, left ~105 kg m-2 unclosed against 192 kg m-2 reported).

The humidity-on BESSI configuration (WP11) gets a file of its own carrying RHZ
and PS, both uniform: chion_grid.x has no humidity or pressure reader, only
`rh_default` and the sea-level pressure it fills in (DEF_SEA_LEVEL_AIR_PRESSURE),
so the file carries exactly those values for Chion.jl to read -- PS also
overrides Chion.jl's barometric pressure from SH, which chion does not apply.

The snow/rain split is computed HERE and written as two separate variables, so
neither model applies a partition rule of its own.
"""

using NCDatasets
using Dates

const YEAR_LENGTH = 365.0

"""chion_grid.x's air pressure (chion_defs DEF_SEA_LEVEL_AIR_PRESSURE) [Pa]."""
const SEA_LEVEL_PRESSURE = 101325.0

"""
One independent column of the BESSI comparison. Each scenario is a column of a
single (ncol x 1) grid, so all of them are exercised by one run of each model --
which also exercises the multi-column path rather than only the ncol=1 path.
"""
struct Scenario
    name::String
    what::String          # which code paths this column is here to exercise
    t2m_mean::Float64     # [K]
    t2m_amp::Float64      # [K] half-amplitude of the annual cycle
    pr::Float64           # [kg m-2 s-1] precipitation rate
    sw_mean::Float64      # [W m-2]
    sw_amp::Float64       # [W m-2]
end

"""
The four BESSI coverage requirements from docs/PLAN.md WP16.

The settings come from the D18 retuning of par/chion_column.nml, which recorded
which regimes actually fire: colder or drier builds layers but never melts,
warmer ablates to bare ice and never splits a layer.
"""
const BESSI_SCENARIOS = [
    Scenario("cold_dry",
             "layer split/merge and densification, no melt",
             245.0, 8.0, 6.0e-5, 60.0, 55.0),
    Scenario("melting",
             "energy solve, melt, percolation, refreezing (percolation-zone firn)",
             258.15, 16.0, 1.2e-4, 165.0, 145.0),
    Scenario("bare_recover",
             "ablation to bare ice, the early-return bare path, then recovery",
             268.0, 18.0, 5.0e-5, 200.0, 160.0),
    Scenario("ntot_capacity",
             "Ntot capacity: bottom merge and the snow-depth cap",
             250.0, 10.0, 6.0e-4, 60.0, 55.0),
]

"""
Static per-column inputs for ITM: latitude [deg N], surface height [m], ice
thickness [m] and annual positive degree days [K d]. They select ITM's
background albedo (ocean where `zs <= 0`, land where `h_ice == 0`, ice
otherwise), its critical snow depth (desert/tundra/forest by `pdd`) and its
latitude offset. BESSI and PDD forcing omits them (uniform LAT/SH only).
"""
struct Geometry
    lat::Float64
    zs::Float64
    h_ice::Float64
    pdd::Float64
end

"""
The ITM columns: each reaches one of ITM's background-albedo branches and a
distinct part of the budget. All start at H_snow = H_snow_max (both models'
cold start), so the cap is active from the first step wherever snow
accumulates. Coverage is asserted in validate.jl.
"""
const ITM_SCENARIOS = [
    (Scenario("ice_ablation",
              "ice background; snow melts out, then bare-ice melt; rain refreezing",
              275.0, 12.0, 1.5e-5, 200.0, 180.0),
     Geometry(67.0, 500.0, 500.0, 800.0)),
    (Scenario("ice_accum_cap",
              "ice background; cold and snowy, H_snow_max cap exports to ice",
              245.0, 10.0, 2.0e-4, 120.0, 110.0),
     Geometry(72.0, 3000.0, 3000.0, 0.0)),
    (Scenario("land_seasonal",
              "land background (tundra PDDs), seasonal snow",
              273.0, 15.0, 3.0e-5, 180.0, 160.0),
     Geometry(62.0, 300.0, 0.0, 500.0)),
    (Scenario("ocean",
              "ocean background (zs <= 0), forest critical depth",
              278.0, 8.0, 4.0e-5, 180.0, 160.0),
     Geometry(60.0, 0.0, 0.0, 1500.0)),
]

seasonal(day, phase=200.0) = cos(2pi * (day - phase) / YEAR_LENGTH)

air_temperature(s::Scenario, day) = s.t2m_mean + s.t2m_amp * seasonal(day)
shortwave(s::Scenario, day) = max(s.sw_mean + s.sw_amp * seasonal(day), 0.0)

"""
    write_forcing(path, scenarios; nstep, dt_days, t_snow_max, geometry,
                  relative_humidity)

Write the shared forcing file. With `relative_humidity` [1] given, uniform RHZ
and PS (`SEA_LEVEL_PRESSURE`) are added for the humidity-on configuration.

The time axis is written in two forms, because the two readers accept different
ones and neither accepts both: a CF numeric axis for `chion_grid.x`, and
YYYY/MM/DD/HH integer variables for Chion.jl. Both encode the same instants, so
`infer_dt_days` and `times(2)-times(1)` agree on `dt_days`. See the long note at
the write site for why the CF axis alone is not enough.
"""
function write_forcing(path::AbstractString, scenarios::Vector{Scenario};
                       nstep::Int, dt_days::Float64=1.0, t_snow_max::Float64=273.15,
                       geometry::Union{Nothing,Vector{Geometry}}=nothing,
                       relative_humidity::Union{Nothing,Float64}=nothing)
    ncol = length(scenarios)
    nx, ny = ncol, 1

    times = collect(0:(nstep - 1)) .* dt_days
    TT = zeros(nstep, ny, nx)
    SF = zeros(nstep, ny, nx)
    RF = zeros(nstep, ny, nx)
    SWD = zeros(nstep, ny, nx)

    for (i, s) in enumerate(scenarios), k in 1:nstep
        # The forcing sample at index k applies over the step [t_k, t_k + dt),
        # evaluated at the interval midpoint. Named `tday`, not `day`, so it
        # cannot shadow Dates.day, which is used below.
        tday = times[k] + 0.5 * dt_days
        t = air_temperature(s, tday)
        TT[k, 1, i] = t
        SWD[k, 1, i] = shortwave(s, tday)
        # Snow/rain split fixed here so neither model partitions on its own.
        if t < t_snow_max
            SF[k, 1, i] = s.pr
        else
            RF[k, 1, i] = s.pr
        end
    end

    isfile(path) && rm(path)
    NCDataset(path, "c") do ds
        defDim(ds, "x", nx); defDim(ds, "y", ny); defDim(ds, "time", nstep)

        v = defVar(ds, "x", Float64, ("x",)); v[:] = collect(0:(nx - 1)) .* 10.0
        v.attrib["units"] = "km"
        v = defVar(ds, "y", Float64, ("y",)); v[:] = collect(0:(ny - 1)) .* 10.0
        v.attrib["units"] = "km"
        # THE TIME AXIS IS WRITTEN TWICE, deliberately. See the note below.
        v = defVar(ds, "time", Float64, ("time",)); v[:] = times
        v.attrib["units"] = "days since 2000-01-01 00:00:00"
        v.attrib["calendar"] = "proleptic_gregorian"

        # Chion.jl cannot read the axis above. `_read_time_values`
        # (dataloaders.jl:17) reads `ds[time_name].var[:]` -- the UNDECODED
        # accessor -- and then tests `eltype(raw) <: DateTime`. NCDatasets only
        # returns DateTime from the decoded accessor `ds[time_name]`, so that
        # test can never be true for a CF-compliant numeric time axis, however
        # correctly it is written. The function then falls through to a
        # synthetic axis of ONE DAY PER RECORD (dataloaders.jl:38), silently
        # discarding the file's real timing with no warning.
        #
        # That is not cosmetic. chion_grid.x reads the raw numbers and infers
        # dt correctly, so at dt = 30 d the two models stepped 30x apart AND
        # took different PDD branches -- Chion.jl's PISM/Calov-Greve integral is
        # auto-selected on 27 <= dt_days <= 32, so at a phantom dt = 1 it used
        # the simple form instead. Every PDD field was wrong by a factor ~30.
        #
        # The YYYY/MM/DD/HH branch (dataloaders.jl:30) is the only one that
        # actually fires, so the harness supplies it. Reported upstream; see
        # docs/porting_notes.md.
        epoch = DateTime(2000, 1, 1, 0)
        stamps = [epoch + Millisecond(round(Int, t * 86_400_000)) for t in times]
        for (nm, f) in (("YYYY", year), ("MM", month), ("DD", day), ("HH", hour))
            vv = defVar(ds, nm, Int32, ("time",))
            vv[:] = Int32[f(s) for s in stamps]
            vv.attrib["long_name"] = "Calendar $nm, for Chion.jl's _read_time_values"
        end

        humid = relative_humidity === nothing ? () :
            (("RHZ", fill(relative_humidity, size(TT)), "1", "Relative humidity"),
             ("PS", fill(SEA_LEVEL_PRESSURE, size(TT)), "Pa", "Air pressure"))
        for (name, data, units, long) in (
                ("TT", TT, "K", "Air temperature"),
                ("SF", SF, "kg m-2 s-1", "Snowfall rate"),
                ("RF", RF, "kg m-2 s-1", "Rainfall rate"),
                ("SWD", SWD, "W m-2", "Downward shortwave radiation"),
                humid...)
            # Declared ("x","y","time") in Julia order so the file carries
            # (time,y,x) in C order -- the layout chion_grid.x reads.
            vv = defVar(ds, name, Float64, ("x", "y", "time"))
            vv[:, :, :] = permutedims(data, (3, 2, 1))
            vv.attrib["units"] = units
            vv.attrib["long_name"] = long
        end

        # Static fields. Uniform unless per-column geometry is given (ITM),
        # which also adds ITM's HI and PDDA; Chion.jl's load_forcing_file
        # reads neither, so they cannot leak into the BESSI or PDD runs.
        statics = geometry === nothing ?
            (("mask", fill(1.0, ncol), "1"),
             ("LAT", fill(70.0, ncol), "degrees_north"),
             ("SH", fill(1500.0, ncol), "m")) :
            (("mask", fill(1.0, ncol), "1"),
             ("LAT", [g.lat for g in geometry], "degrees_north"),
             ("SH", [g.zs for g in geometry], "m"),
             ("HI", [g.h_ice for g in geometry], "m"),
             ("PDDA", [g.pdd for g in geometry], "K d"))
        for (name, val, units) in statics
            vv = defVar(ds, name, Float64, ("x", "y"))
            vv[:, :] = reshape(val, nx, ny)
            vv.attrib["units"] = units
        end

        ds.attrib["title"] = "chion WP16 validation forcing"
        ds.attrib["scenarios"] = join((s.name for s in scenarios), ", ")
    end
    return (path=path, times=times, ncol=ncol, nstep=nstep, dt_days=dt_days)
end
