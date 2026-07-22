"""
Scenario definitions and forcing-file generation for the WP16 harness.

ONE forcing file drives both models. That is the whole point: any difference in
the output is then attributable to the models, not to two generators that drift
apart. chion reads it with `chion_grid.x` and Chion.jl with `load_forcing_file`,
and the two readers are configured to agree field for field (see runners.jl).

WHAT IS DELIBERATELY ABSENT FROM THE FILE
-----------------------------------------
Only TT, SF, RF, SWD, LAT, SH and mask are written. `load_forcing_file` picks up
LWD / SHF / LHF / RHZ by their default names whenever they are present, while
`chion_grid.x` reads none of them, so writing any of those would silently hand
Chion.jl forcing that chion never sees. Omitting RHZ additionally keeps humidity
forcing off, which is required for the BESSI mass-closure identity to hold at
all -- with humidity on, the vapour diagnostic is unclosed by ~105 kg m-2
against 192 kg m-2 reported (upstream Chion.jl defect 1).

The snow/rain split is computed HERE and written as two separate variables, so
neither model applies a partition rule of its own.
"""

using NCDatasets
using Dates

const YEAR_LENGTH = 365.0

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

seasonal(day, phase=200.0) = cos(2pi * (day - phase) / YEAR_LENGTH)

air_temperature(s::Scenario, day) = s.t2m_mean + s.t2m_amp * seasonal(day)
shortwave(s::Scenario, day) = max(s.sw_mean + s.sw_amp * seasonal(day), 0.0)

"""
    write_forcing(path, scenarios; nstep, dt_days, t_snow_max)

Write the shared forcing file.

The time axis is written in two forms, because the two readers accept different
ones and neither accepts both: a CF numeric axis for `chion_grid.x`, and
YYYY/MM/DD/HH integer variables for Chion.jl. Both encode the same instants, so
`infer_dt_days` and `times(2)-times(1)` agree on `dt_days`. See the long note at
the write site for why the CF axis alone is not enough.
"""
function write_forcing(path::AbstractString, scenarios::Vector{Scenario};
                       nstep::Int, dt_days::Float64=1.0, t_snow_max::Float64=273.15)
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

        for (name, data, units, long) in (
                ("TT", TT, "K", "Air temperature"),
                ("SF", SF, "kg m-2 s-1", "Snowfall rate"),
                ("RF", RF, "kg m-2 s-1", "Rainfall rate"),
                ("SWD", SWD, "W m-2", "Downward shortwave radiation"))
            # Declared ("x","y","time") in Julia order so the file carries
            # (time,y,x) in C order -- the layout chion_grid.x reads.
            vv = defVar(ds, name, Float64, ("x", "y", "time"))
            vv[:, :, :] = permutedims(data, (3, 2, 1))
            vv.attrib["units"] = units
            vv.attrib["long_name"] = long
        end

        for (name, val, units) in (("mask", 1.0, "1"),
                                   ("LAT", 70.0, "degrees_north"),
                                   ("SH", 1500.0, "m"))
            vv = defVar(ds, name, Float64, ("x", "y"))
            vv[:, :] = fill(val, nx, ny)
            vv.attrib["units"] = units
        end

        ds.attrib["title"] = "chion WP16 validation forcing"
        ds.attrib["scenarios"] = join((s.name for s in scenarios), ", ")
    end
    return (path=path, times=times, ncol=ncol, nstep=nstep, dt_days=dt_days)
end
