# ERA5 atmospheric forcing for the MOM6 + WW3 Gulf of Mexico / Caribbean cases

This document describes how hourly ERA5 data were turned into DATM forcing for the CESM case
`hurcn_gomcar_ocn_wav_era5.001`. That case is a copy of the JRA55-forced case
`hurcn_gomcar_ocn_wav.001`, used to simulate Hurricane Laura (19 Aug – 1 Sep 2020).

The approach feeds ERA5 through DATM's existing JRA55 forcing path (`CORE_IAF_JRA`). No CESM source
code changes or rebuild are needed.

## 1. Why reuse the JRA55 forcing path

In the JRA55 case, DATM runs in `CORE_IAF_JRA` mode
(`CESM/components/cdeps/datm/datm_datamode_jra_mod.F90`). That mode reads 8 surface variables, and
from them it derives everything the CMEPS mediator needs to compute bulk air–sea fluxes on the
ocean grid:

- air density, from pressure, temperature and humidity
- the rain/snow split, at the freezing temperature
- the four shortwave components: 0.28 direct visible, 0.31 diffuse visible, 0.24 direct near-infrared,
  0.17 diffuse near-infrared
- net shortwave, using the Large & Yeager (2008) ocean albedo

CDEPS also has a dedicated ERA5 mode, but it exports UFS-style fields. It lacks `Sa_u`, `Sa_tbot`,
`Sa_shum`, `Sa_pbot` and `Sa_dens`, so it can't be used with the CESM mediator's flux calculation.
Giving ERA5 the JRA55 variable names and units keeps all of the JRA55 processing.

## 2. Convert ERA5 to the JRA55 format

Script: `make_era5_datm.py`

```bash
module load conda/latest && conda activate npl
cd /glade/work/shanru/crocodile2026/era5_datm
python make_era5_datm.py --start 2020-08-17 --end 2020-09-03 \
       --outdir /glade/work/shanru/crocodile2026/era5_datm/hurcn_gomcar_2020
```

- **Source:** the ERA5 archive on Glade, `/glade/campaign/collections/gdex/data/d633000`
  (doi:10.5065/BH6N-5N20):
  - `e5.oper.an.sfc`: hourly analyses
  - `e5.oper.fc.sfc.meanflux`: hourly mean rates, stored by forecast start time and forecast hour
- **Region:** 255–325°E, 0–45°N at 0.25° (181 × 281 points), written with latitude running south
  to north. Set by `LON_RANGE` and `LAT_RANGE` in the script.
- **Period:** a few days wider than the run (432 hours), so DATM always has data on both sides of
  each model time for interpolation.

| JRA55 name | Units | ERA5 source | Handling |
|---|---|---|---|
| `u_10`, `v_10` | m/s | 10 m u and v wind (`VAR_10U`, `VAR_10V`) | copied as is |
| `t_10` | K | 2 m temperature (`VAR_2T`) | copied as is; DATM treats it as the reference-height temperature, like JRA55 |
| `q_10` | kg/kg | 2 m dewpoint (`VAR_2D`) + surface pressure (`SP`) | computed as below |
| `slp` | Pa | mean sea-level pressure (`MSL`) | copied as is |
| `prec` | kg/m²/s | mean total precipitation rate (`MTPR`) | time stamp moved 0.5 h earlier |
| `lwdn` | W/m² | mean downward longwave (`MSDWLWRF`) | time stamp moved 0.5 h earlier |
| `swdn` | W/m² | mean downward shortwave (`MSDWSWRF`) | time stamp moved 1 h earlier |

Specific humidity is computed from the 2 m dewpoint \(T_d\) (°C) and surface pressure \(p\) (Pa):

\[
e = 611.2\,\exp\!\left(\frac{17.67\,T_d}{T_d + 243.5}\right), \qquad
q = \frac{0.622\,e}{p - 0.378\,e}
\]

**Time stamps of the mean rates:** each ERA5 mean rate is an average over the hour ending at its
time stamp. The script moves the stamps so that DATM can read every file with `offset = 0`:

- `prec` and `lwdn` go to the middle of the averaging hour, which suits linear interpolation.
- `swdn` goes to the start of the averaging hour, which suits DATM's sun-angle (`coszen`)
  interpolation.

**Output:** one file per variable, `ERA5.<var>.20200817_20200903.nc`, with the variable under its
JRA55 name, 1-D `lat` and `lon`, and time in "days since 2020-01-01" on the Gregorian calendar.

## 3. Build the forcing-grid mesh

Script: `make_era5_mesh.py`

```bash
python make_era5_mesh.py hurcn_gomcar_2020/ERA5.u_10.20200817_20200903.nc \
       hurcn_gomcar_2020/ERA5_gomcar_025_ESMFmesh.nc
```

This writes an ESMF mesh file with 50,861 cells. Cell corners sit at the 0.25° half-points, corners
are listed SW, SE, NE, NW (1-based), and every cell is active. Cells are ordered latitude outer,
longitude inner, matching how CDEPS flattens a (lat, lon) data array. If the cell order and the data
files disagree, the forcing ends up in the wrong place without any error.

## 4. Create the ERA5 case

All CIME commands are run on Derecho, because the case belongs to that machine (`ssh derecho '...'`
from Casper).

```bash
cd /glade/work/shanru/crocodile2026/CESM/cime/scripts
./create_clone --case  /glade/work/shanru/crocodile2026/croc_cases/hurcn_gomcar_ocn_wav_era5.001 \
               --clone /glade/work/shanru/crocodile2026/croc_cases/hurcn_gomcar_ocn_wav.001 --keepexe

cd /glade/work/shanru/crocodile2026/croc_cases/hurcn_gomcar_ocn_wav_era5.001
./xmlchange RUNDIR=/glade/derecho/scratch/shanru/hurcn_gomcar_ocn_wav_era5.001/run
./xmlchange DOUT_S_ROOT=/glade/derecho/scratch/shanru/archive/hurcn_gomcar_ocn_wav_era5.001
MESH=/glade/work/shanru/crocodile2026/era5_datm/hurcn_gomcar_2020/ERA5_gomcar_025_ESMFmesh.nc
./xmlchange ATM_DOMAIN_MESH=$MESH,LND_DOMAIN_MESH=$MESH
```

- **`--keepexe`:** reuses the original case's build directory (`EXEROOT`), so nothing is rebuilt.
- **`RUNDIR` / `DOUT_S_ROOT`:** must be changed. The clone otherwise keeps the original case's
  directories and would overwrite its run.
- **`ATM_DOMAIN_MESH`:** puts the atmosphere on the native 0.25° ERA5 grid instead of the coarser
  JRA55 global grid (TL319, about 0.56°), so the hurricane's wind field isn't smoothed out. The
  mediator interpolates from this grid to the ocean grid.

### Point the forcing inputs at ERA5

In `user_nl_datm_streams`, override each of the 8 JRA55 forcing inputs, one block per input. For
example, the `U_10` input:

```
CORE_IAF_JRA_1p5_2023.NCEP.U_10:meshfile = /glade/work/shanru/crocodile2026/era5_datm/hurcn_gomcar_2020/ERA5_gomcar_025_ESMFmesh.nc
CORE_IAF_JRA_1p5_2023.NCEP.U_10:datafiles = /glade/work/shanru/crocodile2026/era5_datm/hurcn_gomcar_2020/ERA5.u_10.20200817_20200903.nc
CORE_IAF_JRA_1p5_2023.NCEP.U_10:year_first = 2020
CORE_IAF_JRA_1p5_2023.NCEP.U_10:year_last = 2020
CORE_IAF_JRA_1p5_2023.NCEP.U_10:year_align = 2020
CORE_IAF_JRA_1p5_2023.NCEP.U_10:taxmode = limit
CORE_IAF_JRA_1p5_2023.NCEP.U_10:offset = 0
```

The same block is repeated for the other 7 inputs, each pointing to its matching file
(`ERA5.v_10`, `t_10`, `q_10`, `slp`, `prec`, `swdn`, `lwdn`). `taxmode = limit` makes DATM stop
with an error instead of silently cycling the data if the model runs outside the data period.

Do not re-run CrocoDash's `configure_forcings` on this case: it would overwrite these overrides.

### Set up and submit

The `tutorial` queue (PBS reservation) isn't in CIME's list of queues, so switch to `main` for
`case.setup` and back afterwards:

```bash
./xmlchange JOB_QUEUE=main --force
./case.setup
./xmlchange JOB_QUEUE=tutorial --force
./preview_namelists
grep -c era5_datm CaseDocs/datm.streams.xml   # 16 matches: ERA5 mesh and data file for each of the 8 inputs
grep -c jra55     CaseDocs/datm.streams.xml   # should be 0
./case.submit
```

## 5. MOM6 time step fix

The first attempt stopped in MOM6 at the first time step, with "extreme surface values" and
"btstep eta dropped below bathyT". This happened at shallow, steep cells on the Bahamas banks, and
it hit the JRA55 case too. Halving the time steps in `user_nl_mom` of both cases fixed it:

```
DT = 300.0
DT_THERM = 600.0
```

The previous values were 900 s and 1800 s. Smoothing the topography would be a more lasting fix.

## 6. Wave (WW3) history output

The `.002` reruns of both cases (`hurcn_gomcar_ocn_wav.002`, `hurcn_gomcar_ocn_wav_era5.002`) add
hourly wave output:

```bash
./xmlchange HIST_OPTION=nhours,HIST_N=1
```

- **Why not `user_nl_ww3`:** WW3 takes its history frequency only from these run-wide (`ALLCOMP`)
  settings. `history_option` and `history_n` in `user_nl_ww3` are only written to the log.
- **Side effect:** the same settings make the mediator write hourly snapshots of all coupled fields
  (`*.cpl.hi.*`), about 240 MB per hour, or 70–100 GB for the 13-day run. These can be deleted if
  not needed.

## 7. Checks performed

- **Converted files:** values, units and time axes were checked against what DATM expects.
  Temperatures were checked in float64, because a float32 mean can give a misleading result.
- **Namelists:** `CaseDocs/datm.streams.xml` lists only ERA5 files and the ERA5 mesh.
- **Runs:** all four runs (`.001` and `.002` of both forcings) completed the 13 days and were
  archived under `/glade/derecho/scratch/shanru/archive/<case>/`.
- **Laura in the forcing:** compared with IBTrACS in `hurcn_gomcar_ocn_wav_diagnostics.ipynb`.

  | Source | Lowest pressure | Strongest 10 m wind | Mean track error |
  |---|---|---|---|
  | IBTrACS | 937 hPa | 66.9 m/s | — |
  | ERA5 | 978 hPa | 30.0 m/s | 51 km |
  | JRA55 | 995 hPa | 23.3 m/s | 82 km |

  Both forcings track Laura well but are far too weak at its peak. ERA5 is closer than JRA55.

## 8. Reusing the workflow for another storm or period

1. **New forcing files:** set `LON_RANGE` and `LAT_RANGE` in `make_era5_datm.py`, then run it with
   new `--start`, `--end` and `--outdir` values. Leave a margin of 1–2 days on each side of the run.
2. **New mesh:** run `make_era5_mesh.py` on one of the new files. Only needed if the region changed.
3. **Case settings:**
   - In `user_nl_datm_streams`: update `datafiles`, and `year_first`, `year_last` and `year_align`
     if the year changed.
   - With `xmlchange`: update `ATM_DOMAIN_MESH` and `LND_DOMAIN_MESH` if the mesh changed, and
     `RUN_STARTDATE`, `STOP_OPTION` and `STOP_N` for the new run period.
4. **Check before running:** run `./preview_namelists`, then confirm `CaseDocs/datm.streams.xml`
   lists the new files.

## Files

| Path | Contents |
|---|---|
| `/glade/work/shanru/crocodile2026/era5_datm/make_era5_datm.py` | ERA5 to JRA55-format conversion |
| `/glade/work/shanru/crocodile2026/era5_datm/make_era5_mesh.py` | ESMF mesh writer |
| `/glade/work/shanru/crocodile2026/era5_datm/hurcn_gomcar_2020/` | 8 forcing files and the mesh |
| `/glade/work/shanru/crocodile2026/croc_cases/hurcn_gomcar_ocn_wav_era5.001/` | ERA5 case |
| `/glade/work/shanru/crocodile2026/croc_cases/hurcn_gomcar_ocn_wav_era5.002/` | ERA5 case with hourly WW3 output |
| `/glade/work/shanru/crocodile2026/workspace/hurcn_gomcar_ocn_wav_diagnostics.ipynb` | SST, forcing and wave comparison |
