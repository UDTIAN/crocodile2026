#!/usr/bin/env python3
"""Build JRA-format DATM forcing files from the NCAR GDEX ERA5 archive (d633000).

The output files are read by DATM in CORE_IAF_JRA mode, which does no unit
conversion and exports exactly the fields the CMEPS mediator needs to compute
bulk air-sea fluxes. Each output file holds one variable named with its JRA
name and in JRA units:

    u_10, v_10  [m/s]      <- 10u, 10v                     (analysis)
    t_10        [K]        <- 2t                           (analysis, 2 m)
    q_10        [kg/kg]    <- 2d + sp                      (analysis, 2 m)
    slp         [Pa]       <- msl                          (analysis)
    prec        [kg/m2/s]  <- mtpr                         (forecast mean rate)
    swdn        [W/m2]     <- msdwswrf                     (forecast mean rate)
    lwdn        [W/m2]     <- msdwlwrf                     (forecast mean rate)

Forecast mean rates are valid over the hour ending at init+forecast_hour.
They are re-stamped so that CDEPS interpolates them correctly with offset=0:
  prec, lwdn : middle of the averaging hour (linear interpolation)
  swdn       : start of the averaging hour  (coszen interpolation)

Latitude is written south-to-north; the ESMF mesh from make_era5_mesh.py uses
the same (lat outer, lon inner) element ordering.

Usage:
  module load conda/latest && conda activate npl
  python make_era5_datm.py
  python make_era5_datm.py --start 2020-08-17 --end 2020-09-03 --outdir /path
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


ERA5_ROOT = Path("/glade/campaign/collections/gdex/data/d633000")
DEFAULT_OUTDIR = Path("/glade/work/shanru/crocodile2026/era5_datm/hurcn_gomcar_2020")

LON_RANGE = (255.0, 325.0)
LAT_RANGE = (0.0, 45.0)

ANALYSIS = {
    "u_10": ("128_165_10u", "VAR_10U"),
    "v_10": ("128_166_10v", "VAR_10V"),
    "t_10": ("128_167_2t", "VAR_2T"),
    "slp": ("128_151_msl", "MSL"),
}
HUMIDITY_INPUTS = {"d2m": ("128_168_2d", "VAR_2D"), "sp": ("128_134_sp", "SP")}
MEANFLUX = {
    "prec": ("235_055_mtpr", "MTPR", -0.5),
    "lwdn": ("235_036_msdwlwrf", "MSDWLWRF", -0.5),
    "swdn": ("235_035_msdwswrf", "MSDWSWRF", -1.0),
}

ATTRS = {
    "u_10": ("10 m eastward wind", "m/s"),
    "v_10": ("10 m northward wind", "m/s"),
    "t_10": ("2 m air temperature (used as 10 m)", "K"),
    "q_10": ("2 m specific humidity (used as 10 m)", "kg/kg"),
    "slp": ("mean sea level pressure", "Pa"),
    "prec": ("total precipitation rate", "kg/m2/s"),
    "swdn": ("surface downward shortwave radiation", "W/m2"),
    "lwdn": ("surface downward longwave radiation", "W/m2"),
}

TIME_UNITS = "days since 2020-01-01 00:00:00"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--start", default="2020-08-17", help="First date (YYYY-MM-DD)")
    p.add_argument("--end", default="2020-09-03", help="Last date, inclusive (YYYY-MM-DD)")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return p.parse_args()


def months_between(start: pd.Timestamp, end: pd.Timestamp) -> list[str]:
    first = (start - pd.Timedelta(days=1)).to_period("M")
    last = end.to_period("M")
    return [str(p).replace("-", "") for p in pd.period_range(first, last, freq="M")]


def subset(da: xr.DataArray) -> xr.DataArray:
    da = da.sel(longitude=slice(*LON_RANGE))
    da = da.sel(latitude=slice(LAT_RANGE[1], LAT_RANGE[0]))
    return da.isel(latitude=slice(None, None, -1))


def open_analysis(code: str, var: str, months: list[str]) -> xr.DataArray:
    files = []
    for m in months:
        files += sorted((ERA5_ROOT / "e5.oper.an.sfc" / m).glob(f"e5.oper.an.sfc.{code}.ll025sc.*.nc"))
    if not files:
        raise FileNotFoundError(f"No analysis files for {code} in {months}")
    ds = xr.open_mfdataset(files, combine="by_coords", data_vars="minimal", coords="minimal", compat="override")
    return subset(ds[var])


def open_meanflux(code: str, var: str, months: list[str], shift_hours: float) -> xr.DataArray:
    files = []
    for m in months:
        files += sorted((ERA5_ROOT / "e5.oper.fc.sfc.meanflux" / m).glob(f"e5.oper.fc.sfc.meanflux.{code}.ll025sc.*.nc"))
    if not files:
        raise FileNotFoundError(f"No meanflux files for {code} in {months}")

    pieces = []
    for f in files:
        ds = xr.open_dataset(f)
        da = subset(ds[var])
        init = pd.to_datetime(ds["forecast_initial_time"].values)
        fhr = ds["forecast_hour"].values.astype(float)
        valid = (init.values[:, None] + (fhr[None, :] * 3600e9).astype("timedelta64[ns]")).ravel()
        da = da.stack(time=("forecast_initial_time", "forecast_hour"))
        da = da.drop_vars(["time", "forecast_initial_time", "forecast_hour"])
        da = da.assign_coords(time=pd.DatetimeIndex(valid) + pd.Timedelta(hours=shift_hours))
        da = da.transpose("time", "latitude", "longitude")
        pieces.append(da)

    da = xr.concat(pieces, dim="time").sortby("time")
    _, keep = np.unique(da["time"].values, return_index=True)
    return da.isel(time=np.sort(keep))


def specific_humidity(d2m: xr.DataArray, sp: xr.DataArray) -> xr.DataArray:
    td_c = d2m - 273.15
    e = 611.2 * np.exp(17.67 * td_c / (td_c + 243.5))
    return 0.622 * e / (sp - 0.378 * e)


def write(name: str, da: xr.DataArray, start: pd.Timestamp, end: pd.Timestamp, outdir: Path) -> Path:
    da = da.sel(time=slice(start, end + pd.Timedelta(hours=23, minutes=59)))
    times = pd.DatetimeIndex(da["time"].values)
    gaps = np.diff(times.values).astype("timedelta64[m]").astype(int)
    if len(gaps) and (gaps.min() != 60 or gaps.max() != 60):
        raise RuntimeError(f"{name}: time axis is not continuous hourly (gaps {gaps.min()}-{gaps.max()} min)")

    long_name, units = ATTRS[name]
    out = xr.Dataset(
        {name: (("time", "lat", "lon"), da.values.astype("float32"), {"long_name": long_name, "units": units})},
        coords={
            "time": ("time", times),
            "lat": ("lat", da["latitude"].values, {"units": "degrees_north", "long_name": "latitude"}),
            "lon": ("lon", da["longitude"].values, {"units": "degrees_east", "long_name": "longitude"}),
        },
        attrs={
            "title": f"ERA5 {name} in JRA DATM format",
            "source": "NCAR GDEX d633000 (ERA5), https://doi.org/10.5065/BH6N-5N20",
            "history": "created by make_era5_datm.py",
        },
    )
    fname = outdir / f"ERA5.{name}.{start:%Y%m%d}_{end:%Y%m%d}.nc"
    out.to_netcdf(
        fname,
        format="NETCDF4",
        encoding={
            name: {"zlib": True, "complevel": 1, "_FillValue": None},
            "time": {"units": TIME_UNITS, "calendar": "gregorian", "dtype": "float64"},
            "lat": {"_FillValue": None},
            "lon": {"_FillValue": None},
        },
    )
    print(f"  wrote {fname}  ({out.sizes['time']} times, {times[0]} to {times[-1]})")
    return fname


def main() -> None:
    args = parse_args()
    start = pd.Timestamp(args.start)
    end = pd.Timestamp(args.end)
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    months = months_between(start, end)
    print(f"ERA5 {start:%Y-%m-%d} to {end:%Y-%m-%d}, months {months}, output {outdir}")

    for name, (code, var) in ANALYSIS.items():
        print(f"{name} <- {code}")
        write(name, open_analysis(code, var, months).load(), start, end, outdir)

    print("q_10 <- 2d + sp")
    d2m = open_analysis(*HUMIDITY_INPUTS["d2m"], months).load()
    sp = open_analysis(*HUMIDITY_INPUTS["sp"], months).load()
    write("q_10", specific_humidity(d2m, sp), start, end, outdir)

    for name, (code, var, shift) in MEANFLUX.items():
        print(f"{name} <- {code} (time shift {shift:+} h)")
        write(name, open_meanflux(code, var, months, shift).load(), start, end, outdir)

    print("Done.")


if __name__ == "__main__":
    main()
