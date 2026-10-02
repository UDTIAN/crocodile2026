#!/usr/bin/env python3
"""Write an ESMF unstructured mesh for a regular lat/lon forcing grid.

Element order is (lat outer, lon inner), matching a (lat, lon) data array
flattened in C order, which is how CDEPS reads stream data onto the mesh.

Usage:
  module load conda/latest && conda activate npl
  python make_era5_mesh.py ERA5.u_10.20200817_20200903.nc ERA5_gomcar_025_ESMFmesh.nc
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr


def edges(centers: np.ndarray) -> np.ndarray:
    mid = 0.5 * (centers[1:] + centers[:-1])
    return np.concatenate([[centers[0] - (mid[0] - centers[0])], mid, [centers[-1] + (centers[-1] - mid[-1])]])


def build_mesh(lat: np.ndarray, lon: np.ndarray) -> xr.Dataset:
    nlat, nlon = lat.size, lon.size
    lat_e = np.clip(edges(lat), -90.0, 90.0)
    lon_e = edges(lon)

    node_lon, node_lat = np.meshgrid(lon_e, lat_e)
    node_coords = np.stack([node_lon.ravel(), node_lat.ravel()], axis=1)

    nid = np.arange((nlat + 1) * (nlon + 1)).reshape(nlat + 1, nlon + 1) + 1
    sw = nid[:-1, :-1].ravel()
    se = nid[:-1, 1:].ravel()
    ne = nid[1:, 1:].ravel()
    nw = nid[1:, :-1].ravel()
    conn = np.stack([sw, se, ne, nw], axis=1).astype(np.int32)

    c_lon, c_lat = np.meshgrid(lon, lat)
    center_coords = np.stack([c_lon.ravel(), c_lat.ravel()], axis=1)

    nelem = nlat * nlon
    return xr.Dataset(
        {
            "nodeCoords": (("nodeCount", "coordDim"), node_coords, {"units": "degrees"}),
            "elementConn": (
                ("elementCount", "maxNodePElement"),
                conn,
                {"long_name": "Node indices that define the element connectivity", "start_index": 1},
            ),
            "numElementConn": (("elementCount",), np.full(nelem, 4, dtype=np.int8)),
            "centerCoords": (("elementCount", "coordDim"), center_coords, {"units": "degrees"}),
            "elementMask": (("elementCount",), np.ones(nelem, dtype=np.int32)),
        },
        attrs={"gridType": "unstructured mesh", "title": f"regular lat/lon mesh {nlat}x{nlon}"},
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("datafile", type=Path, help="NetCDF file with 1-D lat and lon")
    p.add_argument("meshfile", type=Path, help="Output ESMF mesh file")
    args = p.parse_args()

    with xr.open_dataset(args.datafile) as ds:
        lat = ds["lat"].values.astype("float64")
        lon = ds["lon"].values.astype("float64")
    if not (np.all(np.diff(lat) > 0) and np.all(np.diff(lon) > 0)):
        raise ValueError("lat and lon must both be strictly increasing")

    mesh = build_mesh(lat, lon)
    encoding = {v: {"_FillValue": None} for v in mesh.data_vars}
    mesh.to_netcdf(args.meshfile, format="NETCDF4", encoding=encoding)
    print(f"wrote {args.meshfile}: {mesh.sizes['elementCount']} elements ({lat.size} lat x {lon.size} lon)")


if __name__ == "__main__":
    main()
