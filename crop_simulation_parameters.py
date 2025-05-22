import pathlib as pl
import datetime as dt
import pcraster as pcr
import xarray as xr
import numpy as np
import pandas as pd


def store_pcraster(pcr_map: pcr.Field,
                   file: pl.Path,) -> None:
    file_tmp = file.with_suffix('.tmp.nc')
    file_tmp.parent.mkdir(parents=True, exist_ok=True)
    pcr.report(pcr_map, str(file_tmp))
    file_tmp.rename(file)


def store_netcdf(nc_ds: xr.Dataset,
                 file: pl.Path) -> None:
    file_tmp = file.with_suffix('.tmp.nc')
    file_tmp.parent.mkdir(parents=True, exist_ok=True)
    nc_ds.to_netcdf(file_tmp)
    file_tmp.rename(file)


def crop_pcraster_file(file: pl.Path,
                       domain_west: float,
                       domain_north: float,
                       domain_east: float,
                       domain_south: float,
                       file_out: pl.Path) -> None:

    pcr.setclone(str(file))
    nrRows = pcr.clone().nrRows()
    nrCols = pcr.clone().nrCols()
    cellSize = round(pcr.clone().cellSize(), 5)
    west = round(pcr.clone().west())
    north = round(pcr.clone().north())
    east = round(west + (cellSize * nrCols))
    south = round(north - (cellSize * nrRows))

    lats = np.linspace(north - cellSize/2,
                       south + cellSize/2,
                       nrRows)
    lons = np.linspace(west + cellSize/2,
                       east - cellSize/2,
                       nrCols)

    domain_lat_sel = (lats >= domain_south) & (lats <= domain_north)
    domain_lon_sel = (lons >= domain_west) & (lons <= domain_east)
    domain_lats = lats[domain_lat_sel]
    domain_lons = lons[domain_lon_sel]
    domain_nrRows = len(domain_lats)
    domain_nrCols = len(domain_lons)

    pcr_map = pcr.readmap(str(file))
    data_type = pcr_map.dataType()
    if data_type == pcr.Scalar:
        missing_value = -9999
    elif data_type == pcr.Nominal:
        missing_value = -1
    elif data_type == pcr.Ordinal:
        missing_value = -1
    elif data_type == pcr.Boolean:
        missing_value = 255
    elif data_type == pcr.Ldd:
        missing_value = 255
    else:
        raise ValueError(f"Unknown data type: {data_type}")
    pcr_map = pcr.pcr2numpy(pcr_map, missing_value)

    domain_map = pcr_map[domain_lat_sel, :][:, domain_lon_sel]
    pcr.setclone(domain_nrRows,
                 domain_nrCols,
                 cellSize,
                 domain_west,
                 domain_north)
    domain_map = pcr.numpy2pcr(data_type, domain_map, missing_value)

    store_pcraster(pcr_map=domain_map,
                   file=file_out)


def crop_netcdf_file(file: pl.Path,
                     domain_west: float,
                     domain_north: float,
                     domain_east: float,
                     domain_south: float,
                     domain_start: pd.Timestamp,
                     domain_end: pd.Timestamp,
                     file_out: pl.Path) -> None:

    with xr.open_dataset(file) as ds:
        pass

    lat_name = 'lat'
    lon_name = 'lon'
    if lat_name not in ds.dims:
        lat_name = 'latitude'
    if lon_name not in ds.dims:
        lon_name = 'longitude'

    domain_ds = ds.sel({lat_name: slice(domain_north, domain_south),
                        lon_name: slice(domain_west, domain_east)})

    if 'time' not in domain_ds.dims:
        store_netcdf(nc_ds=domain_ds, file=file_out)
        return

    times = ds.time
    domain_times = pd.Series(pd.date_range(start=domain_start,
                                           end=domain_end,
                                           freq='D'))
    if len(times) == 1:
        store_netcdf(nc_ds=domain_ds, file=file_out)
        return

    doys = times.dt.dayofyear
    udoys = np.unique(doys.values)
    domain_doys = domain_times.dt.dayofyear
    domain_udoys = np.unique(domain_doys.values)
    if len(doys) == len(udoys):
        domain_ds = domain_ds.sel(time=doys.isin(domain_udoys))
        store_netcdf(nc_ds=domain_ds, file=file_out)
        return

    start = times[0]
    end = times[-1]
    if domain_end < start:
        domain_end.year = start.year
    if domain_start > end:
        domain_start.year = end.year
    domain_times = pd.Series(pd.date_range(start=domain_start,
                                           end=domain_end,
                                           freq='D'))
    domain_ds = domain_ds.sel(time=times.isin(domain_times))
    store_netcdf(nc_ds=domain_ds, file=file_out)


simulation_dir = pl.Path('simulation')
parameter_dir = pl.Path('parameters')
levels = [0.0, 0.01, 0.05, 0.10, 0.20, 0.30,
          0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]

simulations = simulation_dir.iterdir()
simulations = sorted(simulations)

simulation = next(iter(simulations))
for simulation in simulations:
    print(f'simulation: {simulation}')

    configuration_file = simulation / 'configuration.ini'
    out_dir = simulation / 'parameters'
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(configuration_file, 'r') as f:
        lines = f.readlines()

    domain_file = None
    start = None
    end = None
    opendap = False
    for line in lines:
        if 'inputDir' in line:
            opendap = 'opendap.4tu.nl/' in line
            continue
        if 'cloneMap' in line:
            domain_file = line.split('=')[1].strip()
            domain_file = pl.Path(domain_file)
            domain_file = simulation / domain_file.name
            continue
        if 'startTime' in line:
            start = line.split('=')[1].strip()
            start = dt.datetime.strptime(start, '%Y-%m-%d').date()
            start = pd.Timestamp(start)
            continue
        if 'endTime' in line:
            end = line.split('=')[1].strip()
            end = dt.datetime.strptime(end, '%Y-%m-%d').date()
            end = pd.Timestamp(end)
            continue
    if start is None or end is None:
        raise ValueError('Start or end time not found in configuration file.')
    if domain_file is None:
        raise ValueError('Domain file not found in configuration file.')

    if opendap:
        continue

    pcr.setclone(str(domain_file))
    nrRows = pcr.clone().nrRows()
    nrCols = pcr.clone().nrCols()
    cellSize = round(pcr.clone().cellSize(), 5)
    west = round(pcr.clone().west())
    north = round(pcr.clone().north())
    east = round(west + (cellSize * nrCols))
    south = round(north - (cellSize * nrRows))

    pcr_files = []
    for index, line in enumerate(lines):
        if '.map' not in line:
            continue
        if 'cloneMap' in line:
            continue
        if 'landmask' in line:
            continue
        pcr_file = line.split('=')[1].strip()
        pcr_file = pl.Path(pcr_file)
        if '%04d' in pcr_file.name:
            for level in levels:
                prc_name_level = pcr_file.name.replace('%04d',
                                                       f'{level*100:04.0f}')
                pcr_file_level = pcr_file.parent / prc_name_level
                pcr_files.append(pcr_file_level)
        else:
            pcr_files.append(pcr_file)

    nc_files = []
    for index, line in enumerate(lines):
        if '.nc' not in line:
            continue
        nc_file = line.split('=')[1].strip()
        nc_file = pl.Path(nc_file)
        if '%04d' in nc_file.name:
            for level in levels:
                prc_name_level = nc_file.name.replace('%04d',
                                                      f'{level*100:04.0f}')
                nc_file_level = nc_file.parent / prc_name_level
                nc_files.append(nc_file_level)
        else:
            nc_files.append(nc_file)

    for pcr_file in pcr_files:
        pcr_out = out_dir / pcr_file
        if pcr_out.exists():
            continue
        crop_pcraster_file(file=parameter_dir / pcr_file,
                           domain_west=west,
                           domain_north=north,
                           domain_east=east,
                           domain_south=south,
                           file_out=pcr_out)

    for nc_file in nc_files:
        nc_out = out_dir / nc_file
        if nc_out.exists():
            continue
        crop_netcdf_file(file=parameter_dir / nc_file,
                         domain_west=west,
                         domain_north=north,
                         domain_east=east,
                         domain_south=south,
                         domain_start=start,
                         domain_end=end,
                         file_out=nc_out)
