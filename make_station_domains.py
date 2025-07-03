import pathlib as pl
import pcraster as pcr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

input_dir = pl.Path("input")
save_dir = pl.Path("saves")

ldd_dir = input_dir / "ldd"
station_dir = save_dir / "stations"

ldd_info = {
    "30min": ldd_dir / "lddsound_30min.map",
    "05min": ldd_dir / "lddsound_05min.map",
    "30sec": ldd_dir / "lddsound_30sec_version_202005XX.map",
}
out_dir = pl.Path("saves/domains")

for resolution, ldd_file in ldd_info.items():
    print(f"resolution: {resolution}")

    pcr.setclone(str(ldd_file))
    nrRows = pcr.clone().nrRows()
    nrCols = pcr.clone().nrCols()
    cellSize = round(pcr.clone().cellSize(), 5)
    west = round(pcr.clone().west())
    north = round(pcr.clone().north())
    east = round(west + (cellSize * nrCols))
    south = round(north - (cellSize * nrRows))

    lats = np.linspace(north - cellSize / 2, south + cellSize / 2, nrRows)
    lons = np.linspace(west + cellSize / 2, east - cellSize / 2, nrCols)

    ldd = pcr.readmap(str(ldd_file))

    station_files = station_dir.glob(f"station_*_{resolution}.csv")
    station_files = sorted(station_files)

    if len(station_files) == 0:
        continue

    for station_file in station_files:
        print(f"stations_file: {station_file}")

        station = pd.read_csv(station_file, index_col=0)
        name = station.loc["name"].iloc[0]
        point_lat = station.loc["adjusted_lat"].iloc[0]
        point_lon = station.loc["adjusted_long"].iloc[0]
        point_lat = float(point_lat)
        point_lon = float(point_lon)
        point_y = np.argmin(np.abs(lats - point_lat))
        point_x = np.argmin(np.abs(lons - point_lon))

        pit = np.zeros((nrRows, nrCols), dtype=np.int32)
        pit[point_y, point_x] = 1

        pcr.setclone(str(ldd_file))
        pit = pcr.numpy2pcr(pcr.Nominal, pit, 0)
        subcatchment = pcr.catchment(ldd, pit)
        subcatchment = pcr.pcr2numpy(subcatchment, 0)

        mask = subcatchment == 1
        mask_lat_sel = np.any(mask, axis=1)
        mask_lon_sel = np.any(mask, axis=0)
        mask_lats = lats[mask_lat_sel]
        mask_lons = lons[mask_lon_sel]

        domain_north = np.ceil(mask_lats[0])
        domain_south = np.floor(mask_lats[-1])
        domain_west = np.floor(mask_lons[0])
        domain_east = np.ceil(mask_lons[-1])

        domain_lat_sel = (lats >= domain_south) & (lats <= domain_north)
        domain_lon_sel = (lons >= domain_west) & (lons <= domain_east)

        domain_lats = lats[domain_lat_sel]
        domain_lons = lons[domain_lon_sel]
        domain_nrRows = len(domain_lats)
        domain_nrCols = len(domain_lons)

        domain_mask = mask[domain_lat_sel, :][:, domain_lon_sel].copy()

        domain_plot = mask.astype(np.int32)
        domain_plot[point_y, point_x] = 2
        domain_plot = domain_plot[domain_lat_sel, :][:, domain_lon_sel].copy()
        plt.imshow(domain_plot)
        plt.title(f"domain_mask: {name}")
        plt.colorbar()
        plt.show()

        pcr.setclone(
            domain_nrRows, domain_nrCols, cellSize, domain_west, domain_north
        )
        domain_mask = pcr.numpy2pcr(pcr.Boolean, domain_mask, False)

        mask_out = out_dir / f"domain_{name}_{resolution}_station.map"
        mask_out.parent.mkdir(parents=True, exist_ok=True)
        pcr.report(domain_mask, str(mask_out))
        print(f"mask_out: {mask_out}")
