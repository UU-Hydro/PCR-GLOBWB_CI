import pathlib as pl
import numpy as np
import pandas as pd
import pcraster as pcr

input_dir = pl.Path("input")
save_dir = pl.Path("saves")
distance = 0.5  # degrees

input_files = {
    "30min": ("lddsound_30min.map", "cellarea30min.map"),
    "05min": ("lddsound_05min.map", "cellsize05min.correct.map"),
}

station_dir = save_dir / "stations"
ldd_dir = input_dir / "ldd"
area_dir = input_dir / "area"

station_files = station_dir.glob("station_*.csv")
station_files = sorted(station_files)

resolution = "05min"
ldd_file, area_file = input_files[resolution]
for resolution, (ldd_file, area_file) in input_files.items():
    print(f"resolution: {resolution}")

    ldd_file = ldd_dir / ldd_file
    area_file = area_dir / area_file

    pcr.setclone(str(ldd_file))
    nrRows = pcr.clone().nrRows()
    nrCols = pcr.clone().nrCols()
    cellSize = pcr.clone().cellSize()
    lats = np.linspace(90 - cellSize / 2, -90 + cellSize / 2, nrRows)
    lons = np.linspace(-180 + cellSize / 2, 180 - cellSize / 2, nrCols)
    cellDist = int(distance // cellSize)

    ldd = pcr.readmap(str(ldd_file))
    area = pcr.readmap(str(area_file))
    accumulation = pcr.catchmenttotal(area, ldd)
    accumulation = pcr.pcr2numpy(accumulation, np.nan)

    for station_file in station_files:
        if station_file.stem.endswith(("_30min", "_05min", "_30sec")):
            continue
        print(f"stations_file: {station_file}")

        station = pd.read_csv(station_file, index_col=0)
        name = station.loc["name"].iloc[0]
        lat = station.loc["lat"].iloc[0]
        lon = station.loc["long"].iloc[0]
        area = station.loc["area"].iloc[0]
        lat = float(lat)
        lon = float(lon)
        area = float(area)

        y = np.argmin(np.abs(lats - lat))
        x = np.argmin(np.abs(lons - lon))
        station.loc["y"] = y
        station.loc["x"] = x

        acc = accumulation[
            (y - cellDist) : (y + cellDist + 1),
            (x - cellDist) : (x + cellDist + 1),
        ]
        acc = acc * 1e-6

        diff = np.abs(acc - area)
        diff[np.isnan(diff)] = np.inf
        i = np.argmin(diff)
        yoffset = i // (2 * cellDist + 1) - cellDist
        xoffset = i % (2 * cellDist + 1) - cellDist

        adjusted_y = y + yoffset
        adjusted_x = x + xoffset
        adjusted_lat = lats[adjusted_y]
        adjusted_long = lons[adjusted_x]
        station.loc["adjusted_y"] = adjusted_y
        station.loc["adjusted_x"] = adjusted_x
        station.loc["adjusted_lat"] = adjusted_lat
        station.loc["adjusted_long"] = adjusted_long

        station_out = (
            station_file.parent / f"{station_file.stem}_{resolution}.csv"
        )
        station_out.parent.mkdir(parents=True, exist_ok=True)
        station.to_csv(station_out)
