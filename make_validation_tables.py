import pathlib as pl
import xarray as xr
import pandas as pd

base_dir = pl.Path("simulation")

simulation_dirs = base_dir.glob("validation_*")
simulation_dirs = sorted(simulation_dirs)

simulation_dir = simulation_dirs[0]
for simulation_dir in simulation_dirs:
    print(f"simulation: {simulation_dir.name}")

    station_files = simulation_dir.glob("station_*.csv")
    station_files = sorted(station_files)
    station_file = station_files[0]

    station = pd.read_csv(station_file, index_col=0)
    lat = station.loc["adjusted_lat"].iloc[0]
    lon = station.loc["adjusted_long"].iloc[0]
    lat = float(lat)
    lon = float(lon)

    observed_file = simulation_dir / "discharge.csv"

    observed = pd.read_csv(observed_file, index_col=0, parse_dates=True)

    reference_dir = simulation_dir / "reference"
    commit_dirs = reference_dir.iterdir()
    commit_dirs = sorted(commit_dirs)

    commit_dir = commit_dirs[0]
    for commit_dir in commit_dirs:
        print(f"commit: {commit_dir.name}")

        netcdf_dir = commit_dir / "netcdf"
        simulated_file = netcdf_dir / "discharge_dailyTot_output.nc"
        with xr.open_dataarray(simulated_file) as simulated:
            simulated = simulated.sel(
                lat=lat, lon=lon, method="nearest", drop=True
            )

        simulated = simulated.to_dataframe()
        simulated = simulated.rename(columns={"discharge": "simulated"})
        simulated.index.name = "date"

        table = pd.merge(
            simulated, observed, left_index=True, right_index=True, how="outer"
        )
        table = table.dropna(subset=["simulated"])
        table = table.sort_index()

        table_dir = commit_dir / "table"
        table_out = table_dir / "discharge.csv"
        table_out.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(table_out)
