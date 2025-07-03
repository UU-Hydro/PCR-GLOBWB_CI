import pathlib as pl
import pandas as pd

input_dir = pl.Path("input")
base_dir = pl.Path("simulation")

grdc_dir = input_dir / "grdc"
data_dir = grdc_dir / "data"

simulation_dirs = base_dir.glob("validation_*")
simulation_dirs = sorted(simulation_dirs)

simulation_dir = simulation_dirs[0]
for simulation_dir in simulation_dirs:
    print(f"simulation: {simulation_dir.name}")

    station_files = simulation_dir.glob("station_*.csv")
    station_files = sorted(station_files)
    station_file = station_files[0]

    station = pd.read_csv(station_file, index_col=0)
    index = station.columns[0]

    data_out = simulation_dir / "discharge.csv"

    data_file = data_dir / f"{index}_Q_Day.Cmd.txt"
    data = pd.read_csv(
        data_file,
        delimiter=";",
        encoding="ISO-8859-1",
        comment="#",
    )

    data.columns = ["date", None, "observed"]
    data = data[["date", "observed"]]
    data = data.set_index("date")

    data.loc[data["observed"] < 0, "observed"] = pd.NA
    data = data.dropna()

    data_out.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(data_out)
