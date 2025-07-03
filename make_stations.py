import pathlib as pl
import pandas as pd

input_dir = pl.Path("input/grdc")
save_dir = pl.Path("saves/stations")

grdc_file = input_dir / "GRDC_stations.xlsx"
grdc = pd.read_excel(grdc_file, index_col=0)

station_infos = {
    "Rhine": (6435060, pd.Timestamp("2000-01-01"), pd.Timestamp("2010-12-31")),
    "Colombia": (
        4115201,
        pd.Timestamp("1992-01-01"),
        pd.Timestamp("2002-12-31"),
    ),
    "Orinoco": (
        3206720,
        pd.Timestamp("1979-01-01"),
        pd.Timestamp("1989-12-31"),
    ),
    "Limpopo": (
        1896502,
        pd.Timestamp("1983-01-01"),
        pd.Timestamp("1993-12-31"),
    ),
    "Mekong": (
        2469260,
        pd.Timestamp("1983-01-01"),
        pd.Timestamp("1993-12-31"),
    ),
    "Murray": (
        5404271,
        pd.Timestamp("1989-01-01"),
        pd.Timestamp("1999-12-31"),
    ),
    "Kolyma": (
        2998510,
        pd.Timestamp("1998-01-01"),
        pd.Timestamp("2008-12-31"),
    ),
    "Mackenzie": (
        4208025,
        pd.Timestamp("2000-01-01"),
        pd.Timestamp("2010-12-31"),
    ),
}

for name, (grdc_id, start_date, end_date) in station_infos.items():
    print(f"name: {name}")
    station = grdc.loc[grdc_id].copy()
    station["name"] = name
    station["start_date"] = start_date
    station["end_date"] = end_date

    stations_out = save_dir / f"station_{name}.csv"
    stations_out.parent.mkdir(parents=True, exist_ok=True)
    station.to_csv(stations_out)
