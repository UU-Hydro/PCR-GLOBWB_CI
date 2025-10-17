import pathlib as pl
import pcraster as pcr
import xarray as xr
import numpy as np
import pandas as pd

# landsurface influx: precipitation,irrGrossDemand,surfaceWaterInf
# landsurface stores: interceptStor,snowFreeWater,snowCoverSWE,topWaterLayer,storUppTotal,storLowTotal,storGroundwater
# landsurface outflux: actualET,runoff,nonFossilGroundwaterAbstraction
# routing influx: nonIrrReturnFlow,runoff
# routing stores: surfaceWaterStorage
# routing outflux: waterBodyActEvaporation,discharge,surfaceWaterInf,surfaceWaterAbstraction
water_balance_components = {
    "landsurface": {
        "influx": ["precipitation", "irrGrossDemand", "surfaceWaterInf"],
        "stores": [
            "interceptStor",
            "snowFreeWater",
            "snowCoverSWE",
            "topWaterLayer",
            "storUppTotal",
            "storLowTotal",
            "storGroundwater",
        ],
        "outflux": ["actualET", "runoff", "nonFossilGroundwaterAbstraction"],
    },
    "routing": {
        "influx": ["nonIrrReturnFlow", "runoff"],
        "stores": ["surfaceWaterStorage"],
        "outflux": [
            "waterBodyActEvaporation",
            "discharge",
            "surfaceWaterInf",
            "surfaceWaterAbstraction",
        ],
    },
}

base_dir = pl.Path("simulation")

simulation_dirs = base_dir.glob("validation_*")
simulation_dirs = sorted(simulation_dirs)

simulation_dir = simulation_dirs[0]
for simulation_dir in simulation_dirs:
    print(f"simulation: {simulation_dir.name}")

    reference_dir = simulation_dir / "reference"
    commit_dirs = reference_dir.iterdir()
    commit_dirs = sorted(commit_dirs)

    commit_dir = commit_dirs[0]
    for commit_dir in commit_dirs:
        print(f"commit: {commit_dir.name}")

        log_dir = commit_dir / "log"
        configuration_files = log_dir.glob("*.ini")
        configuration_file = next(configuration_files, None)
        if configuration_file is None:
            raise FileNotFoundError(
                f"No configuration file found in {log_dir}"
            )

        with open(configuration_file) as f:
            configuration = f.read()
        input_dir = configuration.split("inputDir")[1]
        input_dir = input_dir.split("=")[1]
        input_dir = input_dir.split("\n")[0].strip()
        area_file = configuration.split("cellAreaMap")[1]
        area_file = area_file.split("=")[1]
        area_file = area_file.split("\n")[0].strip()

        area_file = pl.Path(input_dir) / area_file
        pcr.setclone(str(area_file))
        area = pcr.readmap(str(area_file))
        area = pcr.pcr2numpy(area, np.nan)

        tables = []

        netcdf_dir = commit_dir / "netcdf"

        system, components = next(iter(water_balance_components.items()))
        component, variables = next(iter(components.items()))
        variable = variables[0]

        system = "routing"
        component = "stores"
        variable = "surfaceWaterStorage"

        for system, components in water_balance_components.items():
            for component, variables in components.items():
                for variable in variables:

                    simulated_file = (
                        netcdf_dir / f"{variable}_dailyTot_output.nc"
                    )
                    if not simulated_file.exists():
                        raise FileNotFoundError(
                            f"NetCDF file {simulated_file} does not exist."
                        )
                    with xr.open_dataarray(simulated_file) as simulated:
                        if variable == "discharge":
                            ntimesteps = simulated.sizes["time"]
                            simulated = simulated.isel(time=slice(1, None))
                            simulated = simulated.mean(dim="time")
                            pit_location = simulated.argmax(dim=("lat", "lon"))
                            simulated = simulated.isel(
                                lat=pit_location["lat"],
                                lon=pit_location["lon"],
                            )
                            simulated = (
                                simulated * 86400 * (ntimesteps - 2) * 1e-9
                            )
                        elif component == "stores":
                            simulated = simulated.isel(time=[0, -1])
                            simulated = simulated * area * 1e-9
                            simulated = simulated.diff(dim="time")
                        elif component in ["influx", "outflux"]:
                            simulated = simulated.isel(time=slice(1, None))
                            simulated = simulated.sum(dim="time")
                            simulated = simulated * area * 1e-9
                        else:
                            raise ValueError(
                                f"Unknown component type: {component}"
                            )
                        simulated = simulated.sum()
                        table = {
                            "value": simulated.item(),
                            "variable": variable,
                            "component": component,
                            "system": system,
                        }
                        table = pd.DataFrame(table, index=[0])
                        tables.append(table)

            table = pd.concat(tables)
            table = table.set_index(["system", "component", "variable"])

            table_dir = commit_dir / "table"
            table_out = table_dir / "waterbalance.csv"
            table_out.parent.mkdir(parents=True, exist_ok=True)
            table.to_csv(table_out)
