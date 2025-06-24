import pathlib as pl
import os
import shutil as sh
import git

simulation_dir = pl.Path("simulation")
model_dir = pl.Path("pcrglobwb")
conda_dir = pl.Path("conda")

model_github = (
    "https://{username}:{token}@github.com/UU-Hydro/PCR-GLOBWB_model.git"
)
commit = "2e3a5d85b0a1c264c7cb22dbc864bce87d115053"

# model_subdir = pl.Path("pcrglobwb")
# conda_subdir = pl.Path("conda")
bare_subdir = pl.Path("bare")
reference_subdir = pl.Path("reference")
parameter_subdir = pl.Path("parameters")

username = os.environ.get("GITHUB_USERNAME")
token = os.environ.get("GITHUB_TOKEN")
model_github = model_github.format(username=username, token=token)

simulations = simulation_dir.iterdir()
simulations = sorted(simulations)

print("setup pcrglobwb model")

# Setup
model_conda_file = model_dir / "conda_env" / "model_py3_standard.yml"
model_runner_file = model_dir / "model" / "deterministic_runner.py"
if not model_dir.exists():
    _ = git.Repo.clone_from(
        url=model_github,
        to_path=model_dir,
    )
repo = git.Repo(model_dir)

# Update
repo.remotes.origin.fetch(refspec="refs/heads/master")
repo.remotes.origin.pull(refspec="refs/heads/master")
if commit == "latest":
    commit = repo.git.log(n=1, pretty="format:%H")
    print(f"> commit: {commit}")

# Checkout
repo.git.checkout(commit, force=True)

print("setup conda environment")

# Setup
conda_commit_dir = conda_dir / commit
if not conda_commit_dir.exists():
    status = os.system(
        f"conda env create \
            --prefix {conda_commit_dir} \
                --file {model_conda_file} 1> /dev/null"
    )
    if status != 0:
        raise RuntimeError(f"> Failed: {status}")

simulation = simulations[0]
for simulation in simulations:
    print(f"simulation: {simulation}")

    # Setup
    configuration_file = simulation / "configuration.ini"
    simulation_parameter_dir = simulation / parameter_subdir
    simulation_reference_dir = simulation / reference_subdir
    simulation_commit_dir = simulation_reference_dir / commit
    completed_out = simulation_commit_dir / "completed.txt"
    configuration_out = simulation_commit_dir / "configuration.ini"
    time_out = simulation_commit_dir / "time.txt"
    out_out = simulation_commit_dir / "simulation.out"
    err_out = simulation_commit_dir / "simulation.err"

    # Check
    if completed_out.exists():
        print(f"> Simulation already completed")
        continue

    # Cleanup
    simulation_commit_dir.mkdir(parents=True, exist_ok=True)
    if configuration_out.exists():
        configuration_out.unlink()
    if time_out.exists():
        time_out.unlink()
    if out_out.exists():
        out_out.unlink()
    if err_out.exists():
        err_out.unlink()

    # Configuration
    sh.copy(configuration_file, configuration_out)
    with open(configuration_out, "r") as f:
        configuration = f.read()
    configuration = configuration.format(
        inputDir=simulation_parameter_dir.resolve(),
        outputDir=simulation_commit_dir.resolve(),
    )
    with open(configuration_out, "w") as f:
        f.write(configuration)

    # Simulation
    command = (
        f"conda run --prefix {conda_commit_dir} "
        f"python {model_runner_file} {configuration_out} "
        f"1> {out_out} 2> {err_out}"
    )
    command = f"{{ time {command} ; }} 2> {time_out}"
    status = os.system(command)
    if status != 0:
        raise RuntimeError(f"> Failed: {status}")

    # Completion
    completed_out.touch()
