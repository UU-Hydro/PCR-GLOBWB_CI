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

reference_subdir = pl.Path("reference")
parameter_subdir = pl.Path("parameters")

username = os.environ.get("GITHUB_USERNAME")
token = os.environ.get("GITHUB_TOKEN")
model_github = model_github.format(username=username, token=token)

simulation = simulation_dir / "invalid_configuration"

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

# Setup
simulation_parameter_dir = simulation / parameter_subdir
simulation_reference_dir = simulation / reference_subdir
simulation_commit_dir = simulation_reference_dir / commit
configuration_out = simulation_commit_dir / "configuration.ini"
out_out = simulation_commit_dir / "simulation.out"
err_out = simulation_commit_dir / "simulation.err"

# Cleanup
simulation_commit_dir.mkdir(parents=True, exist_ok=True)

configuration_files = simulation.glob("configuration_*.ini")
configuration_files = sorted(configuration_files)

for configuration_file in configuration_files:
    print(f"configuration: {configuration_file}")

    # Cleanup
    if configuration_out.exists():
        configuration_out.unlink()

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
    status = os.system(command)
    if status == 0:
        raise RuntimeError(f"> Completed: {status}")

# Cleanup
sh.rmtree(simulation_reference_dir)
