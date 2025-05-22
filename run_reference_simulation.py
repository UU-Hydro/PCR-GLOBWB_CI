import pathlib as pl
import os
import shutil as sh
import git

simulation_dir = pl.Path('simulation')
pcrglobwb_github = 'https://{username}:{token}@github.com/UU-Hydro/PCR-GLOBWB_model.git'
pcrglobwb_github = 'https://{username}:{token}@github.com/UU-Hydro/PCR-GLOBWB_3.git'
commit = 'latest'

pcrglobwb_subdir = pl.Path('pcrglobwb')
conda_subdir = pl.Path('conda')
reference_subdir = pl.Path('reference')
parameter_subdir = pl.Path('parameters')

username = os.environ.get('GITHUB_USERNAME')
token = os.environ.get('GITHUB_TOKEN')
pcrglobwb_github = pcrglobwb_github.format(
    username=username,
    token=token
)

simulations = simulation_dir.iterdir()
simulations = sorted(simulations)

simulation = simulations[0]
for simulation in simulations:
    if simulation.name != 'rhine_tinyrun':
        continue
    print(f'simulation: {simulation}')

    print('- setup pcrglobwb model')
    pcrglobwb_dir = simulation / pcrglobwb_subdir
    if pcrglobwb_dir.exists():
        sh.rmtree(pcrglobwb_dir)
    pcrglobwb_repo = git.Repo.clone_from(url=pcrglobwb_github,
                                         to_path=pcrglobwb_dir)
    if commit == 'latest':
        commit = pcrglobwb_repo.git.log(n=1, pretty='format:%H')
    pcrglobwb_repo.git.checkout(commit, force=True)
    print(f'commit: {commit}')

    print('- setup commit reference')
    reference_dir = simulation / reference_subdir
    reference_dir = reference_dir / commit
    completed_file = reference_dir / 'completed.txt'
    if completed_file.exists():
        sh.rmtree(pcrglobwb_dir)
        print(f'> Simulation already completed')
        continue
    reference_dir.mkdir(parents=True, exist_ok=True)

    print('- setup conda environment')
    conda_dir = simulation / conda_subdir
    if conda_dir.exists():
        sh.rmtree(conda_dir)
    environment_file = pcrglobwb_dir / 'conda_env' / 'pcrglobwb_py3_standard.yml'
    if not environment_file.exists():
        raise FileNotFoundError(
            f'Environment file not found: {environment_file}')
    return_value = os.system(
        f'conda env create --prefix {conda_dir} --file {environment_file} 1> /dev/null')
    if return_value != 0:
        raise RuntimeError(
            f'Failed to create conda environment: {conda_dir}')

    print('- setup simulation configuration')
    configuration_file = simulation / 'configuration.ini'
    if not configuration_file.exists():
        raise FileNotFoundError(
            f'Configuration file not found: {configuration_file}')
    configuration_out = reference_dir / 'configuration.ini'
    sh.copy(configuration_file, configuration_out)
    parameter_dir = simulation / parameter_subdir
    if not parameter_dir.exists():
        raise FileNotFoundError(
            f'Parameter directory not found: {parameter_dir}')
    with open(configuration_out, 'r') as f:
        configuration = f.read()
    configuration = configuration.format(inputDir=parameter_dir.resolve(),
                                         outputDir=reference_dir.resolve())
    with open(configuration_out, 'w') as f:
        f.write(configuration)

    print('- run simulation')
    runner_file = pcrglobwb_dir / 'model' / 'deterministic_runner.py'
    if not runner_file.exists():
        raise FileNotFoundError(
            f'Runner file not found: {runner_file}')
    time_file = reference_dir / 'time.txt'
    if time_file.exists():
        time_file.unlink()
    out_file = reference_dir / 'simulation.out'
    if out_file.exists():
        out_file.unlink()
    err_file = reference_dir / 'simulation.err'
    if err_file.exists():
        err_file.unlink()

    return_value = os.system(
        f'{{ time conda run --prefix {conda_dir} python {runner_file} {configuration_out} 1> {out_file} 2> {err_file} ; }} 2> {time_file}')
    if return_value != 0:
        raise RuntimeError(
            f'Failed to run simulation: {runner_file}')

    print('- cleanup')
    sh.rmtree(pcrglobwb_dir)
    sh.rmtree(conda_dir)
    completed_file.touch()
