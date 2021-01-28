FQLM - Julia style
===================
This is the same as the Python version but quite a bit faster.

## Setup
The setup is relatively easy: Make sure you install the current version of Julia (the code was written and tested with `Julia 1.5.3`). This is simply achieved by downloading Julia and append to the path the location to the binary.

Then, to get all dependencies for this code, set up a local environment by entering the Julia REPL and type `]` which leads you to the  `Pkg` REPL. Then type `activate <path/to/root_folder>` where `<path/to/root_folder>` should be replaced by the actual path to the directory `julia_gl` (this directory). If you started the REPL from this folder simply put `.` there. Finally, in order to download and "install" all packages for this project issue `instantiate`. Once this is done, you're ready to go.


## Usage
The usage is exactly the same as for the Python case. Simply run `julia --project=<path-to-root_folder> multilambda-run.jl <config.yml>`. The specification of the project is important, otherwise you're not running with the required dependencies.

### Best practice
Several files are generated in the process of running this code, which are also needed for later stages. Although possible, it's best to not mess with the filenames there and just work with the default values. Simply specify a `working_directory` in the config file, evyerhting will be stored/read from there.
If larger runs are done on the cluster, it might be beneficial to store states and hamiltonians of CPU heavy runs. Those can specifically be targeted by supplying `state_file` and `hamiltonian_file` in the config. If those are not specified, they must be present in the working directory, otherwise the code won't work (or they are generated on the fly, which could take a bit).

## Benchmarks
I ran some quick tests during the development. On a single core, the construction of the `2x2x4` Hamiltonian took about 35 seconds as opposed to 40 *minutes* with the Python implementation. I don't know whether this benchmark will survive for other use cases, but it seems to be quite helpful.
Moreover, the usage of memory is greatly reduced, since it's much more straightforward to specify smaller types whenever those are needed.


## TODO
Some things are left to do here:
- Implement multithreading of Hamiltonian construction (not at all the bottleneck at the moment though)
- State search (currently only Python-found states are read in)
- Better structure for Hilbert Operators (more general?)
