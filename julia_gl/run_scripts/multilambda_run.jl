#===============================================================================

    single_diagonalization.jl - LR, January 2021

    Serves to diagonalize a single value. First test with Julia. In essence,
    this is to be replaced with a multilambda run later.

===============================================================================#
using Printf
include("../src/typedefs.jl")
include("../src/io/io.jl")
include("../src/io/logging.jl")
include("../src/param_checks.jl")
include("../src/io/data_storage.jl")
include("../src/io/python_import.jl")
include("../src/hamiltonian_construction.jl")
include("../src/gl_hamiltonian.jl")

# Read the config file (first argument after the program name).
# param = read_config(ARGS[1])
param = param_checks!(read_config("config.yml"))
make_logger(param)


@info "==================== starting multilambda run ====================" param=param

# Create the lattice structure that we're working on.
latt = LinkLattice(param["L"])

# Set the precision for the Link Type.
# LinkType = 

# First, read the lookup tables (+ inverse lookup table).
(ws, lookup_table, ilookup_table) = read_lookup_tables(param)

# Make file for results.
results_file = _get_result_file(param, ws)

# Then, construct the Hamiltonian.
@info "Setting up Hamiltonian." nthreads = Threads.nthreads() nfock=length(lookup_table)
time = @elapsed hamiltonian = construct_hamiltonian(lookup_table, ilookup_table, latt)

@info " **** Done constructing the Hamiltonian. **** " time=time

# Finally, diagonalize for all lambda values specified.
for lambda in list_from_param("lambda", param)
    @info "---------- Starting to diagonalize the Hamiltonian. ----------" lambda=lambda
    local time = @elapsed local (ev, est) = diagonalize(
        hamiltonian,
        param["n_eigenvalues"],
        param["ev_type"],
        param["J"],
        lambda,
        param["gauge_particles"]
    )
    @info "Done computing the lower spectrum." time=time spectrum=ev

    # Output the data.
    store_data(
        results_file,
        "spectrum_lam_"*(@sprintf "%.6f" lambda),
        ev;
        attrs=Dict("lambda"=>lambda),
        overwrite=get(param, "overwrite", false)
    )

    #TODO: store eigenvalues.
    #TODO: low energy support.
end

@info "============================== fin =============================="
