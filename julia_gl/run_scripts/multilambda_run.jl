#===============================================================================

    single_diagonalization.jl - LR, January 2021

    Serves to diagonalize a single value. First test with Julia. In essence,
    this is to be replaced with a multilambda run later.

===============================================================================#
ENV["JULIA_DEBUG"]="all"
@debug "Debugging activated!"

include("../src/typedefs.jl")
include("../src/io/io.jl")
include("../src/io/logging.jl")
include("../src/param_checks.jl")
include("../src/io/data_storage.jl")
include("../src/gl_hamiltonian.jl")

include("../src/operators/gl_operators.jl")
include("../src/operators/construct_mb_operator.jl")

# Read the config file (first argument after the program name).
# param = param_checks!(read_config(ARGS[1]))
param = param_checks!(read_config("config.yml"))
make_logger(param)



@info "==================== starting multilambda run ====================" param=param

# Create the lattice structure that we're working on.
latt = LinkLattice(param["L"])

# Set the size for the Link Type.
# (when dealing with links, this type should be used to convert the datatype to
# be consistent - otherwise there'll likely be an error or *very* subtle bugs)
const LinkType = latt.S[end]*latt.d >= 63 ? LargeLinkState : SmallLinkState
@info "Set the datatype for link representation" type=LinkType

# This must be included after setting the LinkType, otherwise it will not be known
# to these methods at compile time (or at all, for that matter).
include("../src/io/data_import.jl")
include("../src/hamiltonian_construction.jl")


# First, read the lookup tables (+ inverse lookup table).
(lookup_table, ilookup_table) = read_lookup_tables(param)


# First, try to read the Hamiltonian. If not possible (for whatever reason),
# cosntruct it.
hamiltonian = read_hamiltonian(param)
if isnothing(hamiltonian)
    # Then, construct the Hamiltonian.
    @info "Setting up Hamiltonian." nthreads = Threads.nthreads() nfock=length(lookup_table)
    local time = @elapsed hamiltonian = construct_hamiltonian(lookup_table, ilookup_table, latt)
    @info " **** Done constructing the Hamiltonian. **** " time=time nonzero_entries=length(hamiltonian.data)

    if param["store_hamiltonian"]
        store_data(
            param["hamiltonian_file"],
            param["ws_label"],
            hcat(hamiltonian.col, hamiltonian.row, hamiltonian.data);
            overwrite=param["overwrite"],
            prefix=_nflip_tag(param),
            attrs=Dict("n_fock"=>length(lookup_table))
        )
        @info "Stored Hamiltonian." file=param["hamiltonian_file"]
    end
end

# Make all operators that we need.
@info "Setting up additional Hilbert operators." operators=param["observables"]
hilbert_ops = Dict{String,HilbertOperator}()
for op_label in param["observables"]
    local time = @elapsed h_op = construct_operator(
        gl_operators[join(latt.L, "x")][op_label],
        lookup_table,
        ilookup_table
    )
    hilbert_ops[op_label] = h_op
    @info " *** Constructed $op_label operator." time = time
end

# Finally, diagonalize for all lambda values specified.
for lambda in list_from_param("lambda", param)
    @info "---------- Starting to diagonalize the Hamiltonian. ----------" lambda=lambda
    local time = @elapsed local (ev, est) = diagonalize(hamiltonian, lambda, param)
    @info "Done computing the lower spectrum." time=time spectrum=ev

    # ---
    # Output the data.

    # Eigenvalues.
    store_data(
        param["result_file"],
        "spectrum"*_lambda_tag(lambda),
        ev;
        attrs=Dict("lambda"=>lambda),
        overwrite=param["overwrite"],
        prefix=_nflip_tag(param)
    )

    # Eigenstates.
    store_data(
        param["result_file"],
        "eigenstates"*_lambda_tag(lambda),
        est[:,begin:param["n_eigenstates"]];
        attrs=Dict("lambda"=>lambda),
        overwrite=param["overwrite"],
        prefix=_nflip_tag(param)
    )

    # ---
    # Compute some other observables.
    for (op_label, h_op) in hilbert_ops
        evs = [
            expectation_value(h_op, CType.(est[:,k]))
            for k=1:param["n_eigenvalues"]
        ]
        @info "Computed $op_label eigenvalues." evs=evs

        # Export.
        store_data(
            param["result_file"],
            "$op_label"*_lambda_tag(lambda),
            evs;
            attrs=Dict("lambda"=>lambda),
            overwrite=param["overwrite"],
            prefix=_nflip_tag(param)
        )
    end
end

@info "============================== fin =============================="
