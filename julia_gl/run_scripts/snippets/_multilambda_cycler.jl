#===============================================================================

    multilambda_cycler.jl - LR, May 2021.

    A function that wrapps the actual diagonalization for multiple lambda values.
    This is done completely, with I/O and all. Actually, this is like a script
    that does simply one parameter set (already checked).

===============================================================================#
include("../../src/io/io.jl")
include("../../src/io/data_storage.jl")
include("../../src/gl_hamiltonian.jl")

include("../../src/operators/gl_operators.jl")
include("../../src/operators/construct_mb_operator.jl")


# This must be included after setting the LinkType, otherwise it will not be known
# to these methods at compile time (or at all, for that matter).
include("../../src/io/data_import.jl")
include("../../src/hamiltonian_construction.jl")


function multilambda_cycler(param::Dict{Any,Any})
    """ Cycles through a list of lambda values for fixed other params.
    """
    @info "==================== starting multilambda run ====================" param=param
    local latt = LinkLattice(param["L"])

    # First, read the lookup tables (+ inverse lookup table).
    (lookup_table, ilookup_table) = read_lookup_tables(param)

    # First, try to read the Hamiltonian. If not possible (for whatever reason),
    # cosntruct it.
    hamiltonian = read_hamiltonian(param)
    n_flip_vector = nothing
    if isnothing(hamiltonian)
        # Then, construct the Hamiltonian.
        @info "Setting up Hamiltonian." nthreads = Threads.nthreads() nfock=length(lookup_table)
        local time = @elapsed hamiltonian, n_flip_vector = construct_hamiltonian(lookup_table, ilookup_table, latt)
        @info " **** Done constructing the Hamiltonian. **** " time=time nonzero_entries=length(hamiltonian.data)

        if param["store_hamiltonian"]
            store_data(
                param["hamiltonian_file"],
                conf["has_charges"] ? conf["charge_label"] : conf["ws_label"],
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
            ilookup_table;
            compute_sign=(param["gauge_particles"]=="fermions")
        )
        hilbert_ops[op_label] = h_op
        @info " *** Constructed $op_label operator." time = time
    end

    # Finally, diagonalize for all lambda values specified.
    for lambda in list_from_param("lambda", param)
        @info "---------- Starting to diagonalize the Hamiltonian. ----------" lambda=lambda
        if param["full_diag"]
            @warn "Using full diagonalization!"
        end
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

        # Compute fidelity susceptibility.
        # (resolved for energy levels, for later diagnostic - this could be changed)
        # Note: this only works reliably for non-degenerate ground-states, since
        # we use the gap as the denominator (from the perturbative expression).
        if param["compute_fidelity"]
            if !isnothing(n_flip_vector)
                @info "Computing fidelity susceptibility."
                chi_f = zeros(DType, length(ev))
                hi_psi = n_flip_vector .* est[:,1] # Application of the potential term to the ground state.
                for k = 2:length(ev)
                    chi_f[k] = dot(est[:,k], hi_psi)^2 / (ev[k] - ev[1])^2
                end

                store_data(
                    param["result_file"],
                    "susceptibility"*_lambda_tag(lambda),
                    chi_f;
                    attrs=Dict("lambda"=>lambda),
                    overwrite=param["overwrite"],
                    prefix=_nflip_tag(param)
                )
            else
                @error "Fidelity susceptibility couldn't be computed."
            end
        end
        # ---

        # Compute some other observables.
        for (op_label, h_op) in hilbert_ops
            evs = [
                expectation_value(h_op, CType.(est[:,k]))
                for k=1:length(ev)
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
end
