#===============================================================================

    param_checks.jl - LR, January 2021

    This ensures a functional config specifcation and it is the central place
    to set default values.

===============================================================================#
function param_checks!(conf::Dict{Any,Any})::Dict{Any,Any}
    """ Checks if crucial parameters are present and sets default values. Also,
        the default paths for I/O are read/constructed in a unified manner such
        that later everything can be assumed to exist.
    """
    if !haskey(conf, "L")
        error("Fatal: No lattice specified!")
    end

    if !haskey(conf, "J")
        conf["J"] = -1.0
        @warn "Set coupling J to default value." J=conf["J"]
    end

    # Lambda is required.
    if !(haskey(conf, "lambda") || haskey(conf, "lambda_range") || haskey(conf, "lambda_list"))
        error("Fatal: No coupling lambda specified!")
    end

    if !haskey(conf, "gauge_particles")
        error("Fatal: gauge_particles are required!")
    end

    if !haskey(conf, "observables")
        conf["observables"] = []
    end

    # Low energy run.
    if !haskey(conf, "low_energy_run")
        conf["low_energy_run"] = false
        @info "Running in regular mode by default."
    else
        if !haskey(conf, "maximum_excitation_level")
            error("Fatal: No excitation level specified for low-energy run.")
        end

        if !haskey(conf, "notification_level")
            conf["notification_level"] = 0
            @info "Push notifications are off."
        end
    end

    # Winding sector.
    if !haskey(conf, "winding_sector")
        conf["ws_label"] = "all-ws"
    else
        conf["ws_label"] = _winding_tag(conf["winding_sector"], latt=LinkLattice(conf["L"]))
    end
    @info "Set the winding sector label." ws_label=conf["ws_label"]


    # ------------------------------------------------
    # Some stuff for diagonalization.
    if !haskey(conf, "n_eigenvalues")
        conf["n_eigenvalues"] = 15
        @info "Setting default # of eigenvalues." n_eigenvalues=conf["n_eigenvalues"]
    end

    if !haskey(conf, "n_eigenstates")
        conf["n_eigenstates"] = 0
        @info "Setting default # of eigenstates." n_eigenstatees=conf["n_eigenstates"]
    else
        if conf["n_eigenstates"] > conf["n_eigenvalues"]
            conf["n_eigenstates"] = conf["n_eigenvalues"]
            @warn "Can't have more eigenstates than eigenvalues, using maximal #." n_eigenstates = conf["n_eigenstates"]
        end
    end

    if !haskey(conf, "ev_type")
        conf["ev_type"] = "SA"
        @info "Setting default style of eigenvalues."  ev_type = conf["ev_type"]
    end

    if !haskey(conf, "full_diag")
        conf["full_diag"] = false
        @info "Full diagonalization disabled by default."
    end

    # ------------------------------------------------
    # I/O business.

    # Overwrite.
    if !haskey(conf, "overwrite")
        conf["overwrite"] = false
        @info "Overwriting results is disabled by default."
    end

    # Logfile.
    if !haskey(conf, "logfile")
        conf["logfile"] = "logfile.log"
        @info "Using default logfilename." logfile = conf["logfile"]
    end

    # Working directory.
    # [default: ./workdir_gl/ ]
    if !haskey(conf, "working_directory")
        conf["working_directory"] = "./workdir_gl/"
        if !isdir(conf["working_directory"])
            mkdir(conf["working_directory"])
            @warn "No working directory specified, created directory at default path." working_directory = conf["working_directory"]
        else
            @info "No working directory specified, using default path." working_directory = conf["working_directory"]
        end
    end

    # State file.
    # (works for both, low-energy and full runs)
    if !haskey(conf, "state_file")
        conf["state_file"] = (
            conf["working_directory"] *
            "/" *
            (conf["low_energy_run"] ? "le_" : "winding_") *
            "states_" *
            _size_tag(conf["L"]) *
            ".hdf5"
        )
        @info "Using default file for GL states." state_file = conf["state_file"]
    end

    # Hamiltonian file.
    if !haskey(conf, "store_hamiltonian")
        conf["store_hamiltonian"] = false
        @info "Hamiltonian not stored by default."
    end
    if !haskey(conf, "hamiltonian_file")
        @info "Using default file for Hamiltonian." hamiltonian_file = conf["hamiltonian_file"] = nothing
        conf["hamiltonian_file"] = (
            conf["working_directory"] *
            "/"*(conf["low_energy_run"] ? "le_" : "")*
            "hamiltonian_"*
            _size_tag(conf["L"]) *
            ".hdf5"
        )
    end

    # Result file.
    if !haskey(conf, "result_file")
        conf["result_file"] = (
            conf["working_directory"] *
            "/"*(conf["low_energy_run"] ? "le_" : "")*
            "results_"*
            conf["gauge_particles"] * "_" *
            conf["ws_label"] * "_" *
            _size_tag(conf["L"]) *
            ".hdf5"
        )
        @info "Using default file for results." result_file = conf["result_file"]
    end

    return conf
end

# ==============================================================================
# Auxiliary methods to ensure consistent naming for I/O stuff.

function _size_tag(L)
    stag = ""
    for l in L
        stag *= "$(l)x"
    end
    return stag[begin:end-1]
end


function _winding_shift(latt::LinkLattice, ws::Array{Int,1})::Array{Int,1}
    """ Maps between the representation of winding numbers.
    """
    if latt.d == 2
        shift = reverse(latt.L) .÷ 2
    elseif latt.d == 3
        shift = [
            latt.L[2]*latt.L[3] ÷ 2;
            latt.L[1]*latt.L[3] ÷ 2;
            latt.L[1]*latt.L[2] ÷ 2
        ]
    else
        error("Dimension not implemented!")
    end
    return ws + shift
end


function _winding_tag(ws::Array{Int,1}; labels=["x", "y", "z"], latt::Union{LinkLattice,Nothing}=nothing)
    """ Returns the naming convention of the winding datasets.

        The shift is a lattice configuration L = [Lx,Ly,Lz], such
        that the HDF5 datasets may be resolved.
    """
    if !isnothing(latt)
        ws = _winding_shift(latt, ws)
    end
    wtag = ""
    for k=1:length(ws)
        wtag *= "w$(labels[k])_$(ws[k])-"
    end
    return wtag[begin:end-1]
end
