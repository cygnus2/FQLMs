#===============================================================================

    param_checks.jl - LR, January 2021

    This ensures a functional config specifcation and it is the central place
    to set default values.

===============================================================================#
function param_checks!(conf::Dict{Any,Any}; state_run::Bool=false)::Dict{Any,Any}
    """ Checks if crucial parameters are present and sets default values. Also,
        the default paths for I/O are read/constructed in a unified manner such
        that later everything can be assumed to exist.
    """
    if !haskey(conf, "L")
        error("Fatal: No lattice specified!")
    end

    # Set a few default values.
    defaults = Dict(
        "J" => -1.0,
        "winding_sector" => "all-ws",

        "low_energy_run" => false,
        "observables" => [],

        "n_eigenvalues" => 15,
        "n_eigenstates" => 0,
        "ev_type" => "SA",
        "full_diag" => false,
        "compute_fidelity" => false,

        "logfile" => "logfile.log",
        "overwrite" => false,
        "working_directory" => "./workdir_gl/",
        "store_hamiltonian" => false,
        "read_hamiltonian" => true
    )
    for (k,v) in defaults
        if !haskey(conf, k)
            @info "Setting default value for '$k'." value=v
            conf[k] = v
        end
    end


    # ----------------------------------------------------
    # Check if essential parameters are given.
    required = nothing
    if state_run
        required = [
            ["L"],
            ["base_states"],
            ["winding_sector"]
        ]
    else
        required = [
            ["L"],
            ["gauge_particles"],
            ["lambda", "lambda_range", "lambda_list"]
        ]
    end
    for keys in required
        if !any([haskey(conf, key) for key in keys])
            error("Fatal: required parameter $keys not specified!")
        end
    end


    # ----------------------------------------------------
    # Low energy run.
    if conf["low_energy_run"]
        if !haskey(conf, "maximum_excitation_level")
            error("Fatal: No excitation level specified for low-energy run.")
        end

        if !haskey(conf, "notification_level")
            conf["notification_level"] = 0
            @info "Push notifications are off."
        end
    end


    # ------------------------------------------------
    # Winding sector.
    conf["ws_label"] = conf["winding_sector"] == "all-ws" ? "all-ws" : _winding_tag(conf["winding_sector"], latt=LinkLattice(conf["L"]))
    @info "Set the winding sector label." ws_label=conf["ws_label"]
    # ------------------------------------------------
    # Some stuff for diagonalization.
    if conf["n_eigenstates"] > conf["n_eigenvalues"]
        conf["n_eigenstates"] = conf["n_eigenvalues"]
        @warn "Can't have more eigenstates than eigenvalues, using maximal #." n_eigenstates = conf["n_eigenstates"]
    end

    # ------------------------------------------------
    # I/O business.

    # Working directory.
    if !isdir(conf["working_directory"])
        mkdir(conf["working_directory"])
        @warn "Created working directory at path." working_directory = conf["working_directory"]
    else
        @info "Using working directory at path." working_directory = conf["working_directory"]
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


    if !state_run

        # Hamiltonian file.
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
