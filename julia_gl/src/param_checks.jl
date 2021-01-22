#===============================================================================

    param_checks.jl - LR, January 2021

    This ensures a functional config specifcation and it is the central place
    to set default values.

===============================================================================#
function param_checks!(conf::Dict{Any,Any})::Dict{Any,Any}
    """ Checks if crucial parameters are present and sets default values.
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

    # Working directory.
    if !haskey(conf, "working_directory")
        conf["working_directory"] = "./workdir_gl/"
        if !isdir(conf["working_directory"])
            mkdir(conf["working_directory"])
            @warn "No working directory specified, created directory at default path." working_directory = conf["working_directory"]
        else
            @info "No working directory specified, using default path." working_directory = conf["working_directory"]
        end
    end

    # Low energy run.
    if !haskey(conf, "low_enegy_run")
        conf["low_enegy_run"] = false
    else
        if !haskey(conf, "maximum_excitation_level")
            error("Fatal: No excitation level specified for low-energy run.")
        end

        if !haskey(conf, "notification_level")
            conf["notification_level"] = 0
            @info "Push notifications are off."
        end
    end

    return conf
end
