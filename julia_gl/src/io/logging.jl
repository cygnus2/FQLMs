#===============================================================================

    logging.jl - LR, January 2021

    Holds stuff for setting up a convenient logger.

===============================================================================#
using Logging, LoggingExtras, Dates


const date_format = "dd-mm-yyyy HH:MM:SS"

# timestamp_logger(logger) = TransformerLogger(logger) do log
#   merge(log, (; message = "$(Dates.format(now(), date_format)) $(log.message)"))
# end
#
# ConsoleLogger(stdout, Logging.Debug) |> timestamp_logger |> global_logger


function make_logger(param)
    """ Set up a logger that persists to file and writes to the console.
    """
    filename = haskey(param, "logfile") ? param["logfile"] : "logfile.log"
    logger = TeeLogger(
        ConsoleLogger(),
        SimpleLogger(open(filename, "w")),
    )
    global_logger(logger)
    global_logger(logger)
end
