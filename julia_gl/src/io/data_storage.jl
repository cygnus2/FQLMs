using HDF5, Dates

function store_data(filename, ds_name, data; attrs=Dict(), prefix="", overwrite=false)
    """ Stores the results in standardized fashion. This should be the only
        place where any output is generated - it allows for a coherent
        tagging of the results/output produced.
    """
    all_attrs = attrs
    # all_attrs["time"] = Dates.format(now(), date_format)[2:end-1]
    # all_attrs["host"] = gethostname()
    # all_attrs['version'] = version

    ds = prefix*"/"*ds_name
    _write_dataset(filename, ds, data; overwrite=overwrite)
    h5writeattr(filename, ds, all_attrs)
end


function _write_dataset(filename::String, path::String, data; overwrite=false)
    """ Definitely writes a dataset to file, regardless whether it exists or not.
    """
    repack = false
    h5open(filename, "cw") do file
        if overwrite
            repack = _unlink_ds_if_exists(file, path)
        end
        write(file, path, data)
    end
    if repack
        @warn "HDF5 file should be repacked!"
    end
end

function _unlink_ds_if_exists(file::HDF5.File, path::String)::Bool
    """ Apparently HDF5 does not want us to delete datasets. This is a bit of a
        dirty hack, which essentially unlinks the dataset so that we can reuse
        the address (loosely speaking). This does not free up disk space, though,
        so we would need to repack the file to do this.

        Theres next to little information on this. The routine h5l_delete() is
        discussed in the HDF5 documentation. However, I have no idea what the
        third parameter does - I set it to 0 and it seemed to work this way.
        (DEFINITELY nees a more profound check at some point)
    """
    if haskey(file, path)
        HDF5.h5l_delete(file, path, 0)
        return true
    end
    return false
end
