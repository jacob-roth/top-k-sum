using CxxWrap
function mod_cmakelists!()
    # Read all lines from the file
    lines = readlines("CMakeLists.txt")
    
    # Modify line 14
    julia_cxxwrapprefix = CxxWrap.prefix_path()
    lines[14] = "list(APPEND CMAKE_PREFIX_PATH \"$(julia_cxxwrapprefix)\")"

    # Write the modified lines back to the file
    open("CMakeLists.txt", "w") do file
        write(file, join(lines, "\n"))
    end
end
mod_cmakelists!()