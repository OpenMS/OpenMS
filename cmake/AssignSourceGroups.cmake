# Function to automatically create source groups for files based on their directory structure
function(assign_source_groups)
    foreach(FILE_PATH ${ARGN})
        # Get the directory part of the file path
        get_filename_component(PARENT_DIR "${FILE_PATH}" DIRECTORY)
        # Replace the base source or include directory path to form the group structure
        string(REPLACE "${CMAKE_SOURCE_DIR}/source/OPENSWATHALGO/" "" GROUP "${PARENT_DIR}")
        string(REPLACE "${CMAKE_SOURCE_DIR}/include/OpenMS/OPENSWATHALGO/" "" GROUP "${GROUP}")
        string(REPLACE "/" "\\" GROUP "${GROUP}")

        # Assign the file to the appropriate source group
        source_group("${GROUP}" FILES "${FILE_PATH}")
    endforeach()
endfunction()
