# Function to add files to a source group, stripping a base path if necessary
function(add_to_source_group base_path group_prefix files)
    set(sources_with_path "")
    foreach(file ${files})
        list(APPEND sources_with_path "${base_path}/${file}")
    endforeach()

    # Create a formatted group name by removing the base path and replacing slashes
    string(REPLACE "${CMAKE_SOURCE_DIR}/" "" group_name "${base_path}")
    string(REPLACE "/" "\\" group_name "${group_name}")

    # Prefix the group name if provided
    if(NOT "${group_prefix}" STREQUAL "")
        set(group_name "${group_prefix}\\${group_name}")
    endif()

    source_group("${group_name}" FILES ${sources_with_path})
endfunction()
