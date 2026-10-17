cmake_minimum_required(VERSION 3.24)

if(NOT BINARY_DIR)
  message(FATAL_ERROR "Pass -DBINARY_DIR=<scratch directory>.")
endif()

function(check_s3 name enabled s3 module expected diagnostic)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${CMAKE_CURRENT_LIST_DIR}" -B "${BINARY_DIR}/${name}"
      "-DWITH_S3=${enabled}" "-DARROW_S3=${s3}" "-DARROW_S3_MODULE=${module}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(expected STREQUAL "success")
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "${name} unexpectedly failed:\n${output}\n${error}")
    endif()
  else()
    if(result EQUAL 0 OR NOT "${output}\n${error}" MATCHES "${diagnostic}")
      message(FATAL_ERROR "${name} did not report the expected failure:\n${output}\n${error}")
    endif()
    if(NOT "${output}\n${error}" MATCHES "WITH_S3=OFF")
      message(FATAL_ERROR "${name} did not explain how to opt out:\n${error}")
    endif()
  endif()
endfunction()

check_s3(builtin ON ON OFF success "")
check_s3(no_s3 ON OFF OFF failure "ARROW_S3=ON")
check_s3(module ON ON ON failure "Separate dynamic Arrow S3 modules")
check_s3(disabled OFF OFF OFF success "")
check_s3(disabled_module OFF ON ON success "")
check_s3(unknown ON "" "" failure "ARROW_S3=ON")

message(STATUS "Arrow S3 configuration checks passed")
