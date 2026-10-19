cmake_minimum_required(VERSION 3.24)
project(OpenMSArrowS3Configuration NONE)

include("${CMAKE_CURRENT_LIST_DIR}/../../../cmake/OpenMSArrowS3.cmake")
openms_check_arrow_s3()
