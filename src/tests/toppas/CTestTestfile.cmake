# CMake generated Testfile for 
# Source directory: /workspaces/OpenMS/src/tests/toppas
# Build directory: /workspaces/OpenMS/src/tests/toppas
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TOPP_ExecutePipeline_1 "/workspaces/OpenMS/bin/ExecutePipeline" "-test" "-in" "/workspaces/OpenMS/src/tests/toppas/ExecutePipeline_1.toppas" "-resource_file" "/workspaces/OpenMS/src/tests/toppas/ExecutePipeline_1.trf" "-out_dir" ".")
set_tests_properties(TOPP_ExecutePipeline_1 PROPERTIES  _BACKTRACE_TRIPLES "/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;65;add_test;/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;0;")
add_test(TOPP_ExecutePipeline_BSA_Quantitation.toppas "/workspaces/OpenMS/bin/ExecutePipeline" "-test" "-in" "/workspaces/OpenMS/share/OpenMS/examples/TOPPAS/BSA_Quantitation.toppas" "-out_dir" "." "-num_jobs" "4")
set_tests_properties(TOPP_ExecutePipeline_BSA_Quantitation.toppas PROPERTIES  SKIP_RETURN_CODE "14" _BACKTRACE_TRIPLES "/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;78;add_test;/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;0;")
add_test(TOPP_ExecutePipeline_Ecoli_Identification.toppas "/workspaces/OpenMS/bin/ExecutePipeline" "-test" "-in" "/workspaces/OpenMS/share/OpenMS/examples/TOPPAS/Ecoli_Identification.toppas" "-out_dir" "." "-num_jobs" "4")
set_tests_properties(TOPP_ExecutePipeline_Ecoli_Identification.toppas PROPERTIES  SKIP_RETURN_CODE "14" _BACKTRACE_TRIPLES "/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;78;add_test;/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;0;")
add_test(TOPP_ExecutePipeline_FDR_NeighborSearch.toppas "/workspaces/OpenMS/bin/ExecutePipeline" "-test" "-in" "/workspaces/OpenMS/share/OpenMS/examples/TOPPAS/FDR_NeighborSearch.toppas" "-out_dir" "." "-num_jobs" "4")
set_tests_properties(TOPP_ExecutePipeline_FDR_NeighborSearch.toppas PROPERTIES  SKIP_RETURN_CODE "14" _BACKTRACE_TRIPLES "/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;78;add_test;/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;0;")
add_test(TOPP_ExecutePipeline_QualityControl.toppas "/workspaces/OpenMS/bin/ExecutePipeline" "-test" "-in" "/workspaces/OpenMS/share/OpenMS/examples/TOPPAS/QualityControl.toppas" "-out_dir" "." "-num_jobs" "4")
set_tests_properties(TOPP_ExecutePipeline_QualityControl.toppas PROPERTIES  SKIP_RETURN_CODE "14" _BACKTRACE_TRIPLES "/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;78;add_test;/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;0;")
add_test(TOPP_ExecutePipeline_merger_tutorial.toppas "/workspaces/OpenMS/bin/ExecutePipeline" "-test" "-in" "/workspaces/OpenMS/share/OpenMS/examples/TOPPAS/merger_tutorial.toppas" "-out_dir" "." "-num_jobs" "4")
set_tests_properties(TOPP_ExecutePipeline_merger_tutorial.toppas PROPERTIES  SKIP_RETURN_CODE "14" _BACKTRACE_TRIPLES "/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;78;add_test;/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;0;")
add_test(TOPP_ExecutePipeline_peakpicker_tutorial.toppas "/workspaces/OpenMS/bin/ExecutePipeline" "-test" "-in" "/workspaces/OpenMS/share/OpenMS/examples/TOPPAS/peakpicker_tutorial.toppas" "-out_dir" "." "-num_jobs" "4")
set_tests_properties(TOPP_ExecutePipeline_peakpicker_tutorial.toppas PROPERTIES  SKIP_RETURN_CODE "14" _BACKTRACE_TRIPLES "/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;78;add_test;/workspaces/OpenMS/src/tests/toppas/CMakeLists.txt;0;")
