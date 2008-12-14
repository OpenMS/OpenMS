###

set(TOPP_BIN_PATH ${PROJECT_SOURCE_DIR}/bin)

### FileConverter tests
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff)
set(TOPPTOOL ${TOPP_BIN_PATH}/FileConverter)
add_test("FileConverter_test1" ${TOPPTOOL} -in FileConverter_1_input.txt -in_type mzData -out FileConverter.tmp -out_type DTA2D)
add_test("FileConverter_test1_output" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_1_output.dta2d)
add_test("FileConverter_test2" ${TOPPTOOL} -in FileConverter_2_input.txt -in_type DTA2D -out FileConverter.tmp -out_type DTA2D)
add_test("FileConverter_test2_output" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_2_output.dta2d)
add_test("FileConverter_test3" ${TOPPTOOL} -in FileConverter_3_input.txt -in_type featureXML -out FileConverter.tmp -out_type DTA2D)
add_test("FileConverter_test3_output" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_3_output.dta2d)
add_test("FileConverter_test4" ${TOPPTOOL} -in FileConverter_4_input.txt -in_type mzXML -out FileConverter.tmp -out_type DTA2D)
add_test("FileConverter_test4_output" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_4_output.dta2d)

if (USE_ANDIMS)
  add_test("FileConverter_test5" ${TOPPTOOL} -in FileConverter_5_input.txt -in_type cdf -out FileConverter.tmp -out_type DTA2D)
	add_test("FileConverter_test5_output" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_5_output.dta2d)
endif()

add_test("FileConverter_test6a" ${TOPPTOOL} -in FileConverter_6_input.dta2d -out FileConverter.tmp -out_type mzXML)
add_test("FileConverter_test6a_output" ${DIFF} -in1 FileConverter2.tmp -in2 FileConverter_6_input.dta2d)

add_test("FileConverter_test6b" ${TOPPTOOL} -in FileConverter_6_input.dta2d -out FileConverter.tmp -out_type mzData)
add_test("FileConverter_test6b_output" ${DIFF} -in1 FileConverter2.tmp -in2 FileConverter_6_input.dta2d)

add_test("FileConverter_test7" ${TOPPTOOL} -in FileConverter_7_input.consensusXML -out FileConverter_7_output.tmp -out_type featureXML)
add_test("FileConverter_test7_output" ${DIFF} -in1 FileConverter_7_output.tmp -in2 FileConverter_7_output.featureXML)

add_test("FileConverter_test8" ${TOPPTOOL} -in FileConverter_8_input.mzML -out FileConverter_8_output.tmp -out_type mzData)
add_test("FileConverter_test8_output" ${DIFF} -in1 FileConverter_8_output.tmp -in2 FileConverter_8_output.mzData)

