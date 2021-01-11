//! [Test]

# DatabaseFilter test:
add_test("UTILS_DatabaseFilter_1" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_1.fasta -accession ${DATA_DIR_TOPP}/DatabaseFilter_1.idXML -out DatabaseFilter_1_out.fasta.tmp)
add_test("UTILS_DatabaseFilter_1_out" ${DIFF} -in1 DatabaseFilter_1_out.fasta.tmp -in2 ${DATA_DIR_TOPP}/DatabaseFilter_1_out.fasta )
set_tests_properties("UTILS_DatabaseFilter_1_out" PROPERTIES DEPENDS "UTILS_DatabaseFilter_1")
add_test("UTILS_DatabaseFilter_2" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_1.fasta -accession ${DATA_DIR_TOPP}/DatabaseFilter_1.idXML -out DatabaseFilter_2_out.fasta.tmp -method blacklist)
add_test("UTILS_DatabaseFilter_2_out" ${DIFF} -in1 DatabaseFilter_2_out.fasta.tmp -in2 ${DATA_DIR_TOPP}/DatabaseFilter_2_out.fasta )
set_tests_properties("UTILS_DatabaseFilter_2_out" PROPERTIES DEPENDS "UTILS_DatabaseFilter_2")

//! [Test]
