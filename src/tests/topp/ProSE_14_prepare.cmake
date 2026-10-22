# Inputs of TOPP_ProSE_14: two copies of MZML in OUT_DIR, ProSE_14_a.mzML and the gzipped ProSE_14_b.mzML.gz, two
# copies of the same name in different folders, ProSE_14_x/same.mzML and ProSE_14_y/same.mzML, and copies of the reserved
# basenames, ProSE_14_reserved/quantms.mzML and ProSE_14_reserved/QuantMS.mzML (-out_qpx), and the gzipped
# ProSE_14_reserved/OpenMS.mzML.gz (-out_parquet).
file(COPY_FILE "${MZML}" "${OUT_DIR}/ProSE_14_a.mzML")
file(COPY_FILE "${MZML}" "${OUT_DIR}/ProSE_14_b.mzML")
file(REMOVE "${OUT_DIR}/ProSE_14_b.mzML.gz")
file(ARCHIVE_CREATE OUTPUT "${OUT_DIR}/ProSE_14_b.mzML.gz" PATHS "${OUT_DIR}/ProSE_14_b.mzML" FORMAT raw COMPRESSION GZip)
file(REMOVE "${OUT_DIR}/ProSE_14_b.mzML")
foreach(_dir ProSE_14_x ProSE_14_y)
  file(MAKE_DIRECTORY "${OUT_DIR}/${_dir}")
  file(COPY_FILE "${MZML}" "${OUT_DIR}/${_dir}/same.mzML")
endforeach()
file(MAKE_DIRECTORY "${OUT_DIR}/ProSE_14_reserved")
foreach(_name quantms QuantMS)
  file(COPY_FILE "${MZML}" "${OUT_DIR}/ProSE_14_reserved/${_name}.mzML")
endforeach()
file(COPY_FILE "${OUT_DIR}/ProSE_14_b.mzML.gz" "${OUT_DIR}/ProSE_14_reserved/OpenMS.mzML.gz")
