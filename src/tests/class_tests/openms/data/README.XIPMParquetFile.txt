XIPMParquetFile_reader_input.xipm and XIPMParquetFile_reader_empty.xipm
are independent reader fixtures. The populated file contains the same two
peak maps as the former XIPMParquetFile_test producer: run 7 (run1.mzML),
one y7^1 transition at m/z 500.2 and one precursor at m/z 600.2, both with
TARGET_RT 100. Arrays use uncompressed little-endian doubles (compression 0).
The empty file has the same schema and no rows.

XIPMParquetRoundTrip_test in src/topp/OpenSwathPeakMapExtractor/test retains
the original consumer-generated fixtures and reader assertions, including
the consumer's compressed array output. This keeps reader coverage enabled
without BUILD_TOPP_TOOLS while testing the producer/reader round trip with it.
