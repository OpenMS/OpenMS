#!/bin/bash
# Fix remaining enum references after migration

# Fix SpectrumSettings.h - member variable initializer
sed -i '' 's/SpectrumType type_ = UNKNOWN;/SpectrumType type_ = SpectrumType::UNKNOWN;/' \
  src/openms/include/OpenMS/METADATA/SpectrumSettings.h

# Fix ChromatogramSettings.h - array size with arithmetic
sed -i '' 's/ChromatogramNames\[SIZE_OF_CHROMATOGRAM_TYPE+1\]/ChromatogramNames[static_cast<size_t>(ChromatogramType::SIZE_OF_CHROMATOGRAM_TYPE)+1]/' \
  src/openms/include/OpenMS/METADATA/ChromatogramSettings.h

echo "Fixed enum references in headers"
