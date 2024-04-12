### the directory name
set(directory source/COMPARISON/SPECTRA)

### list all filenames of the directory here
set(sources_list
BinnedSharedPeakCount.cpp
BinnedSpectralContrastAngle.cpp
BinnedSpectrum.cpp
BinnedSpectrumCompareFunctor.cpp
BinnedSumAgreeingIntensities.cpp
PeakAlignment.cpp
PeakSpectrumCompareFunctor.cpp
SpectraSTSimilarityScore.cpp
SpectrumAlignment.cpp
SpectrumAlignmentScore.cpp
SpectrumCheapDPCorr.cpp
SpectrumPrecursorComparator.cpp
SteinScottImproveScore.cpp
ZhangSimilarityScore.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\COMPARISON\\SPECTRA" FILES ${sources})

