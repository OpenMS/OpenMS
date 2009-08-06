### the directory name
set(directory include/OpenMS/COMPARISON/SPECTRA)

### list all header files of the directory here
set(sources_list_h
BinnedSharedPeakCount.h
BinnedSpectralContrastAngle.h
BinnedSpectrum.h
BinnedSpectrumCompareFunctor.h
BinnedSumAgreeingIntensities.h
CompareFouriertransform.h
PeakAlignment.h
PeakSpectrumCompareFunctor.h
SpectraSTSimilarityScore.h
SpectrumAlignment.h
SpectrumAlignmentScore.h
SpectrumCheapDPCorr.h
SpectrumPrecursorComparator.h
SteinScottImproveScore.h
ZhangSimilarityScore.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\COMPARISON\\SPECTRA" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

