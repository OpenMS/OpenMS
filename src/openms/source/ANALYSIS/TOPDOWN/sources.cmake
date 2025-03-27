### the directory name
set(directory source/ANALYSIS/TOPDOWN)

### list all filenames of the directory here
set(sources_list
        DeconvolvedSpectrum.cpp
        SpectralDeconvolution.cpp
        FLASHDeconvAlgorithm.cpp
        FLASHHelperClasses.cpp
        FLASHExtenderAlgorithm.cpp
        FLASHIda.cpp
        FLASHIdaBridgeFunctions.cpp
        FLASHTnTAlgorithm.cpp
        MassFeatureTrace.cpp
        PeakGroup.cpp
        Qscore.cpp
        Qvalue.cpp
        TopDownIsobaricQuantification.cpp
        FLASHTaggerAlgorithm.cpp
)

### add path to the filenames
set(sources)
foreach (i ${sources_list})
    list(APPEND sources ${directory}/${i})
endforeach (i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\TOPDOWN" FILES ${sources})

