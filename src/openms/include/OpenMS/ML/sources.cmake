### the directory name
set(directory include/OpenMS/ML)

### list all header files of the directory here
set(sources_list_h
    # Constants mirrored from peptdeep plus small helpers over them. Header-only and free of
    # ONNX, and ProSE registers the same 'peptdeep:instrument' values with and without ONNX
    # support, so an ini file stays portable between the two builds.
    PEPTDEEP/PeptDeepUtils.h
)

if (WITH_ONNX)
    list(APPEND sources_list_h
        ONNX/ONNXPredictorBase.h
        PEPTDEEP/PeptDeepCCSInference.h
        PEPTDEEP/PeptDeepInput.h
        PEPTDEEP/PeptDeepMS2Inference.h
        PEPTDEEP/PeptDeepRTInference.h
    )
endif()

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
  list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ML" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
