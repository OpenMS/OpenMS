### the directory name
set(directory include/OpenMS/QC)

### list all header files of the directory here
set(sources_list_h
  Contaminants.h
  FragmentMassError.h
  FWHM.h
  MissedCleavages.h
  Ms2IdentificationRate.h
  Ms2SpectrumStats.h
  MzCalibration.h
  PeptideMass.h
  QCBase.h
  RTAlignment.h
  TIC.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
  list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\QC" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
