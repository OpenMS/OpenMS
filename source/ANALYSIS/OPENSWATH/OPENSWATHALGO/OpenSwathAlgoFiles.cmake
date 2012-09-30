### the directory name
set(directory source/ANALYSIS/OPENSWATH/OPENSWATHALGO)

### list all files of the directory here
set(sources_list
  ALGO/Scoring.C
  ALGO/MRMScoring.C
  ALGO/StatsHelpers.C
  DATAACCESS/SpectrumHelpers.C
)

### add path to the filenames
set(OpenSwathAlgoFiles)
foreach(i ${sources_list})
	list(APPEND OpenSwathAlgoFiles ${directory}/${i})
endforeach(i)

