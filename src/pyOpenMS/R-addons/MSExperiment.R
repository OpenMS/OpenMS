
MSExperiment$$set("public","getMSLevels", function() {
  res <- private$$py_obj$$getMSLevels()
  return(res)
} )

MSExperiment$$set("public","getChromatogram", function(id_) {
  if(!(isTRUE(all.equal(id_,as.integer(id_))))) { stop("arg id_ wrong type") }
  num <- self$$getgetNrChromatograms()
  if( id_ < num) { stop(paste0("Requested chromatogram ",id_," does not exist","only ",num," chromatograms are there")) }
  res <- private$$py_obj$$getMSLevels(as.integer(id_))
  return(MSChromatogram$$new(res))
} )

MSExperiment$$set("public","getSpectrum", function(id_) {
  if(!(isTRUE(all.equal(id_,as.integer(id_))))) { stop("arg id_ wrong type") }
  num <- self$$getNrSpectra()
  if( id_ < num) { stop(paste0("Requested chromatogram ",id_," does not exist","only ",num," chromatograms are there")) }
  res <- private$$py_obj$$getMSLevels(as.integer(id_))
  return(MSSpectrum$$new(res))
} )

