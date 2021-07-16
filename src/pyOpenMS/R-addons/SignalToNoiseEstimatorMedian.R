
SignalToNoiseEstimatorMedianChrom <- R6::R6Class(classname="SignalToNoiseEstimatorMedianChrom",cloneable = FALSE,
    private = list(py_obj=NA),

    public = list(
      initialize = function(in_0){
        if(missing(in_0)){
          private$$py_obj <- Pymod$$SignalToNoiseEstimatorMedianChrom()
        } else {
          if (is.R6(in_0) && class(in_0)[1]=="SignalToNoiseEstimatorMedianChrom" ){
            private$$py_obj <- Pymod$$SignalToNoiseEstimatorMedianChrom(r_to_py(in_0))
          } else {
            stop("arg in_0 wrong type")
          }
        }
      },
    init = function(chromatogram){
      if (!(is.R6(chromatogram) && class(chromatogram)[1]=="MSChromatogram")){ stop("arg chromatogram wrong type") }
        private$$py_obj$$init(r_to_py(chromatogram))
        invisible()
    },
    getSignalToNoise = function(data_point){
      if(!(is.R6(data_point) && class(data_point)[1]=="ChromatogramPeak")) { stop("arg data_point wrong type") }
      private$$py_obj$$getSignalToNoise(r_to_py(data_point))
    }
    )
)