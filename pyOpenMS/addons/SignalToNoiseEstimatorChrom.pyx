# two empty lines are important !


cdef class SignalToNoiseEstimatorMedianChrom:

    cdef shared_ptr[_SignalToNoiseEstimatorMedian[_MSChromatogram[_ChromatogramPeak]]] inst

    def __dealloc__(self):
         self.inst.reset()

    
    def init(self, MSChromatogram chromatogram ):
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
    
        self.inst.get().init((deref(chromatogram.inst.get())))
    
    def _init_0(self):
        self.inst = shared_ptr[_SignalToNoiseEstimatorMedian[_MSChromatogram[_ChromatogramPeak]]](new _SignalToNoiseEstimatorMedian[_MSChromatogram[_ChromatogramPeak]]())
    
    def _init_1(self, SignalToNoiseEstimatorMedianChrom in_0 ):
        assert isinstance(in_0, SignalToNoiseEstimatorMedianChrom), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SignalToNoiseEstimatorMedian[_MSChromatogram[_ChromatogramPeak]]](new _SignalToNoiseEstimatorMedian[_MSChromatogram[_ChromatogramPeak]]((deref(in_0.inst.get()))))
    
    def __init__(self, *args):
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SignalToNoiseEstimatorMedian)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSignalToNoise(self, ChromatogramPeak data_point ):
        assert isinstance(data_point, ChromatogramPeak), 'arg data_point wrong type'
    
        self.inst.get().getSignalToNoise((deref(data_point.inst.get()))) 
