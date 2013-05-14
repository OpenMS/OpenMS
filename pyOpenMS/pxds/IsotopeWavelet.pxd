from MSSpectrum cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeWavelet "OpenMS::IsotopeWavelet":
        IsotopeWavelet(IsotopeWavelet) nogil except + #wrap-ignore
        # IsotopeWavelet * init(DoubleReal max_m, UInt max_charge) nogil except +
        # IsotopeWavelet * getInstance() nogil except +
        void destroy() nogil except +
        DoubleReal getValueByMass(DoubleReal t, DoubleReal m, UInt z, Int mode) nogil except +
        DoubleReal getValueByLambda(DoubleReal lambda_, DoubleReal tz1) nogil except +
        DoubleReal getValueByLambdaExtrapol(DoubleReal lambda_, DoubleReal tz1) nogil except +
        DoubleReal getValueByLambdaExact(DoubleReal lambda_, DoubleReal tz1) nogil except +
        UInt getMaxCharge() nogil except +
        void setMaxCharge(UInt max_charge) nogil except +
        DoubleReal getTableSteps() nogil except +
        DoubleReal getInvTableSteps() nogil except +
        void setTableSteps(DoubleReal table_steps) nogil except +
        DoubleReal getLambdaL(DoubleReal m) nogil except +
        # IsotopeDistribution::ContainerType  getAveragine(DoubleReal m, UInt *size) nogil except +
        Size getGammaTableMaxIndex() nogil except +
        Size getExpTableMaxIndex() nogil except +
        float myPow(float a, float b) nogil except +
        UInt getMzPeakCutOffAtMonoPos(DoubleReal mass, UInt z) nogil except +
        UInt getNumPeakCutOff(DoubleReal mass, UInt z) nogil except +
        UInt getNumPeakCutOff(DoubleReal mz) nogil except +

