from MSSpectrum cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeWavelet "OpenMS::IsotopeWavelet":
        IsotopeWavelet(IsotopeWavelet) nogil except + #wrap-ignore
        # IsotopeWavelet * init(double max_m, UInt max_charge) nogil except +
        # IsotopeWavelet * getInstance() nogil except +
        void destroy() nogil except +
        double getValueByMass(double t, double m, UInt z, Int mode) nogil except +
        double getValueByLambda(double lambda_, double tz1) nogil except +
        double getValueByLambdaExtrapol(double lambda_, double tz1) nogil except +
        double getValueByLambdaExact(double lambda_, double tz1) nogil except +
        UInt getMaxCharge() nogil except +
        void setMaxCharge(UInt max_charge) nogil except +
        double getTableSteps() nogil except +
        double getInvTableSteps() nogil except +
        void setTableSteps(double table_steps) nogil except +
        double getLambdaL(double m) nogil except +
        # IsotopeDistribution::ContainerType  getAveragine(double m, UInt *size) nogil except +
        Size getGammaTableMaxIndex() nogil except +
        Size getExpTableMaxIndex() nogil except +
        float myPow(float a, float b) nogil except +
        UInt getMzPeakCutOffAtMonoPos(double mass, UInt z) nogil except +
        UInt getNumPeakCutOff(double mass, UInt z) nogil except +
        UInt getNumPeakCutOff(double mz) nogil except +

