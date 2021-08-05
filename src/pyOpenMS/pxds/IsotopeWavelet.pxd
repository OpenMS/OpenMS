from MSSpectrum cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeWavelet "OpenMS::IsotopeWavelet":

        #protected
        IsotopeWavelet() nogil except + #wrap-ignore
        #protected
        IsotopeWavelet(IsotopeWavelet &) nogil except + #wrap-ignore

        # IsotopeWavelet * init(double max_m, UInt max_charge) nogil except +
        # IsotopeWavelet * getInstance() nogil except +
        void destroy() nogil except + # wrap-doc:Deletes the singleton instance
        double getValueByMass(double t, double m, UInt z, Int mode) nogil except +
            # wrap-doc:
                #   Returns the value of the isotope wavelet at position `t`. Usually, you do not need to call this function
                #   -----
                #   Note that this functions returns the pure function value of psi and not the normalized (average=0)
                #   value given by Psi
                #   -----
                #   :param t: The position at which the wavelet has to be drawn (within the coordinate system of the wavelet)
                #   :param m: The m/z position within the signal (i.e. the mass not de-charged) within the signal
                #   :param z: The charge `z` we want to detect
                #   :param mode: Indicates whether positive mode (+1) or negative mode (-1) has been used for ionization

        double getValueByLambda(double lambda_, double tz1) nogil except +
            # wrap-doc:
                #   Returns the value of the isotope wavelet at position `t` via a fast table lookup
                #   -----
                #   Usually, you do not need to call this function
                #   Please use `sampleTheWavelet` instead
                #   Note that this functions returns the pure function value of psi and not the normalized (average=0)
                #   value given by Psi
                #   -----
                #   :param lambda: The mass-parameter lambda
                #   :param tz1: t (the position) times the charge (z) plus 1

        double getValueByLambdaExtrapol(double lambda_, double tz1) nogil except +
            # wrap-doc:
                #   Returns the value of the isotope wavelet at position `t`
                #   -----
                #   This function is usually significantly slower than the table lookup performed in @see getValueByLambda
                #   Nevertheless, it might be necessary to call this function due to extrapolating reasons caused by the
                #   alignment of the wavelet
                #   -----
                #   Usually, you do not need to call this function
                #   Please use `sampleTheWavelet` instead
                #   Note that this functions returns the pure function value of psi and not the normalized (average=0)
                #   value given by Psi
                #   -----
                #   :param lambda: The mass-parameter lambda
                #   :param tz1: t (the position) times the charge (z) plus 1

        double getValueByLambdaExact(double lambda_, double tz1) nogil except +# TODO
        UInt getMaxCharge() nogil except + # wrap-doc:Returns the largest charge state we will consider
        void setMaxCharge(UInt max_charge) nogil except + # wrap-doc:Sets the `max_charge` parameter
        double getTableSteps() nogil except + # wrap-doc:Returns the table_steps_ parameter
        double getInvTableSteps() nogil except + # wrap-doc:Returns the inv_table_steps_ parameter
        void setTableSteps(double table_steps) nogil except + # wrap-doc:Sets the `table_steps` parameter
        double getLambdaL(double m) nogil except + # wrap-doc:Returns the mass-parameter lambda (linear fit)
        # IsotopeDistribution::ContainerType  getAveragine(double m, UInt *size) nogil except +
        Size getGammaTableMaxIndex() nogil except + # wrap-doc:Returns the largest possible index for the pre-sampled gamma table
        Size getExpTableMaxIndex() nogil except + # wrap-doc:Returns the largest possible index for the pre-sampled exp table
        float myPow(float a, float b) nogil except + # wrap-doc:Internally used function; uses register shifts for fast computation of the power function
        UInt getMzPeakCutOffAtMonoPos(double mass, UInt z) nogil except +# TODO 
        UInt getNumPeakCutOff(double mass, UInt z) nogil except +# TODO
        UInt getNumPeakCutOff(double mz) nogil except +# TODO
