from MSSpectrum cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeWavelet "OpenMS::IsotopeWavelet":

        #protected
        IsotopeWavelet() except + nogil  #wrap-ignore
        #protected
        IsotopeWavelet(IsotopeWavelet &) except + nogil  #wrap-ignore

        # IsotopeWavelet * init(double max_m, UInt max_charge) except + nogil 
        # IsotopeWavelet * getInstance() except + nogil 
        void destroy() except + nogil  # wrap-doc:Deletes the singleton instance
        double getValueByMass(double t, double m, UInt z, Int mode) except + nogil 
            # wrap-doc:
                #  Returns the value of the isotope wavelet at position `t`. Usually, you do not need to call this function\n
                #  
                #  Note that this functions returns the pure function value of psi and not the normalized (average=0)
                #  value given by Psi
                #  
                #  
                #  :param t: The position at which the wavelet has to be drawn (within the coordinate system of the wavelet)
                #  :param m: The m/z position within the signal (i.e. the mass not de-charged) within the signal
                #  :param z: The charge `z` we want to detect
                #  :param mode: Indicates whether positive mode (+1) or negative mode (-1) has been used for ionization

        double getValueByLambda(double lambda_, double tz1) except + nogil 
            # wrap-doc:
                #  Returns the value of the isotope wavelet at position `t` via a fast table lookup\n
                #  
                #  Usually, you do not need to call this function
                #  Please use `sampleTheWavelet` instead
                #  Note that this functions returns the pure function value of psi and not the normalized (average=0)
                #  value given by Psi
                #  
                #  
                #  :param lambda: The mass-parameter lambda
                #  :param tz1: t (the position) times the charge (z) plus 1

        double getValueByLambdaExtrapol(double lambda_, double tz1) except + nogil 
            # wrap-doc:
                #  Returns the value of the isotope wavelet at position `t`\n
                #  
                #  This function is usually significantly slower than the table lookup performed in @see getValueByLambda
                #  Nevertheless, it might be necessary to call this function due to extrapolating reasons caused by the
                #  alignment of the wavelet\n
                #  
                #  Usually, you do not need to call this function
                #  Please use `sampleTheWavelet` instead
                #  Note that this functions returns the pure function value of psi and not the normalized (average=0)
                #  value given by Psi
                #  
                #  
                #  :param lambda: The mass-parameter lambda
                #  :param tz1: t (the position) times the charge (z) plus 1

        double getValueByLambdaExact(double lambda_, double tz1) except + nogil # TODO
        UInt getMaxCharge() except + nogil  # wrap-doc:Returns the largest charge state we will consider
        void setMaxCharge(UInt max_charge) except + nogil  # wrap-doc:Sets the `max_charge` parameter
        double getTableSteps() except + nogil  # wrap-doc:Returns the table_steps_ parameter
        double getInvTableSteps() except + nogil  # wrap-doc:Returns the inv_table_steps_ parameter
        void setTableSteps(double table_steps) except + nogil  # wrap-doc:Sets the `table_steps` parameter
        double getLambdaL(double m) except + nogil  # wrap-doc:Returns the mass-parameter lambda (linear fit)
        # IsotopeDistribution::ContainerType  getAveragine(double m, UInt *size) except + nogil 
        Size getGammaTableMaxIndex() except + nogil  # wrap-doc:Returns the largest possible index for the pre-sampled gamma table
        Size getExpTableMaxIndex() except + nogil  # wrap-doc:Returns the largest possible index for the pre-sampled exp table
        float myPow(float a, float b) except + nogil  # wrap-doc:Internally used function; uses register shifts for fast computation of the power function
        UInt getMzPeakCutOffAtMonoPos(double mass, UInt z) except + nogil # TODO 
        UInt getNumPeakCutOff(double mass, UInt z) except + nogil # TODO
        UInt getNumPeakCutOff(double mz) except + nogil # TODO