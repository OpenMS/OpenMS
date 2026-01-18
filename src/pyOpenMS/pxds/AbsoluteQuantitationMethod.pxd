from Types cimport *
from Param cimport *

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethod:
        # wrap-doc:
        #  Holds information about a quantitation method for absolute quantitation using
        #  Isotope Dilution Mass Spectrometry (IDMS). This includes calibration curve
        #  parameters, limits of detection (LOD), limits of quantitation (LOQ), and
        #  the transformation model used for concentration calculation

        AbsoluteQuantitationMethod() except + nogil  # compiler
        AbsoluteQuantitationMethod(AbsoluteQuantitationMethod &) except + nogil  # compiler

        void setLLOD(double llod) except + nogil  # wrap-doc:Sets the lower limit of detection (LLOD)
        void setULOD(double ulod) except + nogil  # wrap-doc:Sets the upper limit of detection (ULOD)
        double getLLOD() except + nogil  # wrap-doc:Returns the lower limit of detection (LLOD)
        double getULOD() except + nogil  # wrap-doc:Returns the upper limit of detection (ULOD)

        void setLLOQ(double lloq) except + nogil  # wrap-doc:Sets the lower limit of quantitation (LLOQ)
        void setULOQ(double uloq) except + nogil  # wrap-doc:Sets the upper limit of quantitation (ULOQ)
        double getLLOQ() except + nogil  # wrap-doc:Returns the lower limit of quantitation (LLOQ)
        double getULOQ() except + nogil  # wrap-doc:Returns the upper limit of quantitation (ULOQ)

        bool checkLOD(double value) except + nogil  # wrap-doc:Checks if the given value is within the limits of detection (LOD)
        bool checkLOQ(double value) except + nogil  # wrap-doc:Checks if the given value is within the limits of quantitation (LOQ)

        void setComponentName(const String& component_name) except + nogil  # wrap-doc:Sets the component name (unique identifier for the analyte)
        void setISName(const String& IS_name) except + nogil  # wrap-doc:Sets the internal standard name used for ratio calculation
        void setFeatureName(const String& feature_name) except + nogil  # wrap-doc:Sets the feature name (e.g., peak_apex_int or peak_area)
        String getComponentName() except + nogil  # wrap-doc:Returns the component name (unique identifier for the analyte)
        String getISName() except + nogil  # wrap-doc:Returns the internal standard name used for ratio calculation
        String getFeatureName() except + nogil  # wrap-doc:Returns the feature name (e.g., peak_apex_int or peak_area)

        void setConcentrationUnits(const String& concentration_units) except + nogil  # wrap-doc:Sets the concentration units for the component
        String getConcentrationUnits() except + nogil  # wrap-doc:Returns the concentration units for the component

        void setTransformationModel(const String& transformation_model) except + nogil  # wrap-doc:Sets the transformation model type (e.g., linear, b_spline)
        void setTransformationModelParams(Param transformation_model_param) except + nogil  # wrap-doc:Sets the transformation model parameters
        String getTransformationModel() except + nogil  # wrap-doc:Returns the transformation model type (e.g., linear, b_spline)
        Param getTransformationModelParams() except + nogil  # wrap-doc:Returns the transformation model parameters

        void setNPoints(Int n_points) except + nogil  # wrap-doc:Sets the number of points used in the calibration curve
        void setCorrelationCoefficient(double correlation_coefficient) except + nogil  # wrap-doc:Sets the Pearson correlation coefficient of the calibration curve fit
        Int getNPoints() except + nogil  # wrap-doc:Returns the number of points used in the calibration curve
        double getCorrelationCoefficient() except + nogil  # wrap-doc:Returns the Pearson correlation coefficient of the calibration curve fit
