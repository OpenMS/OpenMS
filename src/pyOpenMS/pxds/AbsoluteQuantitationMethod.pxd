from Types cimport *
from Param cimport *

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethod:

        AbsoluteQuantitationMethod() except + nogil  # compiler
        AbsoluteQuantitationMethod(AbsoluteQuantitationMethod &) except + nogil  # compiler

        void setLLOD(double llod) except + nogil 
        void setULOD(double ulod) except + nogil 
        double getLLOD() except + nogil 
        double getULOD() except + nogil 

        void setLLOQ(double lloq) except + nogil 
        void setULOQ(double uloq) except + nogil 
        double getLLOQ() except + nogil 
        double getULOQ() except + nogil 

        bool checkLOD(double value) except + nogil 
        bool checkLOQ(double value) except + nogil 
        
        void setComponentName(const String& component_name) except + nogil 
        void setISName(const String& IS_name) except + nogil 
        void setFeatureName(const String& feature_name) except + nogil 
        String getComponentName() except + nogil 
        String getISName() except + nogil 
        String getFeatureName() except + nogil 

        void setConcentrationUnits(const String& concentration_units) except + nogil 
        String getConcentrationUnits() except + nogil 

        void setTransformationModel(const String& transformation_model) except + nogil 
        void setTransformationModelParams(Param transformation_model_param) except + nogil 
        String getTransformationModel() except + nogil 
        Param getTransformationModelParams() except + nogil 

        void setNPoints(Int n_points) except + nogil 
        void setCorrelationCoefficient(double correlation_coefficient) except + nogil 
        Int getNPoints() except + nogil 
        double getCorrelationCoefficient() except + nogil 
