from Types cimport *
from Param cimport *

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethod:

        AbsoluteQuantitationMethod()  nogil except +
        AbsoluteQuantitationMethod(AbsoluteQuantitationMethod)  nogil except + #wrap-ignore

        void setLLOD(double llod) nogil except +
        void setULOD(double ulod) nogil except +
        double getLLOD() nogil except +
        double getULOD() nogil except +

        void setLLOQ(double lloq) nogil except +
        void setULOQ(double uloq) nogil except +
        double getLLOQ() nogil except +
        double getULOQ() nogil except +

        bool checkLOD(double value) nogil except +
        bool checkLOQ(double value) nogil except +
        
        void setComponentName(const String& component_name) nogil except +
        void setISName(const String& IS_name) nogil except +
        void setFeatureName(const String& feature_name) nogil except +
        String getComponentName() nogil except +
        String getISName() nogil except +
        String getFeatureName() nogil except +

        void setConcentrationUnits(const String& concentration_units) nogil except +
        String getConcentrationUnits() nogil except +

        void setTransformationModel(const String& transformation_model) nogil except +
        void setTransformationModelParams(Param transformation_model_param) nogil except +
        String getTransformationModel() nogil except +
        Param getTransformationModelParams() nogil except +

        void setNPoints(Int n_points) nogil except +
        void setCorrelationCoefficient(double correlation_coefficient) nogil except +
        Int getNPoints() nogil except +
        double getCorrelationCoefficient() nogil except +
