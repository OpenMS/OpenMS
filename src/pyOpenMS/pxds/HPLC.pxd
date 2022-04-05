from String cimport *
from Gradient cimport *

cdef extern from "<OpenMS/METADATA/HPLC.h>" namespace "OpenMS":

    cdef cppclass HPLC:

        HPLC() nogil except + # wrap-doc:Representation of a HPLC experiment
        HPLC(HPLC &) nogil except +

    
        String getInstrument() nogil except + # wrap-doc:Returns a reference to the instument name
    
        void setInstrument(String instrument) nogil except + # wrap-doc:Sets the instument name

    
        String getColumn() nogil except + # wrap-doc:Returns a reference to the column description
    
        void setColumn(String column) nogil except + # wrap-doc:Sets the column description

    
        Int getTemperature() nogil except + # wrap-doc:Returns the temperature (in degree C)
    
        void setTemperature(Int temperature) nogil except + # wrap-doc:Sets the temperature (in degree C)

    
        UInt getPressure() nogil except + # wrap-doc:Returns the pressure (in bar)
    
        void setPressure(UInt pressure) nogil except + # wrap-doc:Sets the pressure (in bar)

    
        UInt getFlux() nogil except + # wrap-doc:Returns the flux (in microliter/sec)
    
        void setFlux(UInt flux) nogil except + # wrap-doc:Sets the flux (in microliter/sec)

    
        String getComment() nogil except + # wrap-doc:Returns the comments
    
        void setComment(String comment) nogil except + # wrap-doc:Sets the comments

    
        Gradient getGradient() nogil except + # wrap-doc:Returns a mutable reference to the used gradient
    
        void setGradient(Gradient gradient) nogil except + # wrap-doc:Sets the used gradient


