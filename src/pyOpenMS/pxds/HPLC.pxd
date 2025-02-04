from String cimport *
from Gradient cimport *

cdef extern from "<OpenMS/METADATA/HPLC.h>" namespace "OpenMS":

    cdef cppclass HPLC:

        HPLC() except + nogil  # wrap-doc:Representation of a HPLC experiment
        HPLC(HPLC &) except + nogil 

    
        String getInstrument() except + nogil  # wrap-doc:Returns a reference to the instument name
    
        void setInstrument(String instrument) except + nogil  # wrap-doc:Sets the instument name

    
        String getColumn() except + nogil  # wrap-doc:Returns a reference to the column description
    
        void setColumn(String column) except + nogil  # wrap-doc:Sets the column description

    
        Int getTemperature() except + nogil  # wrap-doc:Returns the temperature (in degree C)
    
        void setTemperature(Int temperature) except + nogil  # wrap-doc:Sets the temperature (in degree C)

    
        UInt getPressure() except + nogil  # wrap-doc:Returns the pressure (in bar)
    
        void setPressure(UInt pressure) except + nogil  # wrap-doc:Sets the pressure (in bar)

    
        UInt getFlux() except + nogil  # wrap-doc:Returns the flux (in microliter/sec)
    
        void setFlux(UInt flux) except + nogil  # wrap-doc:Sets the flux (in microliter/sec)

    
        String getComment() except + nogil  # wrap-doc:Returns the comments
    
        void setComment(String comment) except + nogil  # wrap-doc:Sets the comments

    
        Gradient getGradient() except + nogil  # wrap-doc:Returns a mutable reference to the used gradient
    
        void setGradient(Gradient gradient) except + nogil  # wrap-doc:Sets the used gradient


