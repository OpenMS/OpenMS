from String cimport *
from Gradient cimport *

cdef extern from "<OpenMS/METADATA/HPLC.h>" namespace "OpenMS":

    cdef cppclass HPLC:

        HPLC() nogil except +
        HPLC(HPLC) nogil except + # wrap-ignore

        # returns a reference to the instument name
        String getInstrument() nogil except +
        # sets the instument name
        void setInstrument(String instrument) nogil except +

        # returns a reference to the column description
        String getColumn() nogil except +
        # sets the column description
        void setColumn(String column) nogil except +

        # returns the temperature (in degree C)
        Int getTemperature() nogil except +
        # sets the temperature (in degree C)
        void setTemperature(Int temperature) nogil except +

        # returns the pressure (in bar)
        UInt getPressure() nogil except +
        # sets the pressure (in bar)
        void setPressure(UInt pressure) nogil except +

        # returns the flux (in microliter/sec)
        UInt getFlux() nogil except +
        # sets the flux (in microliter/sec)
        void setFlux(UInt flux) nogil except +

        # returns the comments
        String getComment() nogil except +
        # sets the comments
        void setComment(String comment) nogil except +

        # returns a mutable reference to the used gradient
        Gradient getGradient() nogil except +
        # sets the used gradient
        void setGradient(Gradient gradient) nogil except +


