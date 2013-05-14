
cdef extern from "<OpenMS/METADATA/ExperimentalSettings.h>" namespace "OpenMS":

    cdef cppclass ExperimentalSettings:

        ExperimentalSettings() nogil except +
        ExperimentalSettings(ExperimentalSettings) nogil except + # wrap-ignore
