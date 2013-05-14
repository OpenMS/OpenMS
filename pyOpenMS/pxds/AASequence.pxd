from libcpp cimport bool
from Types cimport *
from String cimport *
from Residue cimport *

cdef extern from "<OpenMS/CHEMISTRY/AASequence.h>" namespace "OpenMS":

    cdef cppclass AASequence:

        AASequence() nogil except +
        AASequence(AASequence) nogil except + # wrap-ignore

        AASequence(char *) nogil except +

        String toString()  nogil except +
        String toUnmodifiedString()  nogil except +

        void setModification(Size index, String modification) nogil except +
        void setNTerminalModification(String modification) nogil except +
        void setCTerminalModification(String modification) nogil except +
        void setStringSequence(String modification) nogil except +

        String getNTerminalModification() nogil except +
        String getCTerminalModification() nogil except +

        AASequence operator+(AASequence)    nogil except +
        AASequence iadd(AASequence)   nogil except + # wrap-as:operator+=




"""
TODO: Found function in cpp but not in pxd: function public getResidue
TODO: Found function in cpp but not in pxd: function public getResidue
TODO: Found function in cpp but not in pxd: function public getFormula
TODO: Found function in cpp but not in pxd: function public getAverageWeight
TODO: Found function in cpp but not in pxd: function public getMonoWeight
TODO: Found function in cpp but not in pxd: function public operator[]
TODO: Found function in cpp but not in pxd: function public operator[]
TODO: Found function in cpp but not in pxd: function public size
TODO: Found function in cpp but not in pxd: function public getPrefix
TODO: Found function in cpp but not in pxd: function public getSuffix
TODO: Found function in cpp but not in pxd: function public getSubsequence
TODO: Found function in cpp but not in pxd: function public getNumberOf
TODO: Found function in cpp but not in pxd: function public getAAFrequencies
TODO: Found function in cpp but not in pxd: function public isValid
TODO: Found function in cpp but not in pxd: function public has
TODO: Found function in cpp but not in pxd: function public has
TODO: Found function in cpp but not in pxd: function public hasSubsequence
TODO: Found function in cpp but not in pxd: function public hasSubsequence
TODO: Found function in cpp but not in pxd: function public hasPrefix
TODO: Found function in cpp but not in pxd: function public hasPrefix
TODO: Found function in cpp but not in pxd: function public hasSuffix
TODO: Found function in cpp but not in pxd: function public hasSuffix
TODO: Found function in cpp but not in pxd: function public hasNTerminalModification
TODO: Found function in cpp but not in pxd: function public hasCTerminalModification
TODO: Found function in cpp but not in pxd: function public isModified
TODO: Found function in cpp but not in pxd: function public isModified
TODO: Found function in cpp but not in pxd: function public operator<
TODO: Found function in cpp but not in pxd: function public begin
TODO: Found function in cpp but not in pxd: function public begin
TODO: Found function in cpp but not in pxd: function public end
TODO: Found function in cpp but not in pxd: function public end
TODO: Found function in cpp but not in pxd: function public empty

"""
