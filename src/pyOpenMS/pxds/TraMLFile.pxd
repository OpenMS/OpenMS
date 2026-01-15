from String cimport *
from StringList cimport *
from TargetedExperiment cimport *

cdef extern from "<OpenMS/FORMAT/TraMLFile.h>" namespace "OpenMS":

    cdef cppclass TraMLFile:
        # wrap-doc:
        #  File adapter for TraML files
        #
        #  Provides methods to load and store targeted experiment data in TraML format.
        #  TraML (Transition Markup Language) files store SRM/MRM transitions for targeted
        #  proteomics and metabolomics experiments.
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    targeted_exp = TargetedExperiment()
        #    TraMLFile().load("transitions.TraML", targeted_exp)
        #    for transition in targeted_exp.getTransitions():
        #      print(transition.getPrecursorMZ(), transition.getProductMZ())

        TraMLFile() except + nogil
        TraMLFile(TraMLFile &) except + nogil  # compiler

        void load(String filename,
                  TargetedExperiment & id) except + nogil  # wrap-doc:Loads a targeted experiment from a TraML file

        void store(String filename,
                  TargetedExperiment & id) except + nogil  # wrap-doc:Stores a targeted experiment in a TraML file

        bool isSemanticallyValid(String filename, StringList & errors,
                                 StringList & warnings) except + nogil
            # wrap-doc:
            #  Checks if a file is valid with respect to the mapping file and the controlled vocabulary
            #
            #  :param filename: File name of the file to be checked
            #  :param errors: Errors during the validation are returned in this output parameter
            #  :param warnings: Warnings during the validation are returned in this output parameter
            #  :returns: True if the file is valid, False otherwise
