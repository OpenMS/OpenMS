from Types cimport *
from String cimport *
from NASequence cimport *
from IdentificationData cimport *
from EnzymaticDigestion cimport *

cdef extern from "<OpenMS/CHEMISTRY/RNaseDigestion.h>" namespace "OpenMS":

    cdef cppclass RNaseDigestion(EnzymaticDigestion):
        # wrap-inherits:
        #    EnzymaticDigestion
        #
        # wrap-doc:
        #     Class for the enzymatic digestion of RNA
        #     -----
        #     Usage:
        #         from pyopenms import *
        #         oligo = NASequence.fromString("pAUGUCGCAG");
        #         -
        #         dig = RNaseDigestion()
        #         dig.setEnzyme("RNase_T1")
        #         -
        #         result = []
        #         dig.digest(oligo, result)
        #         for fragment in result:
        #           print (fragment)

      RNaseDigestion() nogil except + # compiler
      RNaseDigestion(RNaseDigestion &) nogil except + # compiler

      void setEnzyme(String name) nogil except + # wrap-doc:Sets the enzyme for the digestion (by name)

      void digest(NASequence & rna, libcpp_vector[ NASequence ] & output) nogil except +

      void digest(NASequence & rna, libcpp_vector[ NASequence ] & output, Size min_length, Size max_length) nogil except +
          # wrap-doc:
          #     Performs the enzymatic digestion of a (potentially modified) RNA
          #     -----
          #     :param rna: Sequence to digest
          #     :param output: Digestion productsq
          #     :param min_length: Minimal length of reported products
          #     :param max_length: Maximal length of reported products (0 = no restriction)
          #     :returns: Number of discarded digestion products (which are not matching length restrictions)

      void digest(IdentificationData & id_data) nogil except +

      void digest(IdentificationData & id_data, Size min_length, Size max_length) nogil except +
          # wrap-doc:
          #     Performs the enzymatic digestion of all RNA parent molecules in IdentificationData (id_data)
          #     -----
          #     :param id_data: IdentificationData object which includes sequences to digest
          #     :param min_length: Minimal length of reported products
          #     :param max_length: Maximal length of reported products (0 = no restriction)
          #     :returns: Number of discarded digestion products (which are not matching length restrictions)
