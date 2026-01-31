from libcpp.vector cimport vector as libcpp_vector
from libc.stdint cimport int64_t
from cython.operator cimport dereference as deref
from XICParquetFile cimport XICChromatogram as _XICChromatogram
import numpy as np


    def get_data_dict(self, explode=False, precursor_id=-1, transition_id=-1, modified_sequence="", precursor_charge=-1,
                      product_charge=-1, ms_level=-1, run_id=-1, filter=""):
        """
        Return chromatogram data as a dict of numpy arrays.

        If explode=True, returns long format with rt/intensity rows.
        Otherwise, rt and intensity are stored as object arrays of lists.
        """
        cdef libcpp_vector[_XICChromatogram] chroms
        self.inst.get().getChromatograms(chroms,
                                         <int64_t>precursor_id,
                                         <int64_t>transition_id,
                                         deref((convString(modified_sequence)).get()),
                                         <int64_t>precursor_charge,
                                         <int64_t>product_charge,
                                         <int64_t>ms_level,
                                         <int64_t>run_id,
                                         deref((convString(filter)).get()))

        cdef list run_id_list = []
        cdef list source_file_list = []
        cdef list ms_level_list = []
        cdef list precursor_id_list = []
        cdef list transition_id_list = []
        cdef list modified_sequence_list = []
        cdef list precursor_charge_list = []
        cdef list product_charge_list = []
        cdef list detecting_transition_list = []
        cdef list precursor_decoy_list = []
        cdef list product_decoy_list = []
        cdef list transition_ordinal_list = []
        cdef list transition_type_list = []
        cdef list annotation_list = []
        cdef list rt_list = []
        cdef list intensity_list = []

        cdef size_t i, j
        cdef _XICChromatogram chrom
        for i in range(chroms.size()):
            chrom = chroms[i]
            if explode:
                if chrom.rt.size() == 0:
                    continue
                for j in range(chrom.rt.size()):
                    run_id_list.append(chrom.run_id)
                    source_file_list.append(convOutputString(chrom.source_file))
                    ms_level_list.append(chrom.ms_level)
                    precursor_id_list.append(chrom.precursor_id if chrom.has_precursor_id else None)
                    transition_id_list.append(chrom.transition_id if chrom.has_transition_id else None)
                    modified_sequence_list.append(convOutputString(chrom.modified_sequence))
                    precursor_charge_list.append(chrom.precursor_charge if chrom.has_precursor_charge else None)
                    product_charge_list.append(chrom.product_charge if chrom.has_product_charge else None)
                    detecting_transition_list.append(chrom.detecting_transition if chrom.has_detecting_transition else None)
                    precursor_decoy_list.append(chrom.precursor_decoy if chrom.has_precursor_decoy else None)
                    product_decoy_list.append(chrom.product_decoy if chrom.has_product_decoy else None)
                    transition_ordinal_list.append(chrom.transition_ordinal if chrom.has_transition_ordinal else None)
                    transition_type_list.append(convOutputString(chrom.transition_type))
                    annotation_list.append(convOutputString(chrom.annotation))
                    rt_list.append(chrom.rt[j])
                    intensity_list.append(chrom.intensity[j])
            else:
                run_id_list.append(chrom.run_id)
                source_file_list.append(convOutputString(chrom.source_file))
                ms_level_list.append(chrom.ms_level)
                precursor_id_list.append(chrom.precursor_id if chrom.has_precursor_id else None)
                transition_id_list.append(chrom.transition_id if chrom.has_transition_id else None)
                modified_sequence_list.append(convOutputString(chrom.modified_sequence))
                precursor_charge_list.append(chrom.precursor_charge if chrom.has_precursor_charge else None)
                product_charge_list.append(chrom.product_charge if chrom.has_product_charge else None)
                detecting_transition_list.append(chrom.detecting_transition if chrom.has_detecting_transition else None)
                precursor_decoy_list.append(chrom.precursor_decoy if chrom.has_precursor_decoy else None)
                product_decoy_list.append(chrom.product_decoy if chrom.has_product_decoy else None)
                transition_ordinal_list.append(chrom.transition_ordinal if chrom.has_transition_ordinal else None)
                transition_type_list.append(convOutputString(chrom.transition_type))
                annotation_list.append(convOutputString(chrom.annotation))
                rt_list.append([chrom.rt[k] for k in range(chrom.rt.size())])
                intensity_list.append([chrom.intensity[k] for k in range(chrom.intensity.size())])

        return {
            "run_id": np.array(run_id_list, dtype=np.int64),
            "source_file": np.array(source_file_list, dtype=object),
            "ms_level": np.array(ms_level_list, dtype=np.int64),
            "precursor_id": np.array(precursor_id_list, dtype=object),
            "transition_id": np.array(transition_id_list, dtype=object),
            "modified_sequence": np.array(modified_sequence_list, dtype=object),
            "precursor_charge": np.array(precursor_charge_list, dtype=object),
            "product_charge": np.array(product_charge_list, dtype=object),
            "detecting_transition": np.array(detecting_transition_list, dtype=object),
            "precursor_decoy": np.array(precursor_decoy_list, dtype=object),
            "product_decoy": np.array(product_decoy_list, dtype=object),
            "transition_ordinal": np.array(transition_ordinal_list, dtype=object),
            "transition_type": np.array(transition_type_list, dtype=object),
            "annotation": np.array(annotation_list, dtype=object),
            "rt": np.array(rt_list, dtype=object),
            "intensity": np.array(intensity_list, dtype=object),
        }

    def to_pandas(self, explode=False, **kwargs):
        """
        Return chromatogram data as a pandas DataFrame.

        If explode=True, returns long format with rt/intensity rows.
        """
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError("pandas is required for this method. Install with `pip install pandas`.") from e
        return pd.DataFrame(self.get_data_dict(explode=explode, **kwargs))

    def get_chromatograms(self, explode=True, **kwargs):
        """
        Return filtered chromatograms as a pandas DataFrame (exploded by default).
        """
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError("pandas is required for this method. Install with `pip install pandas`.") from e
        return pd.DataFrame(self.get_data_dict(explode=explode, **kwargs))
