from UniqueIdInterface cimport setUniqueId as _setUniqueId
import numpy as np
cimport numpy as np
from cython.parallel cimport prange


    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))


    def get_intensity_arr(self):
        """Generates a numpy ndarray with feature intensities from each sample in long format (over files).

        For labelled analyses channel intensities will be in one row, therefore resulting in a semi-long/block format.
        Resulting DataFrame can be joined with result from get_metadata_df by their index 'id'.

        Returns:
        np.ndarray: intensity ndarray
        """

        cdef bool labelfree
        labelfree = self.getExperimentType() == "label-free"
        cdef libcpp_map[int, _ColumnHeader] filemeta
        filemeta = self.inst.get().getColumnHeaders()
        s = self.size()

        ## TODO We could check that all labels are the same, otherwise a tabular export is a bit undefined.
        ## But it should work, since we will build a union of all labels, and only the labels that have a value
        ## in the current line (corresponding to a file) will be populated.

        labels = list(set([header.label for header in filemeta.values()]))
        files = list(set([header.filename for header in filemeta.values()]))
        label_to_idx = {k: v for v, k in enumerate(labels)}
        file_to_idx = {k: v for v, k in enumerate(files)}

        def gen(cmap: _ConsensusMap, fun):
            for f in cmap:
                yield from fun(f)

        if not labelfree:

            cdef np.uint64_t[:] idcol
            idcol = np.empty(s, dtype=np.uint64)

            cdef np.float64_t[:, ::1] labelcols
            labelcols = np.empty((s, len(labels)), dtype=np.float64)

            cdef np.ndarray[np.unicode_[300], ndim=1] filecol
            filecol = np.empty(s, dtype='U300')

            with nogil:
                cdef libcpp_vector[_FeatureHandle] subfeatures
                cdef libcpp_vector[_FeatureHandle].iterator fit
                cdef _ConsensusMap cmap
                cdef _FeatureHandle elm
                cmap = self.inst.get()
                for i in prange(s):
                    subfeatures = cmap[i].getFeatureList()
                    fit = subfeatures.begin()

                    int j = 0
                    while fit != subfeatures.end():
                        elm = &deref(fit)
                        
                        idcol[i+j] = cmap[i].getUniqueId()
                        labelcols[i]
                        j++



#### OLD CODE
            def extract_row_blocks_channel_wide_file_long(f: _ConsensusFeature):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                filerows = _defaultdict(lambda: [0] * len(labels))
                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row = filerows[header.filename]
                    row[label_to_idx[header.label]] = fh.getIntensity()
                return (f.getUniqueId(), filerows)

            def extract_rows_channel_wide_file_long(f: _ConsensusFeature):
                uniqueid, rowdict = extract_row_blocks_channel_wide_file_long(f)
                for file, row in rowdict.items():
                    row.append(file)
                    yield tuple([uniqueid] + row)

            if len(labels) == 1:
                labels[0] = "intensity"


            #### TODO: SHOULDNT WE MULTIPLY COUNT BY NUMBER OF FILES? Cause each file will produce a row per feature.

            intyarr = np.fromiter(iter=gen(self, extract_rows_channel_wide_file_long), dtype=dtypes, count=self.size())

            return intyarr

        else:
            # Specialized for LabelFree which has to have only one channel
            def extract_row_blocks_channel_long_file_wide_LF(f: _ConsensusFeature):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                row = [0.] * len(files)

                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row[file_to_idx[header.filename]] = fh.getIntensity()

                yield tuple([f.getUniqueId()] + row)

            dtypes = [('id', np.uint64)] + list(zip(files, [np.float64] * len(files)))

            intyarr = np.fromiter(iter=gen(self, extract_row_blocks_channel_long_file_wide_LF), dtype=dtypes, count=self.size())

            return intyarr
