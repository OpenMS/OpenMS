from libcpp cimport bool
from FeatureMap cimport *

cdef extern from "<OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>" namespace "OpenMS":

    cdef enum MergeIntensityMode:
        # wrap-attach:
        #    FeatureOverlapFilter
        SUM  # Sum intensities of merged features
        MAX  # Keep maximum intensity

cdef extern from "<OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>" namespace "OpenMS":

    cdef cppclass FeatureOverlapFilter:
        # wrap-doc:
        #    Filters and merges overlapping features in a FeatureMap.
        #
        #    This class provides functionality for detecting and handling overlapping features
        #    based on various criteria including convex hull overlap, trace-level overlap,
        #    and centroid-based distance thresholds.

        FeatureOverlapFilter() except + nogil
        FeatureOverlapFilter(FeatureOverlapFilter &) except + nogil  # compiler

# wrap static methods
cdef extern from "<OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>" namespace "OpenMS::FeatureOverlapFilter":

    void mergeOverlappingFeatures(FeatureMap & feature_map,
                                  double max_rt_diff,
                                  double max_mz_diff,
                                  bool require_same_charge,
                                  bool require_same_im,
                                  MergeIntensityMode intensity_mode,
                                  bool write_meta_values) except + nogil
        # wrap-attach:
        #    FeatureOverlapFilter
        # wrap-doc:
        #    Merge overlapping features based on centroid distances.
        #
        #    Identifies features whose centroids are within the specified RT and m/z tolerances
        #    and merges them into a single representative feature. The feature with highest
        #    intensity is kept as the representative.
        #
        #    When write_meta_values is True, the following meta values are stored on merged features:
        #    - "merged_centroid_rts": RT positions of all features that were merged
        #    - "merged_centroid_mzs": m/z positions of all features that were merged
        #    - "merged_centroid_IMs": FAIMS CV values (if present on the original features)
        #    - "FAIMS_merge_count": number of FAIMS CV values that were merged
        #
        #    :param feature_map: The feature map to process (modified in place)
        #    :param max_rt_diff: Maximum RT difference in seconds for considering features as overlapping
        #    :param max_mz_diff: Maximum m/z difference in Da for considering features as overlapping
        #    :param require_same_charge: If True, only merge features with identical charge states
        #    :param require_same_im: If True, only merge features with identical FAIMS CV values. Features without FAIMS_CV are treated as a separate group
        #    :param intensity_mode: How to combine intensities (MergeIntensityMode.SUM or MergeIntensityMode.MAX)
        #    :param write_meta_values: If True, write merge tracking meta values to merged features
