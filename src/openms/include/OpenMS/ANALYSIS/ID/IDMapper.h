// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <algorithm>
#include <limits>

namespace OpenMS
{
  /**
    @brief Annotates an MSExperiment, FeatureMap or ConsensusMap with peptide identifications

    ProteinIdentifications are assigned to the whole map.

    The retention time and mass-to-charge ratio of the PeptideIdentification as the respective "MZ" and "RT" members.

    m/z-matching on peptide side can be done either with the precursor m/z value of the peptide identification or with the theoretical masses of the peptide hits (see "mz_reference" parameter).

    See the documentation of the individual @p annotate methods for more in-depth information.

    @htmlinclude OpenMS_IDMapper.parameters
  */
  class OPENMS_DLLAPI IDMapper :
    public DefaultParamHandler
  {
public:
    enum Measure {MEASURE_PPM = 0, MEASURE_DA};

    /// Default constructor
    IDMapper();

    /// Copy C'Tor
    IDMapper(const IDMapper& cp);

    /// Assignment
    IDMapper& operator=(const IDMapper& rhs);

    /**
      @brief Mapping method for peak maps

      The identifications stored in a PeptideIdentification instance can be added to the
      corresponding spectrum.
      Note that a PeptideIdentication is added to ALL spectra which are within the allowed RT and MZ boundaries.

      @param map MSExperiment to receive the identifications
      @param peptide_ids PeptideIdentification for the MSExperiment
      @param protein_ids ProteinIdentification for the MSExperiment
      @param clear_ids Reset peptide and protein identifications of each scan before annotating
      @param map_ms1 Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)

      @exception Exception::MissingInformation is thrown if entries of @p peptide_ids do not contain 'MZ' and 'RT' information.
    */
    void annotate(PeakMap& map, const std::vector<PeptideIdentification>& peptide_ids, const std::vector<ProteinIdentification>& protein_ids, const bool clear_ids = false, const bool map_ms1 = false);

    /**
      @brief Mapping method for peak maps

      Add peptide identifications stored in a feature map to their
      corresponding spectrum.

      This function converts the feature map to a vector of peptide identifications (all peptide IDs from each feature are taken)
      and calls the respective annotate() function.
      RT and m/z are taken from the peptides, or (if missing) from the feature itself.

      @param map MSExperiment to receive the identifications
      @param fmap FeatureMap with PeptideIdentifications for the MSExperiment
      @param clear_ids Reset peptide and protein identifications of each scan before annotating
      @param map_ms1 attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
    */
    void annotate(PeakMap& map, FeatureMap fmap, const bool clear_ids = false, const bool map_ms1 = false);

    /**
      @brief Mapping method for feature maps

      If @em all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls by default. If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes) if @p use_centroid_rt or @p use_centroid_mz are true.

      In any case, tolerance in RT and m/z dimension is applied according to the global parameters @p rt_tolerance and @p mz_tolerance. Tolerance is understood as "plus or minus x", so the matching range is actually increased by twice the tolerance value.

      If several features (incl. tolerance) overlap the position of a peptide identification, the identification is annotated to all of them.

      @param map FeatureMap to receive the identifications
      @param ids PeptideIdentification for the ConsensusFeatures
      @param protein_ids ProteinIdentification for the ConsensusMap
      @param use_centroid_rt Whether to use the RT value of feature centroids even if convex hulls are present
      @param use_centroid_mz Whether to use the m/z value of feature centroids even if convex hulls are present
      @param spectra Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index.

      @exception Exception::MissingInformation is thrown if entries of @p ids do not contain 'MZ' and 'RT' information.
    */
    void annotate(FeatureMap& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool use_centroid_rt = false, bool use_centroid_mz = false, const PeakMap& spectra = PeakMap());

    /**
      @brief Mapping method for consensus maps

      If several consensus features lie inside the allowed deviation, the peptide identifications
      are mapped to all the consensus features.

      @param map ConsensusMap to receive the identifications
      @param ids PeptideIdentification for the ConsensusFeatures
      @param protein_ids ProteinIdentification for the ConsensusMap
      @param measure_from_subelements Do distance estimate from FeatureHandles instead of Centroid
      @param annotate_ids_with_subelements Store map index of FeatureHandle in peptide identification?
      @param spectra Whether precursors not contained in the identifications are annotated with 
                     an empty PeptideIdentification object containing the scan index.

      @exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'
    */
    void annotate(ConsensusMap& map, const std::vector<PeptideIdentification>& ids, 
                  const std::vector<ProteinIdentification>& protein_ids, 
                  bool measure_from_subelements = false, 
                  bool annotate_ids_with_subelements = false, 
                  const PeakMap& spectra = PeakMap());


    /**
      @brief Result of a partitioning by identification state with mapPrecursorsToIdentifications().
    */
    struct SpectraIdentificationState
    {
      std::vector<Size> no_precursors;
      std::vector<Size> identified;
      std::vector<Size> unidentified;
    }; 

    /**
      @brief Mapping of peptide identifications to spectra
             This helper function partitions all spectra into those that had: 
              - no precursor (e.g. MS1 spectra),
              - at least one identified precursor, 
              - or only unidentified precursor.
      @param spectra The mass spectra
      @param ids The peptide identifications
      @param mz_tol Tolerance used to map to precursor m/z
      @param rt_tol Tolerance used to map to spectrum retention time

      Note: mz/tol and rt_tol should, in principle, be zero (or close to zero under numeric inaccuracies). 

      @return A struct of vectors holding spectra indices of the partitioning.
    */
    static SpectraIdentificationState mapPrecursorsToIdentifications(const PeakMap& spectra, 
                                                                     const std::vector<PeptideIdentification>& ids, 
                                                                     double mz_tol = 0.001, 
                                                                     double rt_tol = 0.001)
    {
      SpectraIdentificationState ret;
      for (Size spectrum_index = 0; spectrum_index < spectra.size(); ++spectrum_index)
      {
        const MSSpectrum& spectrum = spectra[spectrum_index];
        if (!spectrum.getPrecursors().empty())
        {
          bool identified(false);
          const std::vector<Precursor>& precursors = spectrum.getPrecursors();

          // check if precursor has been identified
          for (Size i_p = 0; i_p < precursors.size(); ++i_p)
          {
            // check by precursor mass and spectrum RT
            double mz_p = precursors[i_p].getMZ();
            double rt_s = spectrum.getRT();
          
            for (Size i_id = 0; i_id != ids.size(); ++i_id)
            {
              const PeptideIdentification& pid = ids[i_id];

              // do not count empty ids as identification of a spectrum
              if (pid.getHits().empty()) continue;

              double mz_id = pid.getMZ();
              double rt_id = pid.getRT();

              if ( fabs(mz_id - mz_p) < mz_tol && fabs(rt_s - rt_id) < rt_tol )
              {
                identified = true;
                break; 
              }
            }
          }   
          if (!identified) 
          {
            ret.unidentified.push_back(spectrum_index);
          }
          else
          {
            ret.identified.push_back(spectrum_index);
          }
        }
        else
        {
          ret.no_precursors.push_back(spectrum_index);
        }
      }
      return ret;
    }


protected:
    void updateMembers_() override;

    /// Allowed RT deviation
    double rt_tolerance_;
    /// Allowed m/z deviation
    double mz_tolerance_;
    /// Measure used for m/z
    Measure measure_;
    /// Ignore charge states during matching?
    bool ignore_charge_;

    /// compute absolute Da tolerance, for a given m/z,
    /// when @p measure is MEASURE_DA, the value is unchanged,
    /// for MEASURE_PPM it is computed according to currently allowed ppm tolerance
    double getAbsoluteMZTolerance_(const double mz) const;

    /// check if distance constraint is fulfilled (using @p rt_tolerance_, @p mz_tolerance_ and @p measure_)
    bool isMatch_(const double rt_distance, const double mz_theoretical, const double mz_observed) const;

    /// helper function that checks if all peptide hits are annotated with RT and MZ meta values
    void checkHits_(const std::vector<PeptideIdentification>& ids) const;

    /// get RT, m/z and charge value(s) of a PeptideIdentification
    /// - multiple m/z values are returned if "mz_reference" is set to "peptide" (one for each PeptideHit)
    /// - one m/z value is returned if "mz_reference" is set to "precursor"
    void getIDDetails_(const PeptideIdentification& id, double& rt_pep, DoubleList& mz_values, IntList& charges, bool use_avg_mass = false) const;

    /// increase a bounding box by the given RT and m/z tolerances
    void increaseBoundingBox_(DBoundingBox<2>& box);

    /// try to determine the type of m/z value reported for features, return
    /// whether average peptide masses should be used for matching
    bool checkMassType_(const std::vector<DataProcessing>& processing) const;

  };

} // namespace OpenMS

