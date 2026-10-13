// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt, Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/CONCEPT/HashUtils.h>

#include <functional>
#include <set>
#include <vector>

namespace OpenMS
{
  /**
   * @brief Class for storing MS run data with identifications of its spectra
   *
   * This class stores an MSExperiment (containing spectra) along with identification data (IdentificationData):
   * the identifications of the spectra, with their runs (search settings, databases, proteins). A spectrum links the
   * identifications of it (see getSpectrumIDQueries()), as a feature links identifications, and can link none or
   * several (e.g. of several searches).
   */
  class OPENMS_DLLAPI AnnotatedMSRun
  {
  public:
    using SpectrumType = MSExperiment::SpectrumType;
    using ChromatogramType = MSExperiment::ChromatogramType;


    /// Default constructor
    AnnotatedMSRun() = default;

    /**
     * @brief Move constructor for efficiently loading a MSExperiment without a deep copy
     * @param[in] experiment The MSExperiment to move into this object
     */
    explicit AnnotatedMSRun(MSExperiment&& experiment) : data(std::move(experiment))
    {};

    /// Move constructor
    AnnotatedMSRun(AnnotatedMSRun&&) = default;

    /// Copy constructor
    AnnotatedMSRun(const AnnotatedMSRun&) = default;
    AnnotatedMSRun& operator=(const AnnotatedMSRun&) = default;
    AnnotatedMSRun& operator=(AnnotatedMSRun&&) = default;

    /// Destructor
    ~AnnotatedMSRun() = default;

    /// Equality operator (spectra, identification data and the links of the spectra)
    bool operator==(const AnnotatedMSRun& rhs) const;

    /// Inequality operator
    bool operator!=(const AnnotatedMSRun& rhs) const
    {
      return !(*this == rhs);
    }

    /**
     * @brief Get the MSExperiment
     * @return A reference to the MSExperiment
     */
    MSExperiment& getMSExperiment();
    
    /**
     * @brief Get the MSExperiment (const version)
     * @return A const reference to the MSExperiment
     */
    const MSExperiment& getMSExperiment() const;

    /**
     * @brief Set the MSExperiment
     * @param[in] experiment The MSExperiment to set
     */
    void setMSExperiment(MSExperiment&& experiment);

    /**
     * @brief Set the MSExperiment
     * @param[in] experiment The MSExperiment to set
     */
    void setMSExperiment(const MSExperiment& experiment);

    /// @name Identifications as identification data
    ///@{
    /// The identifications of the spectra, with their runs (search settings, databases, proteins)
    IdentificationData& getIdentificationData()
    {
      return id_data_;
    }
    const IdentificationData& getIdentificationData() const
    {
      return id_data_;
    }

    /**
     * @brief The identifications that the spectrum at @p index links (e.g. those that IDMapper maps to it)
     *
     * A spectrum links identifications with all their matches. An identification can be linked by several spectra
     * (e.g. when it is mapped by retention time), and one that no spectrum links stays in the identification data.
     * Links follow the spectra by index: code that reorders or removes spectra keeps them in step.
     */
    const std::set<IdentificationData::QueryReference>& getSpectrumIDQueries(Size index) const;

    /// Link the spectrum at @p index to the identification @p query
    void addSpectrumIDQuery(Size index, const IdentificationData::QueryReference& query);

    /// Remove the identification links of all spectra
    void clearSpectrumIDQueries();

    /**
     * @brief The identifications that the spectrum at @p index links, each with all its matches, in the order of
     * their export (see IdentificationData::resolveLinks())
     *
     * @throw Exception::MissingInformation if a link refers to an identification that the identification data does not have
     */
    std::vector<IdentificationData::QueryMatches> getSpectrumIdentifications(Size index) const;

    /// Whether a spectrum links an identification
    bool hasSpectrumIdentifications() const;

    /// Remove the identification data and the links of all spectra
    void clearIdentifications();

    /**
     * @brief Erase the matches for which @p remove returns true; their identifications stay, also without matches
     *
     * @return The number of erased matches
     */
    Size eraseMatches(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&,
                                               const IdentificationData::Match&)>& remove);

    /**
     * @brief Erase the identifications for which @p remove returns true, with their matches and the links of the spectra to them
     *
     * @return The number of erased identifications
     */
    Size eraseIdentifications(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&)>& remove);
    ///@}

  private:

    MSExperiment data;
    IdentificationData id_data_;
    /// The identifications that each spectrum links, by spectrum index (spectra after the last one with links have none)
    std::vector<std::set<IdentificationData::QueryReference>> spectrum_queries_;
  };
} // namespace OpenMS

// Hash function specialization for AnnotatedMSRun
namespace std
{
  /**
   * @brief Hash function for OpenMS::AnnotatedMSRun.
   *
   * Hashes all fields used in operator==: data (MSExperiment), the identification data and the links of the spectra.
   * Since the component types don't have std::hash specializations, this hash uses
   * identifying properties (sizes and identifiers) to create a fast, consistent hash.
   */
  template<>
  struct hash<OpenMS::AnnotatedMSRun>
  {
    std::size_t operator()(const OpenMS::AnnotatedMSRun& run) const noexcept
    {
      // Start with hash of the MSExperiment size (number of spectra)
      std::size_t seed = OpenMS::hash_int(run.getMSExperiment().size());

      // Hash the number of chromatograms
      OpenMS::hash_combine(seed, OpenMS::hash_int(run.getMSExperiment().getChromatograms().size()));

      // Hash the identification runs and the links of the spectra
      OpenMS::hash_combine(seed, OpenMS::hash_int(run.getIdentificationData().getRuns().size()));
      for (OpenMS::Size index = 0; index < run.getMSExperiment().size(); ++index)
      {
        OpenMS::hash_combine(seed, OpenMS::hash_int(run.getSpectrumIDQueries(index).size()));
      }

      // Hash identifying properties from spectra
      for (const auto& spectrum : run.getMSExperiment().getSpectra())
      {
        OpenMS::hash_combine(seed, OpenMS::hash_float(spectrum.getRT()));
        OpenMS::hash_combine(seed, OpenMS::hash_int(spectrum.getMSLevel()));
        OpenMS::hash_combine(seed, OpenMS::hash_int(spectrum.size()));
      }

      return seed;
    }
  };
} // namespace std
