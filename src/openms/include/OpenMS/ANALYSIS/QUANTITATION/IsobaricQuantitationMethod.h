// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/KERNEL/Peak2D.h>

#include <array>
#include <memory>
#include <string_view>
#include <utility>
#include <vector>

namespace OpenMS
{

  // Forward declarations
  template <typename Value>
  class Matrix;

  /**
    @brief Abstract base class describing an isobaric quantitation method in terms of the used channels and an isotope correction matrix.
  */
  class OPENMS_DLLAPI IsobaricQuantitationMethod :
    public DefaultParamHandler
  {
public:

    /// Identifies a concrete isobaric quantitation method. UNKNOWN is used as a sentinel (disabled/none).
    enum class MethodType
    {
      UNKNOWN = 0,
      TMT_6PLEX,
      TMT_10PLEX,
      TMT_11PLEX,
      TMT_16PLEX,
      TMT_18PLEX,
      TMT_32PLEX,
      TMT_35PLEX,
      ITRAQ_4PLEX,
      ITRAQ_8PLEX,
      SIZE_OF_METHODTYPE ///< keep last
    };

    /// Returns the canonical string identifier for a given MethodType (e.g. "tmt6plex"); "unknown" for MethodType::UNKNOWN.
    /// @throws Exception::IllegalArgument for an out-of-range value (MethodType::SIZE_OF_METHODTYPE or beyond).
    static std::string_view methodTypeName(MethodType mt);

    /// Returns a human-readable display name for a given MethodType (e.g. "TMT 6-plex"); "none" for MethodType::UNKNOWN.
    /// @throws Exception::IllegalArgument for an out-of-range value (MethodType::SIZE_OF_METHODTYPE or beyond).
    static std::string_view methodDisplayName(MethodType mt);

    /// Returns the MethodType corresponding to @p name (canonical identifier as returned by methodTypeName()), or MethodType::UNKNOWN if not found.
    static MethodType methodTypeFromName(std::string_view name);

    /**
      @brief Factory: creates a new instance of the concrete quantitation method identified by @p mt.

      @param mt The quantitation method to instantiate.
      @return An owning pointer to the newly created method, or @c nullptr if @p mt is MethodType::UNKNOWN (the disabled/none sentinel).
      @throws Exception::IllegalArgument if @p mt is not a concrete method, i.e. MethodType::SIZE_OF_METHODTYPE or any value outside the enum range.
    */
    static std::unique_ptr<IsobaricQuantitationMethod> create(MethodType mt);

    /**
      @brief Summary of an isobaric quantitation channel.
    */
    struct IsobaricChannelInformation
    {
      /// The name of the channel.
      std::string name;
      /// The id of the channel.
      Int id;
      /// Optional description of the channel.
      std::string description;
      /// The expected centroid position of the channel peak in m/z.
      Peak2D::CoordinateType center;
      /// Ids of the affected channels. Must contain 4 or 8 entries, depending on the number of columns in the Thermo data sheet (with or without subchannels). Order has to match the ones from the correction matrix parameter.
      std::vector<Int> affected_channels;

      /// C'tor
      IsobaricChannelInformation(std::string local_name,
                                 const Int local_id,
                                 std::string local_description,
                                 const Peak2D::CoordinateType& local_center,
                                 const std::vector<Int>& affected_channels) :
          name(std::move(local_name)),
          id(local_id),
          description(std::move(local_description)),
          center(local_center),
          affected_channels(affected_channels)
      {
      }
    };

    /// @brief Default c'tor is deleted: a concrete method must be constructed with its MethodType (see the protected c'tor).
    IsobaricQuantitationMethod() = delete;

    /// @brief d'tor
    ~IsobaricQuantitationMethod() override;

    typedef std::vector<IsobaricChannelInformation> IsobaricChannelList;

    /**
      @brief Returns a unique name for the quantitation method (its canonical identifier, e.g. "tmt6plex").

      @return The unique name or identifier of the quantitation method.
    */
    const String& getMethodName() const;

    /// Returns the MethodType enum value of this quantitation method.
    MethodType getMethodType() const { return iso_method_; }

    /**
      @brief Returns information on the different channels used by the quantitation method.

      @return A std::vector containing the channel information for this quantitation method.
    */
    virtual const IsobaricChannelList& getChannelInformation() const = 0;

    /**
      @brief Gives the number of channels available for this quantitation method.

      @return The number of channels available for this quantitation method.
    */
    virtual Size getNumberOfChannels() const = 0;

    /**
      @brief Returns an isotope correction matrix suitable for the given quantitation method.
    */
    virtual Matrix<double> getIsotopeCorrectionMatrix() const = 0;

    /**
      @brief Returns the index of the reference channel in the IsobaricChannelList (see IsobaricQuantitationMethod::getChannelInformation()).
    */
    virtual Size getReferenceChannel() const = 0;

protected:
    /// @brief c'tor for derived classes: sets the underlying param-handler name and records the concrete @p method_type.
    explicit IsobaricQuantitationMethod(MethodType method_type);

    /**
      @brief Helper function to convert a string list containing an isotope correction matrix into a Matrix<double>.

      @param[in] stringlist The StringList to convert.
      @return An isotope correction matrix as Matrix<double>.
    */
    Matrix<double> stringListToIsotopeCorrectionMatrix_(const std::vector<String>& stringlist) const;

private:
    /// The concrete isobaric method this instance represents; set by the c'tor and returned by getMethodType()/getMethodName().
    MethodType iso_method_;
  };
} // namespace

