// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/KERNEL/Peak2D.h>

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

    /**
      @brief Summary of an isobaric quantitation channel.
    */
    struct IsobaricChannelInformation
    {
      /// The name of the channel.
      String name;
      /// The id of the channel.
      Int id;
      /// Optional description of the channel.
      String description;
      /// The expected centroid position of the channel peak in m/z.
      Peak2D::CoordinateType center;
      /// Ids of the affected channels. Must contain 4 or 8 entries, depending on the number of columns in the Thermo data sheet (with or without subchannels). Order has to match the ones from the correction matrix parameter.
      std::vector<Int> affected_channels;

      /// C'tor
      IsobaricChannelInformation(String local_name,
                                 const Int local_id,
                                 String  local_description,
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

    /// @brief c'tor setting the name for the underlying param handler
    IsobaricQuantitationMethod();

    /// @brief d'tor
    ~IsobaricQuantitationMethod() override;

    typedef std::vector<IsobaricChannelInformation> IsobaricChannelList;

    /**
      @brief Returns a unique name for the quantitation method.

      @return The unique name or identifier of the quantitation method.
    */
    virtual const String& getMethodName() const = 0;

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
    /**
      @brief Helper function to convert a string list containing an isotope correction matrix into a Matrix<double>.

      @param stringlist The StringList to convert.
      @return An isotope correction matrix as Matrix<double>.
    */
    Matrix<double> stringListToIsotopeCorrectionMatrix_(const std::vector<String>& stringlist) const;
  };
} // namespace

