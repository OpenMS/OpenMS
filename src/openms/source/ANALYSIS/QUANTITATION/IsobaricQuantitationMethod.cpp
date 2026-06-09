// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTEighteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTThirtyTwoPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTThirtyFivePlexQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <algorithm>
#include <iterator> // std::size

namespace OpenMS
{
  namespace
  {
    // Thunk that instantiates a concrete method as a base-class pointer (stored as a function pointer in the registry).
    template <typename MethodT>
    std::unique_ptr<IsobaricQuantitationMethod> makeIsobaricMethod()
    {
      return std::make_unique<MethodT>();
    }

    using MT = IsobaricQuantitationMethod::MethodType;

    /// One row of the central isobaric-method registry.
    struct MethodRegistryEntry
    {
      MT type;                        ///< enum value; the registry is ordered by (and indexed with) it
      std::string_view name;          ///< canonical identifier (matches getMethodName()), e.g. "tmt6plex"
      std::string_view display_name;  ///< human-readable name, e.g. "TMT 6-plex"
      std::unique_ptr<IsobaricQuantitationMethod> (*factory)(); ///< nullptr for UNKNOWN (no instance)
    };

    // Single source of truth for every isobaric quantitation method (canonical name, display name, factory).
    // To add a method: add a MethodType enum value AND one row here -- the static_asserts below make this mandatory.
    constexpr MethodRegistryEntry METHOD_REGISTRY[] = {
      {MT::UNKNOWN,     "unknown",    "none",                 nullptr},
      {MT::TMT_6PLEX,   "tmt6plex",   "TMT 6-plex",           &makeIsobaricMethod<TMTSixPlexQuantitationMethod>},
      {MT::TMT_10PLEX,  "tmt10plex",  "TMT 10-plex",          &makeIsobaricMethod<TMTTenPlexQuantitationMethod>},
      {MT::TMT_11PLEX,  "tmt11plex",  "TMT 11-plex",          &makeIsobaricMethod<TMTElevenPlexQuantitationMethod>},
      {MT::TMT_16PLEX,  "tmt16plex",  "TMT 16-plex (TMTpro)", &makeIsobaricMethod<TMTSixteenPlexQuantitationMethod>},
      {MT::TMT_18PLEX,  "tmt18plex",  "TMT 18-plex",          &makeIsobaricMethod<TMTEighteenPlexQuantitationMethod>},
      {MT::TMT_32PLEX,  "tmt32plex",  "TMT 32-plex",          &makeIsobaricMethod<TMTThirtyTwoPlexQuantitationMethod>},
      {MT::TMT_35PLEX,  "tmt35plex",  "TMT 35-plex",          &makeIsobaricMethod<TMTThirtyFivePlexQuantitationMethod>},
      {MT::ITRAQ_4PLEX, "itraq4plex", "iTRAQ 4-plex",         &makeIsobaricMethod<ItraqFourPlexQuantitationMethod>},
      {MT::ITRAQ_8PLEX, "itraq8plex", "iTRAQ 8-plex",         &makeIsobaricMethod<ItraqEightPlexQuantitationMethod>},
    };

    // Enforce on every compiler that the registry stays in lockstep with the enum:
    // exactly one row per MethodType, in enum order (so we can index the registry by static_cast<size_t>(MethodType)).
    static_assert(std::size(METHOD_REGISTRY) == static_cast<size_t>(MT::SIZE_OF_METHODTYPE),
                  "METHOD_REGISTRY must contain exactly one row per MethodType. "
                  "Did you add a MethodType without registering it here?");

    constexpr bool methodRegistryIsOrdered()
    {
      for (size_t i = 0; i < std::size(METHOD_REGISTRY); ++i)
      {
        if (static_cast<size_t>(METHOD_REGISTRY[i].type) != i) { return false; }
      }
      return true;
    }
    static_assert(methodRegistryIsOrdered(),
                  "METHOD_REGISTRY rows must be in MethodType order (row i must have type i).");

    // Returns the registry row for @p mt. UNKNOWN is a valid in-range value (row 0); a value outside the
    // enum range (SIZE_OF_METHODTYPE or beyond) is rejected as a programming error.
    const MethodRegistryEntry& registryEntry(MT mt)
    {
      const auto idx = static_cast<size_t>(mt);
      if (idx >= std::size(METHOD_REGISTRY))
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Invalid IsobaricQuantitationMethod::MethodType value " + String(static_cast<int>(mt)) + ".");
      }
      return METHOD_REGISTRY[idx];
    }
  } // anonymous namespace

  std::string_view IsobaricQuantitationMethod::methodTypeName(MethodType mt)
  {
    return registryEntry(mt).name;
  }

  std::string_view IsobaricQuantitationMethod::methodDisplayName(MethodType mt)
  {
    return registryEntry(mt).display_name;
  }

  IsobaricQuantitationMethod::MethodType IsobaricQuantitationMethod::methodTypeFromName(std::string_view name)
  {
    for (const auto& entry : METHOD_REGISTRY)
    {
      if (entry.name == name) { return entry.type; }
    }
    return MethodType::UNKNOWN;
  }

  std::unique_ptr<IsobaricQuantitationMethod> IsobaricQuantitationMethod::create(MethodType mt)
  {
    // registryEntry() throws Exception::IllegalArgument for out-of-range / SIZE_OF_METHODTYPE
    const auto factory = registryEntry(mt).factory;
    return factory ? factory() : nullptr; // UNKNOWN has no factory -> disabled/none sentinel
  }

  IsobaricQuantitationMethod::~IsobaricQuantitationMethod() = default;

  IsobaricQuantitationMethod::IsobaricQuantitationMethod(MethodType method_type) :
    DefaultParamHandler("IsobaricQuantitationMethod"),
    iso_method_(method_type)
  {
  }

  const String& IsobaricQuantitationMethod::getMethodName() const
  {
    // one persistent String per method type (its canonical name), created once on first use,
    // so we can hand out a reference (methodTypeName() only yields a string_view)
    static const std::array<String, static_cast<size_t>(MethodType::SIZE_OF_METHODTYPE)> names = []
    {
      std::array<String, static_cast<size_t>(MethodType::SIZE_OF_METHODTYPE)> n;
      for (size_t i = 0; i < n.size(); ++i)
      {
        n[i] = String(methodTypeName(static_cast<MethodType>(i)));
      }
      return n;
    }();
    return names[static_cast<size_t>(iso_method_)];
  }

  Matrix<double> IsobaricQuantitationMethod::stringListToIsotopeCorrectionMatrix_(const StringList& stringlist) const
  {
    // check the string list
    if (stringlist.size() != getNumberOfChannels())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("IsobaricQuantitationMethod: Invalid string representation of the isotope correction matrix. Expected ") + getNumberOfChannels() + " entries but got " + stringlist.size() + ".");
    }

    // compute frequency matrix based on the deviation matrix
    Matrix<double> channel_frequency(getNumberOfChannels(), getNumberOfChannels(), 0.0);

    // channel index
    Size contributing_channel = 0;

    // fill row-wise
    for (const auto& l : stringlist)
    {
      StringList corrections;
      l.split('/', corrections);

      auto number_of_columns = getChannelInformation()[contributing_channel].affected_channels.size();
      if (corrections.size() != number_of_columns )
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Corrections for channel ID #") + contributing_channel + " must contain " + number_of_columns + " values, but has " + corrections.size() + "!", String(corrections.size()));
      }

      // overwrite line in Matrix with custom values
      Size affected_channel_idx = 0;
      double self_contribution = 100.0;
      double correction;
      Int target_channel;
      for (auto& c : corrections)
      {
        c = c.trim().toUpper();
        if (c != "NA" && c != "-1" && c != "0.0")
        {
          target_channel = getChannelInformation()[contributing_channel].affected_channels[affected_channel_idx];
          try
          {
            correction = c.toDouble();
          }
          catch (Exception::ConversionError& e)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Correction entry #") + affected_channel_idx + " in channel ID " + contributing_channel + " must be one of na/NA/-1 or a floating point number representation!", c);
          }

          if (correction < 0.0 || correction > 100.0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Correction entry #") + affected_channel_idx + " in channel ID " + contributing_channel + " must be a percentage between 0 and 100!", c);
          }
          
          if (target_channel >= 0 && Size(target_channel) < getNumberOfChannels())
          {
            channel_frequency(target_channel, contributing_channel) = correction / 100.0;
          }
          self_contribution -= correction; // count reduced self-contribution even if it does not affect another channel
        }
        affected_channel_idx++;
      }
      // set reduced self contribution
      channel_frequency(contributing_channel, contributing_channel) = self_contribution / 100.0;
      // increment channel index
      ++contributing_channel;
    }
    return channel_frequency;
  }

} // namespace
