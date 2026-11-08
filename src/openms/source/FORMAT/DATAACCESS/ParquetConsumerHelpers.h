// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <arrow/api.h>

#include <cstdint>
#include <limits>
#include <string>
#include <unordered_set>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Private helpers shared by the OpenSWATH Parquet writers
      (MSChromatogramParquetConsumer, MobilogramParquetConsumer, XIPMParquetConsumer).

      Not installed. The append helpers take the column name as a C string: they run once per
      written value, and ParquetFile::appendOrThrow() would build a std::string for every call.
    */
    namespace ParquetConsumerHelpers
    {
      /// Throw Exception::InvalidValue if appending to @p column failed.
      inline void appendOrThrow(const arrow::Status& status, const char* column)
      {
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string("Failed to append value for ") + column, status.ToString());
        }
      }

      /// Throw Exception::InvalidValue if reserving capacity for @p column failed.
      inline void reserveOrThrow(const arrow::Status& status, const char* column)
      {
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string("Failed to reserve capacity for ") + column, status.ToString());
        }
      }

      /// Append @p value, or null if @p has_value is false.
      inline void appendOptionalInt(arrow::Int64Builder& builder, bool has_value, int64_t value, const char* column)
      {
        if (!has_value)
        {
          appendOrThrow(builder.AppendNull(), column);
        }
        else
        {
          appendOrThrow(builder.Append(value), column);
        }
      }

      /// Append @p value, or null if it is empty.
      inline void appendOptionalString(arrow::StringBuilder& builder, const std::string& value, const char* column)
      {
        if (value.empty())
        {
          appendOrThrow(builder.AppendNull(), column);
        }
        else
        {
          appendOrThrow(builder.Append(value), column);
        }
      }

      /// Append @p value as one binary cell; throws if it exceeds Arrow's 2 GB cell limit.
      inline void appendBinary(arrow::BinaryBuilder& builder, const std::string& value, const char* column)
      {
        if (value.size() > static_cast<size_t>(std::numeric_limits<int32_t>::max()))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string("Binary value for column ") + column + " exceeds Arrow int32 size limit (2 GB).",
                                        StringUtils::toStr(value.size()));
        }
        appendOrThrow(builder.Append(reinterpret_cast<const uint8_t*>(value.data()),
                                     static_cast<int32_t>(value.size())), column);
      }

      /// Precursor metadata of a compound of the transition list.
      struct CompoundInfo
      {
        int64_t precursor_id = 0;
        std::string modified_sequence;
        int64_t precursor_charge = 0;
        int64_t precursor_decoy = 0;
      };

      /// Metadata of a transition of the transition list.
      struct TransitionInfo
      {
        int64_t transition_id = 0;
        int64_t precursor_id = 0;
        std::string modified_sequence;
        int64_t precursor_charge = 0;
        int64_t product_charge = 0;
        int64_t detecting_transition = 0;
        int64_t precursor_decoy = 0;
        int64_t product_decoy = 0;
        int64_t transition_ordinal = 0;
        std::string transition_type;
        std::string annotation;
      };

      /**
        @brief The numeric ID for @p text: its value if it is a non-negative integer not used yet,
        otherwise the next free ID (from @p next_id).

        Negative values are "unset" sentinels elsewhere and are never kept as IDs.
      */
      inline int64_t parseOrAssignId(const std::string& text, int64_t& next_id, std::unordered_set<int64_t>& used_ids)
      {
        try
        {
          const int64_t value = StringUtils::toInt64(text);
          if (value >= 0 && used_ids.insert(value).second)
          {
            if (value >= next_id)
            {
              next_id = value + 1;
            }
            return value;
          }
        }
        catch (Exception::ConversionError&)
        {
          // fall through to auto-assigned ID
        }
        while (used_ids.contains(next_id))
        {
          ++next_id;
        }
        used_ids.insert(next_id);
        return next_id++;
      }

      /// Fragment annotation such as "y7^2" (empty if the type is unknown or the ordinal negative).
      inline std::string buildAnnotation(const std::string& transition_type, int64_t ordinal, int64_t charge)
      {
        if (transition_type.empty() || ordinal < 0)
        {
          return "";
        }
        std::string annotation = transition_type + StringUtils::toStr(ordinal);
        if (charge > 0)
        {
          annotation += "^" + StringUtils::toStr(charge);
        }
        return annotation;
      }

      /// The "Precursor_i..." suffix of a precursor trace's native ID (empty if it has none).
      inline std::string buildPrecursorAnnotation(const std::string& native_id)
      {
        const std::string tag = "_Precursor_i";
        const Size pos = native_id.rfind(tag);
        if (pos != std::string::npos)
        {
          return native_id.substr(pos + 1);
        }
        return "";
      }
    } // namespace ParquetConsumerHelpers
  } // namespace Internal
} // namespace OpenMS
