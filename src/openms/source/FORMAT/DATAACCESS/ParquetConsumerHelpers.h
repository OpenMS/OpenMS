// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <arrow/api.h>

#include <charconv>
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
        @brief Assigns the numeric PRECURSOR_ID / TRANSITION_ID values of a writer from the string
        IDs of a transition list.

        A string that is a non-negative integer keeps that value; anything else (non-numeric
        strings, negative numbers, a value already taken by an earlier entry) gets a fresh ID.
        Negative values are "unset" sentinels elsewhere (the readers take a negative ID as "no
        filter"); the first one this assigner replaces is logged as a warning naming @p context and
        @p column. OpenSWATH's own library loading validates IDs as non-negative, so only callers
        that pass a non-canonical LightTargetedExperiment see replacements at all.

        Call reserve() for every string ID before the first assign(): fresh IDs then never take a
        value that an entry later in the list carries explicitly, so explicit IDs are kept even when
        replaced or non-numeric IDs come first.

        Throws Exception::InvalidValue if no fresh ID is left below INT64_MAX (it never wraps to a
        negative value).
      */
      class IdAssigner
      {
      public:
        IdAssigner(const char* context, const char* column, int64_t first_fresh_id = 1) :
          context_(context), column_(column), next_id_(first_fresh_id)
        {
        }

        /// Claim the explicit ID in @p text (if it is a non-negative integer) for its first assign().
        void reserve(const std::string& text)
        {
          int64_t value = 0;
          if (parse_(text, value) && value >= 0)
          {
            reserved_.insert(value);
            if (value >= next_id_ && value < std::numeric_limits<int64_t>::max())
            {
              next_id_ = value + 1; // fresh IDs start above the explicit ones
            }
          }
        }

        /// The ID for @p text: its own value if that is a non-negative integer not assigned yet,
        /// otherwise the next free ID that is neither assigned nor reserved.
        int64_t assign(const std::string& text)
        {
          int64_t value = 0;
          if (parse_(text, value))
          {
            if (value < 0)
            {
              if (!warned_)
              {
                warned_ = true;
                OPENMS_LOG_WARN << context_ << ": negative " << column_ << " '" << text
                                << "' in the transition list replaced by a newly assigned ID (negative IDs mean"
                                   " 'unset' in OpenSWATH Parquet files). Further negative IDs in this file are"
                                   " replaced without a warning." << std::endl;
              }
            }
            else if (assigned_.insert(value).second)
            {
              return value;
            }
          }
          while (assigned_.contains(next_id_) || reserved_.contains(next_id_))
          {
            advance_();
          }
          const int64_t id = next_id_;
          assigned_.insert(id);
          if (id < std::numeric_limits<int64_t>::max())
          {
            ++next_id_; // at INT64_MAX, the next fresh request fails in advance_() instead
          }
          return id;
        }

      private:
        void advance_()
        {
          if (next_id_ == std::numeric_limits<int64_t>::max())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          std::string(context_) + ": no free " + column_ + " left below INT64_MAX",
                                          StringUtils::toStr(next_id_));
          }
          ++next_id_;
        }

        /// Same accepted forms as StringUtils::toInt64() (surrounding whitespace, a leading '+'), but
        /// without throwing for the non-numeric IDs that are normal here.
        static bool parse_(const std::string& text, int64_t& value)
        {
          const auto is_ws = [](char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; };
          const char* first = text.data();
          const char* last = first + text.size();
          while (first != last && is_ws(*first)) ++first;
          while (last != first && is_ws(*(last - 1))) --last;
          if (first != last && *first == '+') ++first;
          if (first == last) return false;
          const auto [ptr, ec] = std::from_chars(first, last, value);
          return ec == std::errc() && ptr == last;
        }

        const char* context_;
        const char* column_;
        bool warned_ = false;
        int64_t next_id_;
        std::unordered_set<int64_t> reserved_;
        std::unordered_set<int64_t> assigned_;
      };

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
