// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <cmath> // for modf() (which is an overloaded function in C++)
#include <vector>

namespace OpenMS
{

  namespace Math
  {

    /**
    @brief Provides access to linearly interpolated values (and
    derivatives) from discrete data points.  Values beyond the given range
    of data points are implicitly taken as zero.

    The input is just a vector of values ("Data").  These are interpreted
    as the y-coordinates at the x-coordinate positions 0,...,data_.size-1.

    The interpolated data can also be <i>scaled</i> and <i>shifted</i> in
    the x-dimension by an <em>affine mapping</em>.  That is, we have "inside" and
    "outside" x-coordinates.  The affine mapping can be specified in two
    ways:
    - using setScale() and setOffset(),
    - using setMapping()
    .
    By default the identity mapping (scale=1, offset=0) is used.

    Using the value() and derivative() methods you can sample linearly
    interpolated values for a given x-coordinate position of the data and
    the derivative of the data.

    @see BilinearInterpolation

    @ingroup Math
    */
    template <typename Key = double, typename Value = Key>
    class LinearInterpolation
    {

public:

      ///\name Typedefs
      //@{
      typedef Value value_type;

      typedef Key key_type;
      typedef std::vector<value_type> container_type;

      typedef value_type      ValueType;
      typedef key_type        KeyType;
      typedef container_type  ContainerType;
      //@}

public:

      /**@brief Constructors and destructor.

      The first argument is the scale which is applied to the arguments of
      value() and derivative() before looking up the interpolated values in
      the container.  The second argument is the offset, which is
      subtracted before everything else.
      */
      LinearInterpolation(KeyType scale = 1., KeyType offset = 0.) :
        scale_(scale),
        offset_(offset),
        inside_(),
        outside_(),
        data_()
      {}

      /// Copy constructor.
      LinearInterpolation(LinearInterpolation const & arg) :
        scale_(arg.scale_),
        offset_(arg.offset_),
        inside_(arg.inside_),
        outside_(arg.outside_),
        data_(arg.data_)
      {}

      /// Assignment operator
      LinearInterpolation & operator=(LinearInterpolation const & arg)
      {
        if (&arg == this)
          return *this;

        scale_   = arg.scale_;
        offset_  = arg.offset_;
        inside_  = arg.inside_;
        outside_ = arg.outside_;
        data_    = arg.data_;

        return *this;
      }

      /// Destructor.
      ~LinearInterpolation() = default;

      // ----------------------------------------------------------------------

      ///@name Interpolated data
      //@{

      /// Returns the interpolated value.
      ValueType value(KeyType arg_pos) const
      {

        typedef typename container_type::difference_type DiffType;

        // apply the key transformation
        KeyType left_key;
        KeyType pos = key2index(arg_pos);
        KeyType frac = std::modf(pos, &left_key);
        DiffType const left = DiffType(left_key);

        // At left margin?
        if (pos < 0)
        {
          if (left /* <= -1 */)
          {
            return 0;
          }
          else        // left == 0
          {
            return data_[0] * (1 + frac);
          }
        }
        else         // pos >= 0
        {
          // At right margin?
          DiffType const back = data_.size() - 1;
          if (left >= back)
          {
            if (left != back)
            {
              return 0;
            }
            else
            {
              return data_[left] * (1 - frac);
            }
          }
          else
          {
            // In between!
            return data_[left + 1] * frac + data_[left] * (1 - frac);
          }
        }
      }

      /** @brief Performs linear resampling.  The arg_value is split up and
      added to the data points around arg_pos.
      */
      void addValue(KeyType arg_pos, ValueType arg_value)
      {

        typedef typename container_type::difference_type DiffType;

        // apply the key transformation
        KeyType left_key;
        KeyType const pos = key2index(arg_pos);
        KeyType const frac = std::modf(pos, &left_key);
        DiffType const left = DiffType(left_key);

        // At left margin?
        if (pos < 0)
        {
          if (left /* <= -1 */)
          {
            return;
          }
          else        // left == 0
          {
            data_[0] += (1 + frac) * arg_value;
            return;
          }
        }
        else         // pos >= 0
        {
          // At right margin?
          DiffType const back = data_.size() - 1;
          if (left >= back)
          {
            if (left != back)
            {
              return;
            }
            else             // left == back
            {
              data_[left] += (1 - frac) * arg_value;
              return;
            }
          }
          else
          {
            // In between!
            data_[left + 1] += frac * arg_value;
            data_[left] += (1 - frac) * arg_value;
            return;
          }
        }
      }

      /** @brief Returns the interpolated derivative.

      Please drop me (= the maintainer) a message if you are using this.
      */
      ValueType derivative(KeyType arg_pos) const
      {

        // apply the key transformation
        KeyType const pos = key2index(arg_pos);

        SignedSize const size_ = data_.size();
        SignedSize const left = int(pos + 0.5);       // rounds towards zero

        if (left < 0)           // quite small
        {
          return 0;
        }
        else
        {
          if (left == 0)             // at the border
          {
            if (pos >= -0.5)               // that is: -0.5 <= pos < +0.5
            {
              return (data_[1] - data_[0]) * (pos + 0.5) + (data_[0]) * (0.5 - pos);
            }
            else             // that is: -1.5 <= pos < -0.5
            {
              return (data_[0]) * (pos + 1.5);
            }
          }
        }
        // "else" case: to the right of the left margin


        KeyType factor = KeyType(left) - pos + KeyType(0.5);

        if (left > size_)           // quite large
        {
          return 0;
        }
        else
        {
          if (left < size_ - 1)             // to the left of the right margin
          {
            // weighted average of derivatives for adjacent intervals
            return (data_[left] - data_[left - 1]) * factor + (data_[left + 1] - data_[left]) * (1. - factor);
          }
          else           // somewhat at the border
          {
            // at the border, first case
            if (left == size_ - 1)
            {
              return (data_[left] - data_[left - 1]) * factor + (-data_[left]) * (1. - factor);
            }
          }
        }
        // else // that is: left == size_

        // We pull the last remaining case out of the "if" tree to avoid a
        // compiler warning ...

        // at the border, second case
        return (-data_[left - 1]) * factor;
      }

      //@}

      // ----------------------------------------------------------------------

      ///@name Discrete (non-interpolated) data
      //@{

      /// Returns the internal random access container from which interpolated values are being sampled.
      ContainerType & getData()
      {
        return data_;
      }

      /// Returns the internal random access container from which interpolated values are being sampled.
      ContainerType const & getData() const
      {
        return data_;
      }

      /** @brief Assigns data to the internal random access container from
      which interpolated values are being sampled.

      SourceContainer must be assignable to ContainerType.
      */
      template <typename SourceContainer>
      void setData(SourceContainer const & data)
      {
        data_ = data;
      }

      /// Returns \c true if getData() is empty.
      bool empty() const
      {
        return data_.empty();
      }

      //@}

      // ----------------------------------------------------------------------

      ///\name Transformation
      //@{

      /// The transformation from "outside" to "inside" coordinates.
      KeyType key2index(KeyType pos) const
      {
        if (scale_)
        {
          pos -= offset_;
          pos /= scale_;
          return pos;
        }
        else
        {
          return 0;
        }
      }

      /// The transformation from "inside" to "outside" coordinates.
      KeyType index2key(KeyType pos) const
      {
        pos *= scale_;
        pos += offset_;
        return pos;
      }

      /// Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".
      KeyType const & getScale() const
      {
        return scale_;
      }

      /** @brief Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".

      <b>Note:</b> Using this invalidates the inside and outside reference
      points.
      */
      void setScale(KeyType const & scale)
      {
        scale_ = scale;
      }

      /// Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data[0]".
      KeyType const & getOffset() const
      {
        return offset_;
      }

      /**@brief Accessor.  "Offset" is the point (in "outside" units) which
      corresponds to "Data[0]".

      <b>Note:</b> Using this invalidates the inside and outside reference
      points.
      */
      void setOffset(KeyType const & offset)
      {
        offset_ = offset;
      }

      /** @brief Specifies the mapping from "outside" to "inside" coordinates by the following data:
       - <code>scale</code>: the difference in outside coordinates between consecutive values in the data vector.
       - <code>inside</code> and <code>outside</code>: these x-axis positions are mapped onto each other.

      For example, when you have a complicated probability distribution
      which is in fact centered around zero (but you cannot have negative
      indices in the data vector), then you can arrange things such that
      inside is the mean of the pre-computed, shifted density values of that
      distribution and outside is the centroid position of, say, a peak in
      the real world which you want to model by a scaled and shifted version
      of the probability distribution.

      */
      void setMapping(KeyType const & scale, KeyType const & inside, KeyType const & outside)
      {
        scale_   = scale;
        inside_  = inside;
        outside_ = outside;
        offset_  = outside - scale * inside;
      }

      /** @brief Specifies the mapping from "outside" to "inside" coordinates by the following data:
      - <code>inside_low</code> and <code>outside_low</code>: these axis positions are mapped onto each other.
      - <code>inside_high</code> and <code>outside_high</code>: these axis positions are mapped onto each other.

      This four argument version is just a convenience overload for the three argument version, which see.
      */
      void setMapping(KeyType const & inside_low, KeyType const & outside_low,
                      KeyType const & inside_high, KeyType const & outside_high)
      {
        if (inside_high != inside_low)
        {
          setMapping((outside_high - outside_low) / (inside_high - inside_low),
                     inside_low, outside_low);
        }
        else
        {
          setMapping(0, inside_low, outside_low);
        }
        return;
      }

      /// Accessor.  See setMapping().
      KeyType const & getInsideReferencePoint() const
      {
        return inside_;
      }

      /// Accessor.  See setMapping().
      KeyType const & getOutsideReferencePoint() const
      {
        return outside_;
      }

      /// Lower boundary of the support, in "outside" coordinates.
      KeyType supportMin() const
      {
        return index2key(KeyType(empty() ? 0 : -1));
      }

      /// Upper boundary of the support, in "outside" coordinates.
      KeyType supportMax() const
      {
        return index2key(KeyType(data_.size()));
      }

      //@}

protected:

      KeyType scale_;
      KeyType offset_;
      KeyType inside_;
      KeyType outside_;

      ContainerType data_;

    };

  }   // namespace Math

} // namespace OpenMS

