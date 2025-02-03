// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Matrix.h>

namespace OpenMS
{

  namespace Math
  {

    /**
         @brief Provides access to bilinearly interpolated values (and
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

         Using the value() and derivative() methods you can sample bilinearly
         interpolated values for a given x-coordinate position of the data and
         the derivative of the data.

         @see LinearInterpolation


         @ingroup Math
    */
    template <typename Key = double, typename Value = Key>
    class BilinearInterpolation
    {

public:

      ///@name Typedefs
      //@{
      typedef Value value_type;

      typedef Key key_type;
      typedef Matrix<value_type> container_type;

      typedef value_type ValueType;
      typedef key_type KeyType;
      typedef container_type ContainerType;
      //@}

public:

      /**@brief Constructors and destructor.
      */
      //@{

      /// Default constructor
      BilinearInterpolation() :
        scale_0_(1),
        offset_0_(0),
        scale_1_(1),
        offset_1_(0),
        inside_0_(0),
        outside_0_(0),
        inside_1_(0),
        outside_1_(0),
        data_()
      {}

      /// Copy constructor
      BilinearInterpolation(BilinearInterpolation const & arg) :
        scale_0_(arg.scale_0_),
        offset_0_(arg.offset_0_),
        scale_1_(arg.scale_1_),
        offset_1_(arg.offset_1_),
        inside_0_(arg.inside_0_),
        outside_0_(arg.outside_0_),
        inside_1_(arg.inside_1_),
        outside_1_(arg.outside_1_),
        data_(arg.data_)
      {}

      /// Assignment operator
      BilinearInterpolation & operator=(BilinearInterpolation const & arg)
      {
        if (&arg == this)
          return *this;

        scale_0_ = arg.scale_0_;
        offset_0_ = arg.offset_0_;
        scale_1_ = arg.scale_1_;
        offset_1_ = arg.offset_1_;
        inside_0_ = arg.inside_0_;
        outside_1_ = arg.outside_1_;
        inside_1_ = arg.inside_1_;
        outside_0_ = arg.outside_0_;
        data_ = arg.data_;
        return *this;
      }

      /// Destructor
      ~BilinearInterpolation()
      {}

      //@}

      // ----------------------------------------------------------------------

      ///@name Interpolated data
      //@{

      /// Returns the interpolated value ("backward resampling")
      ValueType value(KeyType arg_pos_0, KeyType arg_pos_1) const
      {
        // apply the key transformations
        KeyType const pos_0 = key2index_0(arg_pos_0);
        KeyType const pos_1 = key2index_1(arg_pos_1);

        // ???? should use modf() here!

        SignedSize const size_0 = data_.rows();
        SignedSize const lower_0 = SignedSize(pos_0);           // this rounds towards zero
        SignedSize const size_1 = data_.cols();
        SignedSize const lower_1 = SignedSize(pos_1);           // this rounds towards zero

        // small pos_0
        if (pos_0 <= 0)
        {
          if (lower_0 != 0)
          {
            return 0;
          }
          else        // that is: -1 < pos_0 <= 0
          {             // small pos_1
            if (pos_1 <= 0)
            {
              if (lower_1 != 0)
              {
                return 0;
              }
              else            // that is: -1 < pos_1 <= 0
              {
                return data_(0, 0) * (1. + pos_0) * (1. + pos_1);
              }
            }

            // big pos_1
            if (lower_1 >= size_1 - 1)
            {
              if (lower_1 != size_1 - 1)
              {
                return 0;
              }
              else
              {
                return data_(0, lower_1) * (1. + pos_0) * (size_1 - pos_1);
              }
            }

            // medium pos_1
            KeyType const factor_1 = pos_1 - KeyType(lower_1);
            KeyType const factor_1_complement = KeyType(1.) - factor_1;
            return (
                     data_(0, lower_1 + 1) * factor_1 +
                     data_(0, lower_1) * factor_1_complement
                     ) * (1. + pos_0);
          }
        }

        // big pos_0
        if (lower_0 >= size_0 - 1)
        {
          if (lower_0 != size_0 - 1)
          {
            return 0;
          }
          else        // that is: size_0 - 1 <= pos_0 < size_0
          {             // small pos_1
            if (pos_1 <= 0)
            {
              if (lower_1 != 0)
              {
                return 0;
              }
              else            // that is: -1 < pos_1 <= 0
              {
                return data_(lower_0, 0) * (size_0 - pos_0) * (1. + pos_1);
              }
            }

            // big pos_1
            if (lower_1 >= size_1 - 1)
            {
              if (lower_1 != size_1 - 1)
              {
                return 0;
              }
              else
              {
                return data_(lower_0, lower_1) * (size_0 - pos_0) * (size_1 - pos_1);
              }
            }

            // medium pos_1
            KeyType const factor_1 = pos_1 - KeyType(lower_1);
            KeyType const factor_1_complement = KeyType(1.) - factor_1;
            return (
                     data_(lower_0, lower_1 + 1) * factor_1 +
                     data_(lower_0, lower_1) * factor_1_complement
                     )
                   * (size_0 - pos_0);
          }
        }

        // medium pos_0
        {
          KeyType const factor_0 = pos_0 - KeyType(lower_0);
          KeyType const factor_0_complement = KeyType(1.) - factor_0;

          // small pos_1
          if (pos_1 <= 0)
          {
            if (lower_1 != 0)
            {
              return 0;
            }
            else          // that is: -1 < pos_1 <= 0
            {
              return (
                       data_(lower_0 + 1, 0) * factor_0
                       +
                       data_(lower_0, 0) * factor_0_complement
                       )
                     * (1. + pos_1);
            }
          }

          // big pos_1
          if (lower_1 >= size_1 - 1)
          {
            if (lower_1 != size_1 - 1)
            {
              return 0;
            }
            else
            {
              return (
                       data_(lower_0 + 1, lower_1) * factor_0
                       +
                       data_(lower_0, lower_1) * factor_0_complement
                       )
                     * (size_1 - pos_1);
            }
          }
          KeyType const factor_1 = pos_1 - KeyType(lower_1);
          KeyType const factor_1_complement = KeyType(1.) - factor_1;

          // medium pos_0 and medium pos_1 --> "within" the matrix
          return (
                   data_(lower_0 + 1, lower_1 + 1) * factor_0
                   +
                   data_(lower_0, lower_1 + 1) * factor_0_complement
                   )
                 * factor_1
                 +
                 (
                   data_(lower_0 + 1, lower_1) * factor_0
                   +
                   data_(lower_0, lower_1) * factor_0_complement
                 )
                 * factor_1_complement;
        }
      }

      /**@brief Performs bilinear resampling.  The arg_value is split up and
           added to the data points around arg_pos.  ("forward resampling")
      */
      void addValue(KeyType arg_pos_0, KeyType arg_pos_1, ValueType arg_value)
      {

        typedef typename std::ptrdiff_t DiffType;

        // apply key transformation _0
        KeyType const pos_0 = key2index_0(arg_pos_0);
        KeyType lower_0_key;
        KeyType const frac_0 = std::modf(pos_0, &lower_0_key);
        DiffType const lower_0 = DiffType(lower_0_key);

        // Small pos_0 ?
        if (pos_0 < 0)
        {
          if (lower_0)
          {
            return;
          }
          else        // lower_0 == 0
          {             // apply key transformation _1
            KeyType const pos_1 = key2index_1(arg_pos_1);
            KeyType lower_1_key;
            KeyType const frac_1 = std::modf(pos_1, &lower_1_key);
            DiffType const lower_1 = DiffType(lower_1_key);

            // Small pos_1 ?
            if (pos_1 < 0)
            {
              if (lower_1)
              {
                return;
              }
              else            // lower_1 == 0
              {
                data_(0, 0) += arg_value * (1 + frac_0) * (1 + frac_1);
                return;
              }
            }
            else             // pos_1 >= 0
            {
              DiffType const back_1 = data_.cols() - 1;
              // big pos_1
              if (lower_1 >= back_1)
              {
                if (lower_1 != back_1)
                {
                  return;
                }
                else                 // lower_1 == back_1
                {
                  data_(0, lower_1) += arg_value * (1 + frac_0) * (1 - frac_1);
                  return;
                }
              }
              else
              {
                // medium pos_1
                KeyType const tmp_prod = KeyType(arg_value * (1. + frac_0));
                data_(0, lower_1 + 1) += tmp_prod * frac_1;
                data_(0, lower_1) += tmp_prod * (1. - frac_1);
                return;
              }
            }
          }
        }
        else         // pos_0 >= 0
        {
          DiffType const back_0 = data_.rows() - 1;
          if (lower_0 >= back_0)
          {
            if (lower_0 != back_0)
            {
              return;
            }
            else             // lower_0 == back_0
            {

              KeyType const tmp_prod = KeyType(arg_value * (1. - frac_0));

              // apply key transformation _1
              KeyType const pos_1 = key2index_1(arg_pos_1);
              KeyType lower_1_key;
              KeyType const frac_1 = std::modf(pos_1, &lower_1_key);
              DiffType const lower_1 = DiffType(lower_1_key);

              // Small pos_1 ?
              if (pos_1 < 0)
              {
                if (lower_1)
                {
                  return;
                }
                else              // lower_1 == 0
                {
                  data_(lower_0, 0) += tmp_prod * (1 + frac_1);
                  return;
                }
              }
              else               // pos_1 >= 0
              {
                DiffType const back_1 = data_.cols() - 1;
                // big pos_1
                if (lower_1 >= back_1)
                {
                  if (lower_1 != back_1)
                  {
                    return;
                  }
                  else                   // lower_1 == back_1
                  {
                    data_(lower_0, lower_1) += tmp_prod * (1 - frac_1);
                    return;
                  }
                }
                else
                {
                  // medium pos_1
                  data_(lower_0, lower_1 + 1) += tmp_prod * frac_1;
                  data_(lower_0, lower_1) += tmp_prod * (1 - frac_1);
                  return;
                }
              }
            }
          }
          else           // lower_0 < back_0
          {

            // Medium pos_0 !

            // apply key transformation _1
            KeyType const pos_1 = key2index_1(arg_pos_1);
            KeyType lower_1_key;
            KeyType const frac_1 = std::modf(pos_1, &lower_1_key);
            DiffType const lower_1 = DiffType(lower_1_key);

            // Small pos_1 ?
            if (pos_1 < 0)
            {
              if (lower_1)
              {
                return;
              }
              else            // lower_1 == 0
              {
                KeyType const tmp_prod = KeyType(arg_value * (1 + frac_1));
                data_(lower_0 + 1, 0) += tmp_prod * frac_0;
                data_(lower_0, 0) += tmp_prod * (1 - frac_0);
                return;
              }
            }
            else             // pos_1 >= 0
            {
              DiffType const back_1 = data_.cols() - 1;
              // big pos_1
              if (lower_1 >= back_1)
              {
                if (lower_1 != back_1)
                {
                  return;
                }
                else                 // lower_1 == back_1
                {
                  KeyType const tmp_prod = KeyType(arg_value * (1 - frac_1));
                  data_(lower_0 + 1, lower_1) += tmp_prod * frac_0;
                  data_(lower_0, lower_1) += tmp_prod * (1 - frac_0);
                  return;
                }
              }
              else
              {
                // Medium pos_1 !

                // medium pos_0 and medium pos_1 --> "within" the matrix
                KeyType tmp_prod = KeyType(arg_value * frac_0);
                data_(lower_0 + 1, lower_1 + 1) += tmp_prod * frac_1;
                data_(lower_0 + 1, lower_1) += tmp_prod * (1 - frac_1);
                tmp_prod = KeyType(arg_value * (1 - frac_0));
                data_(lower_0, lower_1 + 1) += tmp_prod * frac_1;
                data_(lower_0, lower_1) += tmp_prod * (1 - frac_1);
                return;
              }
            }
          }
        }
      }

      //@}

      // ----------------------------------------------------------------------

      ///@name Discrete (non-interpolated) data
      //@{

      /// Returns the internal random access container storing the data.
      ContainerType & getData()
      {
        return data_;
      }

      /// Returns the internal random access container storing the data.
      ContainerType const & getData() const
      {
        return data_;
      }

      /**@brief Assigns data to the internal random access container storing
           the data.

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
        return data_.getEigenMatrix().size() == 0;
      }

      //@}

      // ----------------------------------------------------------------------

      ///\name Transformation
      //@{

      /// The transformation from "outside" to "inside" coordinates.
      KeyType key2index_0(KeyType pos) const
      {
        if (scale_0_)
        {
          pos -= offset_0_;
          pos /= scale_0_;
          return pos;
        }
        else
        {
          return 0;
        }
      }

      /// The transformation from "inside" to "outside" coordinates.
      KeyType index2key_0(KeyType pos) const
      {
        pos *= scale_0_;
        pos += offset_0_;
        return pos;
      }

      /// The transformation from "outside" to "inside" coordinates.
      KeyType key2index_1(KeyType pos) const
      {
        if (scale_1_)
        {
          pos -= offset_1_;
          pos /= scale_1_;
          return pos;
        }
        else
        {
          return 0;
        }
      }

      /// The transformation from "inside" to "outside" coordinates.
      KeyType index2key_1(KeyType pos) const
      {
        pos *= scale_1_;
        pos += offset_1_;
        return pos;
      }

      /// Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".
      KeyType const & getScale_0() const
      {
        return scale_0_;
      }

      /// Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".
      KeyType const & getScale_1() const
      {
        return scale_1_;
      }

      /**@brief Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".

      <b>Note:</b> Using this invalidates the inside and outside reference
      points.
      */
      void setScale_0(KeyType const & scale)
      {
        scale_0_ = scale;
      }

      /**@brief Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".

      <b>Note:</b> Using this invalidates the inside and outside reference
      points.
      */
      void setScale_1(KeyType const & scale)
      {
        scale_1_ = scale;
      }

      /// Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)".
      KeyType const & getOffset_0() const
      {
        return offset_0_;
      }

      /// Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)".
      KeyType const & getOffset_1() const
      {
        return offset_1_;
      }

      /**@brief Accessor.  "Offset" is the point (in "outside" units) which
           corresponds to "Data(0,0)".

           <b>Note:</b> Using this invalidates the inside and outside reference
           points.
      */
      void setOffset_0(KeyType const & offset)
      {
        offset_0_ = offset;
      }

      /**@brief Accessor.  "Offset" is the point (in "outside" units) which
           corresponds to "Data(0,0)".

           <b>Note:</b> Using this invalidates the inside and outside reference
           points.
      */
      void setOffset_1(KeyType const & offset)
      {
        offset_1_ = offset;
      }

      /**@brief Specifies the mapping from "outside" to "inside" coordinates by the following data:
           - <code>scale</code>: the difference in outside coordinates between consecutive values in the data vector.
           - <code>inside</code> and <code>outside</code>: these axis positions are mapped onto each other.

           For example, when you have a complicated probability distribution
           which is in fact centered around zero (but you cannot have negative
           indices in the data vector), then you can arrange things such that
           inside is the mean of the pre-computed, shifted density values of that
           distribution and outside is the centroid position of, say, a peak in
           the real world which you want to model by a scaled and shifted version
           of the probability distribution.

      */
      void setMapping_0(KeyType const & scale, KeyType const & inside_low, KeyType const & outside_low)
      {
        scale_0_ = scale;
        inside_0_ = inside_low;
        outside_0_ = outside_low;
        offset_0_ = outside_low - scale * inside_low;
        return;
      }

      /**@brief Specifies the mapping from "outside" to "inside" coordinates by the following data:
           - <code>inside_low</code> and <code>outside_low</code>: these axis positions are mapped onto each other.
           - <code>inside_high</code> and <code>outside_high</code>: these axis positions are mapped onto each other.

           This four argument version is just a convenience overload for the three argument version, which see.
      */
      void setMapping_0(KeyType const & inside_low, KeyType const & outside_low,
                        KeyType const & inside_high, KeyType const & outside_high)
      {
        if (inside_high != inside_low)
        {
          setMapping_0((outside_high - outside_low) / (inside_high - inside_low),
                       inside_low, outside_low);
        }
        else
        {
          setMapping_0(0, inside_low, outside_low);
        }
        return;
      }

      /**@brief Specifies the mapping from "outside" to "inside" coordinates by the following data:
           - <code>scale</code>: the difference in outside coordinates between consecutive values in the data vector.
           - <code>inside</code> and <code>outside</code>: these axis positions are mapped onto each other.

           For example, when you have a complicated probability distribution
           which is in fact centered around zero (but you cannot have negative
           indices in the data vector), then you can arrange things such that
           inside is the mean of the pre-computed, shifted density values of that
           distribution and outside is the centroid position of, say, a peak in
           the real world which you want to model by a scaled and shifted version
           of the probability distribution.

      */
      void setMapping_1(KeyType const & scale, KeyType const & inside_low, KeyType const & outside_low)
      {
        scale_1_ = scale;
        inside_1_ = inside_low;
        outside_1_ = outside_low;
        offset_1_ = outside_low - scale * inside_low;
        return;
      }

      /**@brief Specifies the mapping from "outside" to "inside" coordinates by the following data:
           - <code>inside_low</code> and <code>outside_low</code>: these axis positions are mapped onto each other.
           - <code>inside_high</code> and <code>outside_high</code>: these axis positions are mapped onto each other.

           This four argument version is just a convenience overload for the three argument version, which see.
      */
      void setMapping_1(KeyType const & inside_low, KeyType const & outside_low,
                        KeyType const & inside_high, KeyType const & outside_high)
      {
        if (inside_high != inside_low)
        {
          setMapping_1((outside_high - outside_low) / (inside_high - inside_low),
                       inside_low, outside_low);
        }
        else
        {
          setMapping_1(0, inside_low, outside_low);
        }
        return;
      }

      /// Accessor.  See setMapping().
      KeyType const & getInsideReferencePoint_0() const
      {
        return inside_0_;
      }

      /// Accessor.  See setMapping().
      KeyType const & getInsideReferencePoint_1() const
      {
        return inside_1_;
      }

      /// Accessor.  See setMapping().
      KeyType const & getOutsideReferencePoint_0() const
      {
        return outside_0_;
      }

      /// Accessor.  See setMapping().
      KeyType const & getOutsideReferencePoint_1() const
      {
        return outside_1_;
      }

      /// Lower boundary of the support, in "outside" coordinates.
      KeyType supportMin_0() const
      {
        return index2key_0(empty() ? KeyType(0.) : KeyType(-1.));
      }

      /// Lower boundary of the support, in "outside" coordinates.
      KeyType supportMin_1() const
      {
        return index2key_1(empty() ? KeyType(0.) : KeyType(-1.));
      }

      /// Upper boundary of the support, in "outside" coordinates.
      KeyType supportMax_0() const
      {
        return index2key_0(KeyType(data_.rows()));
      }

      /// Upper boundary of the support, in "outside" coordinates.
      KeyType supportMax_1() const
      {
        return index2key_1(KeyType(data_.cols()));
      }

      //@}

protected:

      /**@brief Data members*/
      //@{
      KeyType scale_0_;
      KeyType offset_0_;
      KeyType scale_1_;
      KeyType offset_1_;
      KeyType inside_0_;
      KeyType outside_0_;
      KeyType inside_1_;
      KeyType outside_1_;
      ContainerType data_;
      //@}
    };

  }   // namespace Math

} // namespace OpenMS

