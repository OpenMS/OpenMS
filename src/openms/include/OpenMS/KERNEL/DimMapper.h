// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/RangeManager.h>

#include <array>
#include <memory>
#include <string_view>


namespace OpenMS
{
  /// Enum for different units which can be displayed on a plotting axis
  /// The order is arbitrary.
  enum class DIM_UNIT
  {
    RT = 0,     ///< RT in seconds
    MZ,         ///< m/z
    INT,        ///< intensity
    IM_MS,      ///< ion mobility milliseconds
    IM_VSSC,    ///< volt-second per square centimeter (i.e. 1/K_0)
    FAIMS_CM,   ///< FAIMS compensation voltage
    SIZE_OF_DIM_UNITS
  };
  std::string_view DIM_NAMES[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {"RT [s]", "m/z [Th]", "intensity", "IM [milliseconds]", "IM [vs / cm2]", "FAIMS CV"};


  /**
    @brief A base class for a dimension which represents a certain unit (e.g. RT or m/z).
           Derived classes implement virtual functions, which receive a well-defined data type,
           e.g. a Feature, and return the appropriate value for their dimension (the DimRT class would return the RT of the feature).
           This makes it possible to extract dimensions using a runtime configuration of DimBase instances.
           Very useful when mapping units (RT, m/z) to axis when plotting etc.  
  */
  class OPENMS_DLLAPI DimBase
  {
  public:
    using ValueType = double;
    using ValueTypes = std::vector<ValueType>;

    /// No default c'tor
    DimBase() = delete;

    /// Custom c'tor with unit
    DimBase(DIM_UNIT unit) :
        unit_(unit) 
    {}

    /// Assignment operator
    DimBase& operator=(const DimBase& rhs) = default;

    /// Equality
    bool operator==(const DimBase& rhs) const
    {
      return unit_ == rhs.unit_;
    }

    /// Copy derived objects to avoid slicing when dealing with pointers to DimBase
    virtual std::unique_ptr<DimBase> clone() const = 0;

    virtual ValueType map(const Peak1D& p) const = 0;
    virtual ValueType map(const Peak2D& p) const = 0;
    virtual ValueType map(MSExperiment::ConstAreaIterator it) const = 0;
    
    /// obtain vector of same length as @p spec; one element per peak
    /// @throw Exception::InvalidRange if elements do not support the dimension
    virtual ValueTypes map(const MSSpectrum& spec) const = 0;

    virtual ValueType map(const BaseFeature& bf) const = 0;

    virtual ValueType map(const PeptideIdentification& pi) const = 0;


    /// Return the min/max (range) for a certain dimension
    virtual RangeBase map(const RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& rm) const = 0;

    /// Set the min/max (range) in @p rm for a certain dimension
    virtual void setRange(const RangeBase& in, RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& out) const = 0;

    /// Name of the dimension, e.g. 'RT [s]' 
    std::string_view getDimName() const
    {
      return DIM_NAMES[(int)unit_];
    }

    /// The unit of the dimension
    DIM_UNIT getUnit() const
    {
      return unit_;
    }

  protected:
    DIM_UNIT unit_; ///< the unit of this dimension    
  };


  class OPENMS_DLLAPI DimRT final : public DimBase
  {
  public:
    DimRT() : DimBase(DIM_UNIT::RT) {};

    std::unique_ptr<DimBase> clone() const override
    {
      return std::make_unique<DimRT>();
    }

    ValueType map(const Peak1D& p) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(const Peak2D& p) const override
    {
      return p.getRT();
    }
    ValueTypes map(const MSSpectrum& spec) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(MSExperiment::ConstAreaIterator it) const override
    {
      return it.getRT();
    }

    ValueType map(const BaseFeature& bf) const override
    {
      return bf.getRT();
    }

    ValueType map(const PeptideIdentification& pi) const override
    {
      return pi.getRT();
    }

    RangeBase map(const RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& rm) const override
    {
      return (RangeRT)rm;
    }

    void setRange(const RangeBase& in, RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& rm) const
    {
      rm.RangeRT::operator=(in);
    }

  };
  class OPENMS_DLLAPI DimMZ final : public DimBase
  {
  public:
    DimMZ() : DimBase(DIM_UNIT::MZ) {};
    
    std::unique_ptr<DimBase> clone() const override
    {
      return std::make_unique<DimMZ>();
    }
          
    ValueType map(const Peak1D& p) const override
    {
      return p.getMZ();
    }
    ValueType map(const Peak2D& p) const override
    {
      return p.getMZ();
    }
    ValueType map(MSExperiment::ConstAreaIterator it) const override
    {
      return it->getMZ();
    }

    ValueTypes map(const MSSpectrum& spec) const override
    {
      ValueTypes res;
      res.reserve(spec.size());
      for (const auto& p : spec)
        res.push_back(p.getMZ());
      return res;
    }  
    
    ValueType map(const BaseFeature& bf) const override
    {
      return bf.getMZ();
    }

    ValueType map(const PeptideIdentification& pi) const override
    {
      return pi.getMZ();
    }

    RangeBase map(const RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& rm) const override
    {
      return RangeMZ(rm);
    }

    void setRange(const RangeBase& in, RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& rm) const
    {
      rm.RangeMZ::operator=(in);
    }


  };
  class OPENMS_DLLAPI DimINT final : public DimBase
  {
  public:
    DimINT() : DimBase(DIM_UNIT::INT) {};
    
    std::unique_ptr<DimBase> clone() const override
    {
      return std::make_unique<DimINT>();
    }
    
    ValueType map(const Peak1D& p) const override
    {
      return p.getIntensity();
    }
    ValueType map(const Peak2D& p) const override
    {
      return p.getIntensity();
    }
    ValueType map(MSExperiment::ConstAreaIterator it) const override
    {
      return it->getIntensity();
    }

    ValueTypes map(const MSSpectrum& spec) const override
    {
      ValueTypes res;
      res.reserve(spec.size());
      for (const auto& p : spec)
        res.push_back(p.getIntensity());
      return res;
    }

    ValueType map(const BaseFeature& bf) const override
    {
      return bf.getIntensity();
    }

    ValueType map(const PeptideIdentification& pi) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    
    RangeBase map(const RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& rm) const override
    {
      return RangeIntensity(rm);
    }

    void setRange(const RangeBase& in, RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& rm) const
    {
      rm.RangeIntensity::operator=(in);
    }
  };
  
  /// make axis label of returned Point explicit
  /// e.g. Point[X], Point[Y]
  enum class DIM
  {
    X = 0,
    Y = 1,
    Z = 2
  };

  /**
      @brief Allows dynamical switching (at runtime) between a dimension (RT, m/z, int, IM, etc) 
             and X,Y,Z coordinates. You can set either of them, and query the other.
             The Mapping is stored internally.
             The unit to which the X,Y,Z coordinates currently mapped onto can also be queried (useful for axis labels etc).

      Use N_DIM template parameter to determine the number of axis dimensions (1-3 is currently supported). Usually 2 or 3 make sense.

      @ingroup Visual
  */
  template<int N_DIM>
  class DimMapper
  {
  public:
    using Point = DPosition<N_DIM, DimBase::ValueType>;

    /// No default c'tor (we need dimensions)
    DimMapper() = delete;

    /// Custom C'tor with given dimensions to map to (the order is assumed to be X, Y, Z, ...)
    DimMapper(const DIM_UNIT (&units)[N_DIM])
      :dims_([&]() {
          std::array<std::unique_ptr<const DimBase>, N_DIM> dims_tmp;
          for (int i = 0; i < N_DIM; ++i)
          {
            dims_tmp[i] = create_(units[i]);
          }
          return dims_tmp;
      }()) // immediately evaluated lambda to enable 'dims_' to be const
    {
      static_assert(N_DIM >= 1); // at least one dimension (X)
      static_assert(N_DIM <= 3); // at most three (X, Y, Z)
    }

    /// Copy C'tor
    DimMapper(const DimMapper& rhs) // cannot be defaulted due to unique_ptr
    {
      *this = rhs;
    };

    /// Assignment operator
    DimMapper& operator=(const DimMapper& rhs) // cannot be defaulted due to unique_ptr
    {
      for (int i = 0; i < N_DIM; ++i) dims_[i] = rhs.dims_[i]->clone(); 
      return *this;
    };

    /// Equality 
    bool operator==(const DimMapper& rhs) const
    {
      bool res {true};
      for (int i = 0; i < N_DIM; ++i)
      {
        res &= (*dims_[i] == *rhs.dims_[i]); 
      }
      return res;
    }
    
    /// Inequality
    bool operator!=(const DimMapper& rhs) const
    {
      return !operator==(rhs);
    }

    /// convert an OpenMS datatype (such as Feature) to an N_DIM-dimensional point
    template <typename T>
    Point map(const T& data)
    {
      Point pr;
      for (int i = 0; i < N_DIM; ++i) pr[i] = dims_[i]->map(data);
      return pr;
    }

    /// Convert Range to an N_DIM-dimensional area (min and max for each dimension)
    template<typename ...Ranges>
    DRange<N_DIM> mapRange(const RangeManager<Ranges...>& ranges) const
    {
      DRange<N_DIM> res;
      using Coord = typename DRange<N_DIM>::PositionType;
      for (int i = 0; i < N_DIM; ++i)
      {
        RangeBase mm = dims_[i]->map(ranges);
        if (mm.isEmpty()) continue;
        res.setDimMinMax(i, {mm.getMin(), mm.getMax()});
      }
      return res;
    }

    /// Convert an N_DIM-dimensional area (min and max for each dimension) to a Range.
    /// Empty dimensions in the input @p in, will also be made empty in @p output.
    /// Dimensions not contained in this DimMapper will remain untouched in @p output
    template<typename... Ranges>
    void fromXY(const DRange<N_DIM>& in, RangeManager<Ranges...>& output) const
    {
      for (int i = 0; i < N_DIM; ++i)
      {
        if (in.isEmpty(i))
          dims_[i]->setRange(RangeBase(), output);
        else
          dims_[i]->setRange({in.minPosition()[i], in.maxPosition()[i]}, output);
      }
    }

    /// Convert an N_DIM-Point to a Range.
    /// Dimensions not contained in this DimMapper will remain untouched in @p output
    template<typename... Ranges>
    void fromXY(const Point& in, RangeManager<Ranges...>& output) const
    {
      for (int i = 0; i < N_DIM; ++i)
      {
        dims_[i]->setRange({in[i], in[i]}, output);
      }
    }

    /// obtain unit/name for X/Y/Z dimension.
    const DimBase& getDim(DIM d) const
    {
      assert((int)d <= N_DIM);
      return *dims_[(int)d];
    }

  protected:
    /// a minimal factory
    static std::unique_ptr<const DimBase> create_(DIM_UNIT u)
    {
      switch (u)
      {
        case DIM_UNIT::RT:
          return std::make_unique<DimRT>();
        case DIM_UNIT::MZ:
          return std::make_unique<DimMZ>();
        case DIM_UNIT::INT:
          return std::make_unique<DimINT>();
        default:
          throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    }

    std::array<std::unique_ptr<const DimBase>, N_DIM> dims_; ///< mappers for the X,Y,Z... dimension
  };


  /// The data is stored in two members, one axis-related (X and Y; unit does not matter), and one unit-related (units; no mapping to axis)
  /// You can set either, and the other will be updated accordingly as long as you provide a DimMapper which translates between the two representations.
  template <int N_DIM>
  class Area
  {
  public:
    using UnitRange = RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>;
    /// The Area in X,Y,(Z)... dimension (number of dimensions depends on N_DIM)
    using AreaXYType = DRange<N_DIM>;

    /// No default C'tor
    Area() = delete;
    
    /// Custom C'tor with a mapper
    Area(const DimMapper<N_DIM>* const dims) 
      : mapper_(dims)
    {
    }

    /// Copy C'tor
    Area(const Area& range) = default;

    /// Assignment operator
    Area& operator=(const Area& rhs) = default; 

    bool operator==(const Area& rhs) const
    {
      return data_range_ == rhs.data_range_
            && visible_area_ == rhs.visible_area_
            && (*mapper_ == *rhs.mapper_);
    }
    bool operator!=(const Area& rhs) const
    {
      return !operator==(rhs);
    }

    /**
       @brief Set the area using unit data (RT, m/z, ...)

       @param data Area in units
    */
    const Area& setArea(const UnitRange& data)
    {
      data_range_ = data;
      // update axis view using dims
      visible_area_ = mapper_->mapRange(data);
      return *this;
    }

    /**
       @brief Set the area using axis data (X and Y)

       @param data Area as displayed on the axis
    */
    const Area& setArea(const AreaXYType& data)
    {
      visible_area_ = data;
      // update range view from XY area using dims
      mapper_->fromXY(visible_area_, data_range_);
      return *this;
    }

    const AreaXYType& getAreaXY() const
    {
      return visible_area_;
    }

    const UnitRange& getAreaUnit() const
    {
      return data_range_;
    }

    /**
      @brief Clone the current object, set the area of the clone using axis data (X and Y) and return the clone.

      @param data New area as displayed on the axis
    */
    Area cloneWith(const AreaXYType& data) const
    {
      Area clone(*this);
      clone.setArea(data);
      return clone;
    }

    /**
      @brief Clone the current object, set the area of the clone using unit data (RT, m/z, ...) and return the clone.

      @param data New area in units
    */
    Area cloneWith(const UnitRange& data) const
    {
      Area clone(*this);
      clone.setArea(data);
      return clone;
    }

  private:
    /* two sides of the same coin... */
    UnitRange data_range_;                        ///< range in units
    AreaXYType visible_area_ = AreaXYType::empty; ///< range in terms of axis (X and Y axis)
    /// and a mapper (non-owning pointer) to translate between the two
    const DimMapper<N_DIM>* mapper_;
  };

} // namespace OpenMS
