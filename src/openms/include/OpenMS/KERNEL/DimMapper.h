// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once


#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/CommonEnums.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/MobilityPeak2D.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/RangeManager.h>

#include <array>
#include <memory>


namespace OpenMS
{
  /**
    @brief A base class for a dimension which represents a certain unit (e.g. RT or m/z).
           Derived classes implement virtual functions, which receive a well-defined data type,
           e.g. a Feature, and return the appropriate value for their dimension (the DimRT class would return the RT of the feature).
           This makes it possible to extract dimensions using a runtime configuration of DimBase instances.
           Very useful when mapping units (RT, m/z) to axis when plotting etc.

           The reverse (X-Y coordinates to data type, e.g. Peak1D) is also possible using 'from...()' methods
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

    /// D'tor (needs to be virtual; we are holding pointers to base in DimMapper)
    virtual ~DimBase() noexcept = default;

    /// Equality
    bool operator==(const DimBase& rhs) const
    {
      return unit_ == rhs.unit_;
    }

    /// Copy derived objects to avoid slicing when dealing with pointers to DimBase
    virtual std::unique_ptr<DimBase> clone() const = 0;

    virtual ValueType map(const Peak1D& p) const = 0;
    virtual ValueType map(const Peak2D& p) const = 0;
    virtual ValueType map(const ChromatogramPeak& p) const = 0;
    virtual ValueType map(const MSExperiment::ConstAreaIterator& it) const = 0;
    virtual ValueType map(const MobilityPeak1D& p) const = 0;
    virtual ValueType map(const MobilityPeak2D& p) const = 0;

    /// obtain value from a certain point in a spectrum
    virtual ValueType map(const MSSpectrum& spec, const Size index) const = 0;
    /// obtain value from a certain point in a chromatogram
    virtual ValueType map(const MSChromatogram& chrom, const Size index) const = 0;
    /// obtain value from a certain point in a mobilogram
    virtual ValueType map(const Mobilogram& mb, const Size index) const = 0;
    
    /// obtain vector of same length as @p spec; one element per peak
    /// @throw Exception::InvalidRange if elements do not support the dimension
    virtual ValueTypes map(const MSSpectrum& spec) const = 0;

    /// obtain vector of same length as @p spec; one element per peak
    /// @throw Exception::InvalidRange if elements do not support the dimension
    virtual ValueTypes map(const MSChromatogram& chrom) const = 0;

    virtual ValueType map(const BaseFeature& bf) const = 0;

    virtual ValueType map(const PeptideIdentification& pi) const = 0;

    /// Return the min/max (range) for a certain dimension
    virtual RangeBase map(const RangeAllType& rm) const = 0;

    /// Return the min/max (range) for a certain dimension (i.e. a reference to the base class of @p rm)
    virtual RangeBase& map(RangeAllType& rm) const = 0;

    /// Set the min/max (range) in @p out for a certain dimension
    virtual void setRange(const RangeBase& in, RangeAllType& out) const = 0;


    // from XY to a type

    /// set the dimension of a Peak1D
    virtual void fromXY(const ValueType in, Peak1D& p) const = 0;
    /// set the dimension of a ChromatogramPeak
    virtual void fromXY(const ValueType in, ChromatogramPeak& p) const = 0;
    /// set the dimension of a MobilityPeak1D
    virtual void fromXY(const ValueType in, MobilityPeak1D& p) const = 0;
    /// set the dimension of a MobilityPeak2D
    virtual void fromXY(const ValueType in, MobilityPeak2D& p) const = 0;

    /// Name of the dimension, e.g. 'RT [s]' 
    std::string_view getDimName() const
    {
      return DIM_NAMES[(int)unit_];
    }

    /// Name of the dimension, e.g. 'RT'
    std::string_view getDimNameShort() const
    {
      return DIM_NAMES_SHORT[(int)unit_];
    }

    /// The unit of the dimension
    DIM_UNIT getUnit() const
    {
      return unit_;
    }

    /**
     * \brief Creates a short string representation with "UNIT: value", where value has a predefined precision (see valuePrecision())
     * \param value The value of this Dim to format
     * \return A formatted string, e.g. "RT: 45.32"
     */
    String formattedValue(const ValueType value) const;

    /// like formattedValue() but with a custom unit prefix instead of the default one for the dim, e.g. "myText: 45.32" 
    String formattedValue(ValueType value, const String& prefix) const;

    /// return the recommended precision for the current unit (2 digits for RT, 8 for m/z, etc)
    int valuePrecision() const;

  protected:
    DIM_UNIT unit_; ///< the unit of this dimension    
  };



  class OPENMS_DLLAPI DimRT final : public DimBase
  {
  public:
    DimRT() : DimBase(DIM_UNIT::RT) {}

    std::unique_ptr<DimBase> clone() const override
    {
      return std::make_unique<DimRT>();
    }

    ValueType map(const Peak1D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(const Peak2D& p) const override
    {
      return p.getRT();
    }
    ValueType map(const ChromatogramPeak& p) const override
    {
      return p.getRT();
    }
    ValueType map(const MSSpectrum& spec, const Size /*index*/) const override
    {
      return spec.getRT();
    }
    ValueType map(const MSChromatogram& chrom, const Size index) const override
    {
      return chrom[index].getRT();
    }
    ValueType map(const Mobilogram& mb, const Size /*index*/) const override
    {
      return mb.getRT();
    }

    ValueTypes map(const MSSpectrum&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueTypes map(const MSChromatogram& chrom) const override
    {
      ValueTypes res;
      res.reserve(chrom.size());
      for (const auto& p : chrom)
      {
        res.push_back(p.getRT());
      }
      return res;
    }

    ValueType map(const MSExperiment::ConstAreaIterator& it) const override
    {
      return it.getRT();
    }
    ValueType map(const MobilityPeak1D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(const MobilityPeak2D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    ValueType map(const BaseFeature& bf) const override
    {
      return bf.getRT();
    }

    ValueType map(const PeptideIdentification& pi) const override
    {
      return pi.getRT();
    }

    RangeBase map(const RangeAllType& rm) const override
    {
      return rm.getRangeForDim(MSDim::RT);
    }
    RangeBase& map(RangeAllType& rm) const override
    {
      return rm.getRangeForDim(MSDim::RT);
    }

    void setRange(const RangeBase& in, RangeAllType& out) const override
    {
      out.RangeRT::operator=(in);
    }

    /// set the RT of a Peak1D (throws)
    void fromXY(const ValueType, Peak1D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /// set the RT of a ChromatogramPeak
    void fromXY(const ValueType in, ChromatogramPeak& p) const override
    {
      p.setRT(in);
    }
    /// set the RT of a MobilityPeak1D (throws)
    void fromXY(const ValueType, MobilityPeak1D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    /// set the RT of a MobilityPeak2D (throws)
    void fromXY(const ValueType, MobilityPeak2D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
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
    ValueType map(const ChromatogramPeak&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(const MSExperiment::ConstAreaIterator& it) const override
    {
      return it->getMZ();
    }
    ValueType map(const MobilityPeak1D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(const MobilityPeak2D& p) const override
    {
      return p.getMZ();
    }

    ValueType map(const MSSpectrum& spec, const Size index) const override
    {
      return spec[index].getMZ();
    }
    ValueType map(const MSChromatogram& chrom, const Size /*index*/) const override
    {
      return chrom.getPrecursor().getMZ();
    }
    ValueType map(const Mobilogram& /*mb*/, const Size /*index*/) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    ValueTypes map(const MSSpectrum& spec) const override
    {
      ValueTypes res;
      res.reserve(spec.size());
      for (const auto& p : spec)
      {
        res.push_back(p.getMZ());
      }
      return res;
    }
    ValueTypes map(const MSChromatogram&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }  
    
    ValueType map(const BaseFeature& bf) const override
    {
      return bf.getMZ();
    }

    ValueType map(const PeptideIdentification& pi) const override
    {
      return pi.getMZ();
    }

    RangeBase map(const RangeAllType& rm) const override
    {
      return rm.getRangeForDim(MSDim::MZ);
    }
    RangeBase& map(RangeAllType& rm) const override
    {
      return rm.getRangeForDim(MSDim::MZ);
    }

    void setRange(const RangeBase& in, RangeAllType& out) const override
    {
      out.RangeMZ::operator=(in);
    }

    /// set the MZ of a Peak1D
    void fromXY(const ValueType in, Peak1D& p) const override
    {
      p.setMZ(in);
    }

    /// set the MZ of a ChromatogramPeak (throws)
    void fromXY(const ValueType, ChromatogramPeak&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /// set the MZ of a MobilityPeak1D (throws)
    void fromXY(const ValueType, MobilityPeak1D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    /// set the MZ of a MobilityPeak2D (throws)
    void fromXY(const ValueType in, MobilityPeak2D& p) const override
    {
      p.setMZ(in);
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
    ValueType map(const ChromatogramPeak& p) const override
    {
      return p.getIntensity();
    }
    ValueType map(const MSExperiment::ConstAreaIterator& it) const override
    {
      return it->getIntensity();
    }
    ValueType map(const MobilityPeak1D& p) const override
    {
      return p.getIntensity();
    }
    ValueType map(const MobilityPeak2D& p) const override
    {
      return p.getIntensity();
    }

    ValueType map(const MSSpectrum& spec, const Size index) const override
    {
      return spec[index].getIntensity();
    }
    ValueType map(const MSChromatogram& chrom, const Size index) const override
    {
      return chrom[index].getIntensity();
    }
    ValueType map(const Mobilogram& mb, const Size index) const override
    {
      return mb[index].getIntensity();
    }

    ValueTypes map(const MSSpectrum& spec) const override
    {
      ValueTypes res;
      res.reserve(spec.size());
      for (const auto& p : spec)
      {
        res.push_back(p.getIntensity());
      }
      return res;
    }

    ValueTypes map(const MSChromatogram& chrom) const override
    {
      ValueTypes res;
      res.reserve(chrom.size());
      for (const auto& p : chrom)
      {
        res.push_back(p.getIntensity());
      }
      return res;
    }

    ValueType map(const BaseFeature& bf) const override
    {
      return bf.getIntensity();
    }

    ValueType map(const PeptideIdentification&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    
    RangeBase map(const RangeAllType& rm) const override
    {
      return rm.getRangeForDim(MSDim::INT);
    }
    RangeBase& map(RangeAllType& rm) const override
    {
      return rm.getRangeForDim(MSDim::INT);
    }

    void setRange(const RangeBase& in, RangeAllType& out) const override
    {
      out.RangeIntensity::operator=(in);
    }

    /// set the intensity of a Peak1D
    void fromXY(const ValueType in, Peak1D& p) const override
    {
      p.setIntensity(Peak1D::IntensityType(in));
    }

    /// set the intensity of a ChromatogramPeak
    void fromXY(const ValueType in, ChromatogramPeak& p) const override
    {
      p.setIntensity(ChromatogramPeak::IntensityType(in));
    }
    /// set the intensity of a MobilityPeak1D
    void fromXY(const ValueType in, MobilityPeak1D& p) const override
    {
      p.setIntensity(MobilityPeak1D::IntensityType(in));
    }
    /// set the intensity of a MobilityPeak2D
    void fromXY(const ValueType in, MobilityPeak2D& p) const override
    {
      p.setIntensity(MobilityPeak2D::IntensityType(in));
    }
  };

  class OPENMS_DLLAPI DimIM final : public DimBase
  {
  public:
    DimIM(const DIM_UNIT im_unit) : DimBase(im_unit) {}

    std::unique_ptr<DimBase> clone() const override
    {
      return std::make_unique<DimIM>(*this);
    }

    ValueType map(const Peak1D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(const Peak2D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(const ChromatogramPeak&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueTypes map(const MSSpectrum&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueTypes map(const MSChromatogram&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    ValueType map(const MSExperiment::ConstAreaIterator& it) const override
    {
      return it.getDriftTime();
    }

    ValueType map(const MobilityPeak1D& p) const override
    {
      return p.getMobility();
    }
    ValueType map(const MobilityPeak2D& p) const override
    {
      return p.getMobility();
    }

    ValueType map(const MSSpectrum& spec, const Size /*index*/) const override
    {
      return spec.getDriftTime();
    }
    ValueType map(const MSChromatogram&, const Size) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ValueType map(const Mobilogram& mb, const Size index) const override
    {
      return mb[index].getMobility();
    }

    ValueType map(const BaseFeature&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    ValueType map(const PeptideIdentification&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    RangeBase map(const RangeAllType& rm) const override
    {
      return rm.getRangeForDim(MSDim::IM);
    }
    RangeBase& map(RangeAllType& rm) const override
    {
      return rm.getRangeForDim(MSDim::IM);
    }

    void setRange(const RangeBase& in, RangeAllType& out) const override
    {
      out.RangeMobility::operator=(in);
    }

    /// set the IM of a Peak1D (throws)
    void fromXY(const ValueType, Peak1D&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /// set the IM of a ChromatogramPeak (throws)
    void fromXY(const ValueType, ChromatogramPeak&) const override
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /// set the IM of a MobilityPeak1D
    void fromXY(const ValueType in, MobilityPeak1D& p) const override
    {
      p.setMobility(in);
    }
    /// set the IM of a MobilityPeak2D
    void fromXY(const ValueType in, MobilityPeak2D& p) const override
    {
      p.setMobility(in);
    }
  };
  
  /// Declare Dimensions for referencing dimensions in plotting etc
  /// e.g. Point[X], Point[Y]
  /// The order X,Y,Z,... is important here. Some classes rely upon it.
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
    Point map(const T& data) const
    {
      Point pr;
      for (int i = 0; i < N_DIM; ++i) pr[i] = dims_[i]->map(data);
      return pr;
    }
    /// convert an OpenMS datapoint in a container (such as MSSpectrum) to an N_DIM-dimensional point
    template<typename Container>
    Point map(const Container& data, const Size index) const
    {
      Point pr;
      for (int i = 0; i < N_DIM; ++i)
        pr[i] = dims_[i]->map(data, index);
      return pr;
    }

    /// Convert Range to an N_DIM-dimensional area (min and max for each dimension)
    template<typename ...Ranges>
    DRange<N_DIM> mapRange(const RangeManager<Ranges...>& ranges) const
    {
      DRange<N_DIM> res;
      RangeAllType all;
      all.assign(ranges);
      for (int i = 0; i < N_DIM; ++i)
      {
        RangeBase mm = dims_[i]->map(all);
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
    /// The range will only span a single value in each dimension covered by this mapper.
    /// Dimensions not contained in this DimMapper will remain untouched in @p output
    template<typename... Ranges>
    void fromXY(const Point& in, RangeManager<Ranges...>& output) const
    {
      for (int i = 0; i < N_DIM; ++i)
      {
        dims_[i]->setRange({in[i], in[i]}, output);
      }
    }

    /// Convert an N_DIM-Point to a Peak1D or ChromatogramPeak.
    /// Dimensions not contained in this DimMapper will remain untouched in @p out
    /// @throw Exception::InvalidRange if DimMapper has a dimension not supported by T
    template<typename T>
    void fromXY(const Point& in, T& out) const
    {
      for (int i = 0; i < N_DIM; ++i)
      {
        dims_[i]->fromXY(in[i], out);
      }
    }

    /// Convert an N_DIM-Point to a Range with all known dimensions.
    /// The range will only span a single value in each dimension covered by this mapper.
    /// Dimensions not contained in this DimMapper will remain empty.
    RangeAllType fromXY(const Point& in) const
    {
      RangeAllType output;
      for (int i = 0; i < N_DIM; ++i)
      {
        dims_[i]->setRange({in[i], in[i]}, output);
      }
      return output;
    }

    /// obtain unit/name for X/Y/Z dimension.
    const DimBase& getDim(DIM d) const
    {
      assert((int)d <= N_DIM);
      return *dims_[(int)d];
    }

    bool hasUnit(DIM_UNIT unit) const
    {
      for (int i = 0; i < N_DIM; ++i)
      {
        if (dims_[i]->getUnit() == unit) return true;
      }
      return false;
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
        case DIM_UNIT::FAIMS_CV:
        case DIM_UNIT::IM_MS:
        case DIM_UNIT::IM_VSSC:
          return std::make_unique<DimIM>(u);
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
    /// The Area in X,Y,(Z)... dimension (number of dimensions depends on N_DIM)
    using AreaXYType = DRange<N_DIM>;

    /// No default C'tor
    Area() = delete;
    
    /// Custom C'tor with a mapper (non owning pointer)
    Area(const DimMapper<N_DIM>* const dims) 
      : mapper_(dims)
    {
    }

    /// Copy C'tor
    Area(const Area& range) = default;

    /// Assignment operator - which checks for identical DimMappers and throws otherwise
    Area& operator=(const Area& rhs)
    {
      // check that Dims are identical, otherwise this is very dangerous (the user probably only wanted to update the area, not change its mapping).
      if (mapper_ != rhs.mapper_ && *mapper_ != *rhs.mapper_)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Assignment of Areas using different mappers!");
      }
      data_range_ = rhs.data_range_;
      visible_area_ = rhs.visible_area_;
      return *this;
    } 

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
    const Area& setArea(const RangeAllType& data)
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

    const RangeAllType& getAreaUnit() const
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
    Area cloneWith(const RangeAllType& data) const
    {
      Area clone(*this);
      clone.setArea(data);
      return clone;
    }

    /**
     * \brief Push the area into a sandbox (if its outside the sandbox). See UnitRange::pushInto()
     * \param sandbox The sandbox which delimits the range of this area
     */
    void pushInto(const RangeAllType& sandbox)
    {
      auto a = data_range_;
      a.pushInto(sandbox);
      setArea(a);
    }

    /// empty all dimensions
    void clear()
    {
      setArea(RangeAllType());
    }

  private:
    /* two sides of the same coin... */
    RangeAllType data_range_;                     ///< range in units
    AreaXYType visible_area_ = AreaXYType::empty; ///< range in terms of axis (X and Y axis)
    /// and a mapper (non-owning pointer) to translate between the two
    const DimMapper<N_DIM>* mapper_;
  };

} // namespace OpenMS
