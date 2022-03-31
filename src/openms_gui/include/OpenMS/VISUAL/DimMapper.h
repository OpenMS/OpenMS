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

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>

#include <array>
#include <memory>
#include <string_view>


namespace OpenMS
{
  /// Enum for different units which can be displayed on a plotting axis
  /// The order is arbitrary.
  enum class DIM_UNITS
  {
    RT = 0,     ///< RT in seconds
    MZ,         ///< m/z
    INT,        ///< intensity
    IM_MS,      ///< ion mobility milliseconds
    IM_VSSC,    ///< volt-second per square centimeter (i.e. 1/K_0)
    FAIMS_CM,   ///< FAIMS compensation voltage
    SIZE_OF_DIM_UNITS
  };
  std::string_view DIM_NAMES[(int)DIM_UNITS::SIZE_OF_DIM_UNITS] = {"RT [s]", "m/z [Th]", "intensity", "IM [milliseconds]", "IM [vs / cm2]", "FAIMS CV"};

  class OPENMS_GUI_DLLAPI DimBase
  {
  public:
    using ValueType = float;
    using ValueTypes = std::vector<ValueType>;
    
    DimBase() = delete;
    DimBase(DIM_UNITS unit) :
        unit_(unit) 
    {}

    virtual ValueType map(const Peak1D&p) const = 0;
    virtual ValueType map(const Peak2D& p) const = 0;
    /// obtain vector of same length as @p spec; one element per peak
    /// @throw Exception::InvalidRange if elements do not support the dimension
    virtual ValueTypes map(const MSSpectrum& spec) const = 0;

    std::string_view getDimName() const
    {
      return DIM_NAMES[(int)unit_];
    }

  protected:
    DIM_UNITS unit_;    
  };


  class OPENMS_GUI_DLLAPI DimRT final : public DimBase
  {
  public:
    DimRT() : DimBase(DIM_UNITS::RT) {};
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
  };
  class OPENMS_GUI_DLLAPI DimMZ final : public DimBase
  {
  public:
    DimMZ() : DimBase(DIM_UNITS::MZ) {};
    ValueType map(const Peak1D& p) const override
    {
      return p.getMZ();
    }
    ValueType map(const Peak2D& p) const override
    {
      return p.getMZ();
    }
    ValueTypes map(const MSSpectrum& spec) const override
    {
      ValueTypes res;
      res.reserve(spec.size());
      for (const auto& p : spec)
        res.push_back(p.getMZ());
      return res;
    }  
  };
  class OPENMS_GUI_DLLAPI DimINT final : public DimBase
  {
  public:
    DimINT() : DimBase(DIM_UNITS::INT) {};
    ValueType map(const Peak1D& p) const override
    {
      return p.getIntensity();
    }
    ValueType map(const Peak2D& p) const override
    {
      return p.getIntensity();
    }
    ValueTypes map(const MSSpectrum& spec) const override
    {
      ValueTypes res;
      res.reserve(spec.size());
      for (const auto& p : spec)
        res.push_back(p.getIntensity());
      return res;
    }  
  };
  
  /**
      @brief Allows dynamically (at runtime) switching which dimension (RT, m/z, int, IM, etc) 
             gets mapped onto X,Y,Z coordinates when plotting.

      @ingroup Visual
  */
  template<int N_DIM>
  class DimMapper
  {
  public:
    /// make axis label of returned Point explicit
    /// e.g. Point[X], Point[Y]
    enum DIM
    {
      X = 0,
      Y = 1,
      Z = 1
    };

    using Point = DPosition<N_DIM, DimBase::ValueType>;

    DimMapper(const DIM_UNITS (&units)[N_DIM])
      :dims_([&]() {
          std::array<std::unique_ptr<const DimBase>, N_DIM> dims_tmp;
          for (int i = 0; i < N_DIM; ++i)
          {
            dims_tmp[i] = create_(units[i]);
          }
          return dims_tmp;
      }()) // immediately evaluated lambda to enable 'dims_' to be const
    {
    }

    template <typename T>
    Point map(const T& data)
    {
      Point pr;
      for (int i = 0; i < N_DIM; ++i) pr[i] = dims_[i]->map(data);
      return pr;
    }

  protected:
    /// a minimal factory
    static std::unique_ptr<const DimBase> create_(DIM_UNITS u)
    {
      switch (u)
      {
        case DIM_UNITS::RT:
          return std::make_unique<DimRT>();
        case DIM_UNITS::MZ:
          return std::make_unique<DimMZ>();
        case DIM_UNITS::INT:
          return std::make_unique<DimINT>();
        default:
          throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    }

    const std::array<std::unique_ptr<const DimBase>, N_DIM> dims_; ///< mappers for the X,Y,Z... dimension
  };
} // namespace OpenMS
