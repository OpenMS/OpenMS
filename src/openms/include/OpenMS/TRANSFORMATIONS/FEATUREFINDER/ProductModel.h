// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PRODUCTMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PRODUCTMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

namespace OpenMS
{

  /**
      @brief Class for product models i.e. models with D independent dimensions

      The predicted intensity is simply the product of the intensities in each dimension

      @htmlinclude OpenMS_ProductModel.parameters

  @todo This class provides new member functions, which makes
  Factory<BaseModel<2> >::create("ProductModel2D") pretty much useless!  (Clemens)

      @ingroup FeatureFinder
  */
  template <UInt D>
  class ProductModel;

  /// The class template is only implemented for D=2 because we use Peak2D here.
  template <>
  class OPENMS_DLLAPI ProductModel<2>:
    public BaseModel<2>
  {
public:

    /// Number of dimensions
    enum {D = 2};

    typedef double IntensityType;
    typedef DPosition<D> PositionType;
    typedef BaseModel<D>::SamplesType SamplesType;

    /// Default constructor
    ProductModel() :
      BaseModel<D>(),
      distributions_(D, nullptr)
    {
      this->setName(this->getProductName());

      //Register model info
      for (UInt dim = 0; dim < D; ++dim)
      {
        String name = Peak2D::shortDimensionName(dim);
        this->subsections_.push_back(name);
        this->defaults_.setValue(name, "GaussModel", "Name of the model used for this dimension");
      }

      //defaults
      this->defaults_.setValue("intensity_scaling", 1.0, "Scaling factor used to adjust the model distribution to the intensities of the data");
      this->defaultsToParam_();
    }

    /// copy constructor
    ProductModel(const ProductModel & source) :
      BaseModel<D>(source),
      distributions_(D, nullptr),
      scale_(source.scale_)
    {
      for (UInt dim = 0; dim < D; ++dim)
      {
        // clone source model
        if (source.distributions_[dim])
        {
          ModelDescription<1> desc(source.distributions_[dim]);
          setModel(dim, desc.createModel());
        }
      }
      updateMembers_();
    }

    /// destructor
    ~ProductModel() override
    {
      for (Size dim = 0; dim < D; ++dim)
      {
        delete distributions_[dim];
      }
    }

    /// assignment operator
    virtual ProductModel & operator=(const ProductModel & source)
    {
      if (&source == this) return *this;

      BaseModel<D>::operator=(source);
      scale_ = source.scale_;

      for (UInt dim = 0; dim < D; ++dim)
      {
        if (source.distributions_[dim])
        {
          // clone source model
          ModelDescription<1> desc(source.distributions_[dim]);
          setModel(dim, desc.createModel());
        }
        else
        {
          distributions_[dim] = nullptr;
        }
      }
      updateMembers_();

      return *this;
    }

    /// intensity equals product of intensities in each dimension
    IntensityType getIntensity(const PositionType & pos) const override
    {
      IntensityType intens(scale_);
      for (UInt dim = 0; dim < D; ++dim)
      {
        if (distributions_[dim] == nullptr)
        {
          throw Exception::BaseException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("ProductModel: model for dimension ") + dim + " not set.", "");
        }
        intens *= distributions_[dim]->getIntensity(pos[dim]);
      }
      return intens;
    }

    /// create new ProductModel object (needed by Factory)
    static BaseModel<D> * create()
    {
      return new ProductModel<D>();
    }

    /// Returns the name of the model
    static const String getProductName()
    {
      return String("ProductModel") + D + "D";
    }

    /** @brief set model @p dist for dimension @p dim

            @p dist is supposed to be allocated by new because it will be freed if
            ProductModel is destroyed or the model for that dimensions changes.
            For that reason no model @p dist should be assigned to multiple ProductModels.<br>
            ProductModel parameters are set when calling ProductModel::getParameters().
    */
    ProductModel & setModel(UInt dim, BaseModel<1> * dist)
    {
      OPENMS_PRECONDITION(dim < D, "ProductModel<D>:getModel(Position): index overflow!");
      if (dist == nullptr || dist == distributions_[dim])
      {
        return *this;
      }

      delete distributions_[dim];
      distributions_[dim] = dist;

      // Update model info
      String name = Peak2D::shortDimensionName(dim);
      this->param_.removeAll(name + ':');
      this->param_.insert(name + ':', distributions_[dim]->getParameters());
      this->param_.setValue(name, distributions_[dim]->getName());

      return *this;
    }

    BaseModel<1> * getModel(UInt dim) const
    {
      OPENMS_PRECONDITION(dim < D, "ProductModel<D>:getModel(Position): index overflow!");
      return distributions_[dim];
    }

    /// return the intensity scaling factor
    IntensityType getScale() const
    {
      return scale_;
    }

    /// set the intensity scaling factor
    void setScale(IntensityType scale)
    {
      this->setCutOff(this->getCutOff() / scale_);  // remove scaling from cutoff
      scale_ = scale;
      this->param_.setValue("intensity_scaling", scale);
      this->setCutOff(this->getCutOff() * scale_);  // scale cutoff
    }

    /// get reasonable set of samples from the model (i.e. for printing)
    void getSamples(SamplesType & cont) const override
    {
      cont.clear();
      typedef BaseModel<1>::SamplesType Samples1D;
      std::vector<Samples1D> samples(D);
      // get samples for each dimension
      for (Size dim = 0; dim < D; ++dim)
      {
        distributions_[dim]->getSamples(samples[dim]);
      }

      BaseModel<D>::PeakType peak;
      std::vector<UInt> i(D, 0);  // index vector

      while (i[D - 1] < samples[D - 1].size())
      {
        for (UInt dim = 0; dim < D; ++dim)
        {
          peak.getPosition()[dim] = samples[dim][i[dim]].getPosition()[0];
        }
        fillIntensity(peak);
        cont.push_back(peak);

        ++i[0];
        for (Size dim = 0; dim < D - 1; ++dim)
        {
          if (i[dim] >= samples[dim].size())
          {
            i[dim] = 0;
            ++i[dim + 1];
          }
        }
      }
      return;
    }

protected:
    void updateMembers_() override
    {
      BaseModel<D>::updateMembers_();
      scale_ = (double)(this->param_.getValue("intensity_scaling"));
      for (UInt dim = 0; dim < D; ++dim)
      {
        String name = Peak2D::shortDimensionName(dim);
        if (this->param_.exists(name))
        {
          delete distributions_[dim];
          distributions_[dim] = Factory<BaseModel<1> >::create(this->param_.getValue(name));
          Param copy = this->param_.copy(name + ":", true);
          distributions_[dim]->setParameters(copy);
          if (distributions_[dim]->getName().hasSubstring("IsotopeModel"))
          {
            static_cast<IsotopeModel *>(distributions_[dim])->setSamples(static_cast<IsotopeModel *>(distributions_[dim])->getFormula());
          }
        }
      }
    }

    std::vector<BaseModel<1> *> distributions_;
    IntensityType scale_;
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_PRODUCTMODEL_H
