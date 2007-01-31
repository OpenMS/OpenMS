// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PRODUCTMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PRODUCTMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/KERNEL/DimensionDescription.h>

namespace OpenMS
{
/** @brief Class for product models i.e. models with D independent dimensions
  
	The predicted intensity is simply the product of the intensities in each dimension

	Parameters:
	<table>
	<tr><td>intensity_scaling</td>
			<td>factor used to scale the calculated intensities</td></tr>
	<tr><td>cutoff</td>
			<td>peak with intensity below cutoff is not considered
					 to be part of the model</td></tr>
	<tr><td>DimensionTags: e.g. RT, MZ</td>
			<td>model used for the specified dimension including model parameters</td></tr>
	</table>
	
	@ingroup FeatureFinder
*/
template <UnsignedInt D, typename Traits = KernelTraits,
typename DimensionTags = LCMS_Tag>
class ProductModel
            : public BaseModel<D,Traits>
{

public:
    typedef typename DPeak<D,Traits>::IntensityType IntensityType;
    typedef DPosition<D,Traits> PositionType;
    typedef DPeakArray<D, DPeak<D,Traits> > SamplesType;

    /// Default constructor
    ProductModel()
	    : BaseModel<D,Traits>(),
	    	distributions_(D,0)
    {
    	this->setName(this->getProductName());
    	
    	//Register model info
      for (UnsignedInt dim=0; dim<D; ++dim)
      {
      	String name = DimensionDescription<DimensionTags>::dimension_name_short[dim];
    		this->subsections_.push_back(name);
    		this->defaults_.setValue(name,"GaussModel");
    	}
    	
    	//defaults
      this->defaults_.setValue("intensity_scaling",1.0);
      this->defaultsToParam_();
    }

    /// copy constructor
    ProductModel(const ProductModel& source)
	    : BaseModel<D,Traits>(source),
	    	distributions_(D,0),
	    	scale_(source.scale_)
    {
      for (UnsignedInt dim=0; dim<D; ++dim)
      {
        // clone source model
        if (source.distributions_[dim])
        {
          ModelDescription<1> desc(source.distributions_[dim]);
          setModel(dim,desc.createModel());
        }
      }
      updateMembers_();
    }

    /// destructor
    virtual ~ProductModel()
    {
      for (UnsignedInt dim=0; dim<D; ++dim)
      {
      	delete distributions_[dim];
  		}
    }

    /// assignment operator
    virtual ProductModel& operator = (const ProductModel& source)
    {
    	if (&source == this) return *this;
    	
        BaseModel<D,Traits>::operator = (source);
				scale_ = source.scale_;
				
        for (UnsignedInt dim=0; dim<D; ++dim)
        {
          if (source.distributions_[dim])
          {
            // clone source model
            ModelDescription<1> desc(source.distributions_[dim]);
            setModel(dim,desc.createModel());
          }
      	}
        updateMembers_();
        
        return *this;
    }

    /// intensity equals product of intensities in each dimension
    IntensityType getIntensity(const PositionType& pos) const
    {
      IntensityType intens(scale_);
      for (UnsignedInt dim=0; dim<D; ++dim)
      {
        if (distributions_[dim]==0)
        {
        	throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("ProductModel: model for dimension ") + dim + " not set.","");
				}
        intens *= distributions_[dim]->getIntensity(pos[dim]);
      }
      return intens;
    }

    /// create new ProductModel object (needed by Factory)
    static BaseModel<D,Traits>* create()
    {
        return new ProductModel<D,Traits>();
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
    ProductModel& setModel(const Position dim, BaseModel<1,Traits>* dist)
    {
      OPENMS_PRECONDITION(dim<D, "ProductModel<D>:getModel(Position): index overflow!")
      if (dist==0 || dist==distributions_[dim])
      {
      	return *this;
			}
			
      delete distributions_[dim];
      distributions_[dim] = dist;

			// Update model info
	    String name = DimensionDescription<DimensionTags>::dimension_name_short[dim];
	    this->param_.remove(name.c_str());
	    this->param_.insert(name.c_str(),distributions_[dim]->getParameters());
	    this->param_.setValue(name.c_str(), distributions_[dim]->getName());

      return *this;
    }

    BaseModel<1,Traits>* getModel(const Position dim) const
    {
      OPENMS_PRECONDITION(dim<D, "ProductModel<D>:getModel(Position): index overflow!")
      return distributions_[dim];
    }

    /// return the intensity scaling factor
    const IntensityType& getScale() const
    {
        return scale_;
    }

    /// set the intensity scaling factor
    void setScale(const IntensityType& scale)
    {
      this->setCutOff( this->getCutOff()/scale_ );	// remove scaling from cutoff
      scale_ = scale;
      this->param_.setValue("intensity_scaling",scale);
      this->setCutOff( this->getCutOff()*scale_ );	// scale cutoff
    }

    /// get reasonable set of samples from the model (i.e. for printing)
    void getSamples(SamplesType& cont) const
    {
      cont = SamplesType();
      typedef typename BaseModel<1,Traits>::SamplesType Samples1D;
      std::vector<Samples1D> samples(D);
      // get samples for each dimension
      for (UnsignedInt dim=0; dim<D; ++dim)
      {
      	distributions_[dim]->getSamples(samples[dim]);
      }

      typename BaseModel<D,Traits>::PeakType peak;
      std::vector<Size> i(D,0);  // index vector

      while(i[D-1]<samples[D-1].size())
      {
        for (UnsignedInt dim=0; dim<D; ++dim)
        {
          peak.getPosition()[dim] = samples[dim][ i[dim] ].getPosition()[0];
        }
        fillIntensity(peak);
        cont.push_back(peak);

        ++i[0];
        for (UnsignedInt dim=0; dim<D-1; ++dim)
        {
          if (i[dim]>=samples[dim].size())
          {
            i[dim]=0;
            ++i[dim+1];
          }
        }
      }
    }

	protected:
		void updateMembers_()
		{
			BaseModel<D,Traits>::updateMembers_();
			scale_ = (double)(this->param_.getValue("intensity_scaling"));
	    for (UnsignedInt dim=0; dim<D; ++dim)
      {
        String name = DimensionDescription<DimensionTags>::dimension_name_short[dim];
        DataValue d = this->param_.getValue(name);
        if (d!=DataValue::EMPTY)
        {
         	delete distributions_[dim];
          distributions_[dim] = Factory< BaseModel<1,Traits> >::create(d);
          distributions_[dim]->setParameters( this->param_.copy(name+":",true) );
        }
      }
		}

    std::vector< BaseModel<1,Traits>* > distributions_;
    IntensityType scale_;
};
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_PRODUCTMODEL_H
