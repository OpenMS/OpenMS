// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_FEATURE_H
#define OPENMS_KERNEL_FEATURE_H

#include <OpenMS/KERNEL/RichPeak2D.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>

namespace OpenMS
{

  /**
  	@brief A feature.

  	A feature represents a subset of peaks in a map.  In general, it
  	summarizes all peaks related to a specific peptide or chemical compound
  	and thus reduces partitions of the LCMS dataset to a more meaningful
  	entity.  Picked peaks and raw data points can be converted to features
  	through the FeatureFinder.

  	Features are usually contained in FeatureMap.  Features themselves can
  	either contain features again (composite design pattern) or their
  	constituent peaks.

  	Hierarchical relationships in features (i.e. features containing features
  	containing features...) can be used to express group relationships.  For
  	example, one might group the features corresponding to an ICAT pair into
  	an aggregate ICAT feature.

  	Features are derived from DPeak, as they inherit most of their properties.
  	In particular, a feature has a position and an intensity.  The position of
  	a feature is defined as maximum position of the model for the retention
  	time dimension and the mass of the monoisotopic peak for the m/z
  	dimension.  The intensity of a feature is (proportional to) its total ion
  	count.

  	@ingroup Kernel
  */
  class Feature
    : public RichPeak2D
  {
    public:
      ///@name Type definitions
      //@{
      ///Dimensionality of the feature
      enum { DIMENSION = 2 };
      ///Type of the quality values
      typedef	 DoubleReal QualityType;
      ///Charge type
      typedef Int ChargeType;
      //@}

      /** @name Constructors and Destructor
      */
      //@{
      /// Default constructor
      inline Feature()
          : RichPeak2D(),
          overall_quality_(),
          convex_hulls_(),
          convex_hulls_modified_(true),
          convex_hull_(),
          charge_( 0 )
      {
        std::fill( qualities_, qualities_ + 2, 0 );
      }

      /// Copy constructor
      inline Feature( const Feature& feature )
          : RichPeak2D( feature ),
          overall_quality_( feature.overall_quality_ ),
          model_desc_( feature.model_desc_ ),
          convex_hulls_( feature.convex_hulls_ ),
          convex_hulls_modified_(feature.convex_hulls_modified_),
          convex_hull_( feature.convex_hull_ ),
          charge_( feature.charge_ ),
          identifications_( feature.identifications_ )
      {
        std::copy( feature.qualities_, feature.qualities_ + 2, qualities_ );
      }

      /// Destructor
      ~Feature()
      {
      }
      //@}

      ///	@name Model and Quality methods
      //@{
      /// Non-mutable access to the overall quality
      inline QualityType getOverallQuality() const
      {
        return overall_quality_;
      }
      /// Set the overall quality
      inline void setOverallQuality( QualityType q )
      {
        overall_quality_ = q;
      }

      /// Non-mutable access to the quality in dimension c
      inline QualityType getQuality( UInt index ) const
      {
        OPENMS_PRECONDITION( index < 2, "Feature<2>:getQuality(UInt): index overflow!" );
        return qualities_[ index ];
      }
      /// Set the quality in dimension c
      inline void setQuality( UInt index, QualityType q )
      {
        OPENMS_PRECONDITION( index < 2, "Feature<2>:setQuality(UInt): index overflow!" );
        qualities_[ index ] = q;
      }
			
			/// Non-mutable access to the model description
      inline const ModelDescription<2>& getModelDescription() const
      {
        return model_desc_;
      }
      /// Mutable access to the model description
      inline ModelDescription<2>& getModelDescription()
      {
        return model_desc_;
      }
      /// Set the model description
      inline void setModelDescription( const ModelDescription<2>& q )
      {
        model_desc_ = q;
      }
			//@}
			
      /// Non-mutable access to charge state
      inline const ChargeType& getCharge() const
      {
        return charge_;
      }
      /// Set charge state
      inline void setCharge( const ChargeType& ch )
      {
        charge_ = ch;
      }
			
			///@name Convex hulls and bounding box
      //@{
      /// Non-mutable access to the convex hulls
      inline const std::vector<ConvexHull2D>& getConvexHulls() const
      {
        return convex_hulls_;
      }
      /// Mutable access to the convex hulls of single mass traces
      inline std::vector<ConvexHull2D>& getConvexHulls()
      {
      	convex_hulls_modified_ = true;
        return convex_hulls_;
      }
      /// Set the convex hulls of single mass traces
      inline void setConvexHulls( const std::vector<ConvexHull2D>& hulls )
      {
      	convex_hulls_modified_ = true;
        convex_hulls_ = hulls;
      }
      /**
      	@brief Returns the overall convex hull of the feature (calculated from the convex hulls of the mass traces)
      	
      	@note the bounding box of the feature can be accessed through the returned convex hull
      */
      ConvexHull2D& getConvexHull() const;
      
      ///Returns if the mass trace convex hulls of the feature enclose the position specified by @p rt and @p mz
      bool encloses(DoubleReal rt, DoubleReal mz) const;
      //@}

      /// Assignment operator
      Feature& operator = ( const Feature& rhs );

      /// Equality operator
      bool operator == ( const Feature& rhs ) const;

      /// Compare by getOverallQuality()
      struct OverallQualityLess
            : std::binary_function < Feature, Feature, bool >
      {
        inline bool operator () ( Feature const & left, Feature const & right ) const
        {
          return ( left.getOverallQuality() < right.getOverallQuality() );
        }
        inline bool operator () ( Feature const & left, QualityType right ) const
        {
          return ( left.getOverallQuality() < right );
        }
        inline bool operator () ( QualityType left, Feature const & right ) const
        {
          return ( left < right.getOverallQuality() );
        }
        inline bool operator () ( QualityType left, QualityType right ) const
        {
          return ( left < right );
        }
      };

      /// returns a const reference to the PeptideIdentification vector
      inline const std::vector<PeptideIdentification>& getPeptideIdentifications() const
      {
        return identifications_;
      };

      /// returns a mutable reference to the PeptideIdentification vector
      inline std::vector<PeptideIdentification>& getPeptideIdentifications()
      {
        return identifications_;
      };

      /// sets the PeptideIdentification vector
      inline void setPeptideIdentifications( const std::vector<PeptideIdentification>& identifications )
      {
        identifications_ = identifications;
      };

    protected:
      /// Overall quality measure of the feature
      QualityType overall_quality_;
      /// Quality measures for each dimension
      QualityType qualities_[ 2 ];
      /// Description of the theoretical model the feature was constructed with
      ModelDescription<2> model_desc_;
      /// Array of convex hulls (one for each mass trace)
      std::vector<ConvexHull2D> convex_hulls_;
      /// Flag that indicates if the overall convex hull needs to be recomputed (i.e. mass trace convex hulls were modified)
      mutable bool convex_hulls_modified_;
      /// Overall convex hull of the feature
      mutable ConvexHull2D convex_hull_;
			/// Charge of the peptide represented by this feature.  The default value is 0, which represents an unknown charge state.  
      ChargeType charge_;
      /// Peptide PeptideIdentifications belonging to the feature
      std::vector<PeptideIdentification> identifications_;
  };

} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATURE_H
