// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_STABLEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_STABLEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

namespace OpenMS
{
  /**
   @brief This class implements an element pair finding algorithm.

   It offers a method to determine element pairs across two element maps.
   The corresponding features must be aligned, but may have small position deviations.

   Motivation and definition of the distance measure:

   The similarity value should express our confidence that one element might
   possibly be matched to the other.  Larger quality values are better, the
   maximal similarity is one.  Let \f$\Delta_\textit{RT}\f$ and
   \f$\Delta_\textit{MZ}\f$ be the absolute values of the RT and MZ differences
   in the data (for m/z this might be given in ppm or Da - see parameters).
   Then the distance value is
   \f[
   \big(
   \frac
     {\Delta_\textit{RT}}
     {\big( \textit{max\_pair\_distance:RT} \big)}
   ^\textit{diff\_exponent\_RT}
   +
   \frac
     {\Delta_\textit{MZ}}
     {\big( \textit{max\_pair\_distance:MZ} \big)}
   ^\textit{diff\_exponent\_MZ}
   \big)

   \cdot

   \textit{intensity\_ratio}^\textit{intensity\_exponent}

   \cdot

   \textit{different\_charge\_penalty}

   \f]

   Choosing <i>diff_exponent</i>: This parameter controls the growth rate of the
   penalty for differences.  It is for example possible to search for pairs
   using the absolute distance in RT (which should not be very susceptible to
   outliers) and the squared distance in MZ (where small difference occur
   frequently, but large differences indicate a mismatch).

   Stability criterion:  The similarity of the second nearest neighbor must be
   smaller than the similarity of the nearest neighbor by a certain factor, see
   <i>second_nearest_gap</i>.

   @htmlinclude OpenMS_StablePairFinder.parameters

   @ingroup FeatureGrouping
   */
  class OPENMS_DLLAPI StablePairFinder : public BaseGroupFinder
  {
    public:

      ///Base class
      typedef BaseGroupFinder Base;

      /// Constructor
      StablePairFinder();

      /// Destructor
      virtual
      ~StablePairFinder()
      {
      }

      /// Returns an instance of this class
      static BaseGroupFinder*
      create()
      {
        return new StablePairFinder();
      }

      /// Returns the name of this module
      static const String
      getProductName()
      {
        return "stable";
      }

      /**
       @brief Run the algorithm

       @note Exactly two @em input maps must be provided.

       @exception Exception::IllegalArgument is thrown if the input data is not valid.
       */
      void
      run( const std::vector<ConsensusMap>& input_maps,
           ConsensusMap &result_map );

    protected:

      ///@name Internal helper classes and enums
      //@{
      enum
      {
        MODEL_ = 0,
        SCENE_ = 1
      };
      enum
      {
        RT = Peak2D::RT,
        MZ = Peak2D::MZ
      };
      //@}

      //docu in base class
      virtual void
      updateMembers_();

      /// Computes the distance for a pair of elements.
      DoubleReal distance_( ConsensusFeature const & left,
														ConsensusFeature const & right ) const;

			/**
				 @brief Checks if the peptide IDs of two features are compatible.

				 A feature without identification is always compatible. Otherwise, two features are compatible if the best peptide hits of their identifications have the same sequence.
			*/
			bool compatibleIDs_(const ConsensusFeature& feat1, 
													const ConsensusFeature& feat2) const;

      /// Distances wrt RT and MZ are raised to this power, respectively.
      DoubleReal diff_exponent_[2];

      /// Ratio of intensities is raised to this power.
      DoubleReal intensity_exponent_;

      /// Maximal distance of a matched pair, in both dimension RT and MZ.
      DoubleReal max_pair_distance_[2];

      /// Unit of distance between two m/z values (Da or ppm)
      bool max_pair_distance_mz_as_Da_;

      /// Reciprocal values of #max_pair_distance_ (which see).  ('*' is faster than '/'.)
      DoubleReal max_pair_distance_reciprocal_[2];

      /// The distance of the second nearest neighbors must be this factor larger than for the pair itself.
      DoubleReal second_nearest_gap_;

      /// If charge states are different, distance is multiplied by this factor.
      DoubleReal different_charge_penalty_;

  };

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_STABLEPAIRFINDER_H

/*

 gnuplot history - how the plot was created - please do not delete this receipt

 f(x,intercept,exponent)=1/(1+(abs(x)*intercept)**exponent)
 set terminal postscript enhanced color
 set output "choosingstablepairfinderparams.ps"
 set size ratio .3
 plot [-3:3] [0:1] f(x,1,1), f(x,2,1), f(x,1,2), f(x,2,2)

 */
