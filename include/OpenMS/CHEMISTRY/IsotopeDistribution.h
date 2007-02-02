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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_H

#include <OpenMS/CONCEPT/Types.h>

#include <utility>
#include <vector>

namespace OpenMS
{
  /**
		@ingroup Chemistry
	
		@brief Isotope distribution class

		Holds an isotope distribution with the weight value and according 
		probability. Distribution can be add using the '+' or '+=' operators.
	*/	
  class IsotopeDistribution
  {
		public:

			/// @name typedefs
			//@{
			/// container type, first holds the weight of the isotope, second the probability
			typedef std::vector<std::pair<UnsignedInt, double> > ContainerType;
			typedef ContainerType::iterator iterator;
			typedef ContainerType::iterator Iterator;
			typedef ContainerType::const_iterator const_iterator;
			typedef ContainerType::const_iterator ConstIterator;
			//@}

			/// @name Constructors and Destructors
			//@{
			/** Default constructor, note max_isotope must be set later 
					@see setMaxIsotope(Size max_isotope)
			*/
			IsotopeDistribution();

			/// Detailed constructor which sets the @p max_isotope
			IsotopeDistribution(Size max_isotope);

			/// Copy constructor
			IsotopeDistribution(const IsotopeDistribution& isotope_distribution);

			/// Destructor
			virtual ~IsotopeDistribution();
			//@}
	
			/// @name Accessors
			//@{
			/** @brief sets the maximal isotope with @p max_isotope

					sets the maximal isotope which is included in the distribution 
					and used to limit the calculations. This is useful as distributions
					with numerous isotopes tend to have a lot of numerical zeros at the end
			*/
			void setMaxIsotope(Size max_isotope);
			
			/// returns the currently set maximum isotope
			Size getMaxIsotope() const;
			
			/// Retrieve threshold for isotope probabilities 
			double getTrimRightCutoff()  { return trim_right_cutoff_; }
			
			/// Set threshold for isotope probabilities
			void setTrimRightCutoff(double const cutoff) { trim_right_cutoff_ = cutoff; }

			/// overwrites the container which holds the distribution using @p distribution
			void set(const ContainerType& distribution);

			/// returns the container which holds the distribtion
			const ContainerType& getContainer() const;

			/// returns the weight maximal weight isotope which is stored in the distribution
			UnsignedInt getMax() const;

			/// returns the weight minimal weight isotope which is stored in the distribution
			UnsignedInt getMin() const;

			/// returns the size of the distribtion which is the number of isotopes in the distribution
			Size size() const;

			/// clears the distribution and resets max isotope to 0
			void clear();

			/// Estimate Peptide Isotopdistribution from weight and number of isotopes that should be reported (very fast)
			void estimateFromPeptideWeight(double weight);
			//@}
			
			/// @name Operators
			//@{
			/// Assignment operator
			IsotopeDistribution& operator = (const IsotopeDistribution& isotope_distribution);
			
			/// operator which adds this distribution and the @p isotope_distribution to return IsotopeDisribution (similiar to convolve distributions)
			IsotopeDistribution operator + (const IsotopeDistribution& isotope_distribution) const;

			/// operator which adds @p isotope_distribution to this (similar to convolve distributions)
			IsotopeDistribution& operator += (const IsotopeDistribution& isotope_distribution);

			/// operator which multiplies this distribution by @p factor (similar to @p factor times applying operator '+')
			IsotopeDistribution operator * (UnsignedInt factor) const;

			/// operator which multiplies this distribution by @p factor (similar to @p factor times applying operator '+=')
			IsotopeDistribution& operator *= (UnsignedInt factor);

			/** @brief renormalizes the sum of the probabilities of the isotopes to 1

					The renormalisation is needed as in distributions with a lot of isotopes (and with high max isotope)
					the calculations tend to be inexact.
			*/
	    void renormalize();
	    
	    /** @brief Trims the isotope distribution using the member trim_right_cutoff_ .
	    
	    		If the isotope distribution is calculated for large masses (and with high max isotope)
	    		it might happen that many entries contain only small numbers. This function can be
	    		used to remove these entries.
	    		Do consider to normalize the distribution afterwards.	    
	    */
	    void trimRight(); 
			
			/// equality operator, returns true if the @p isotope_distribution is identical to this, false else
			bool operator == (const IsotopeDistribution& isotope_distribution) const;

			/// inequality operator, returns true if the @p isotope_distribution differs from this, false else
			bool operator != (const IsotopeDistribution& isotope_distribution) const;
			//@}

			/// @name Iterators
			//@{
			inline Iterator begin() { return distribution_.begin(); }

			inline Iterator end() 	{ return distribution_.end(); }

			inline ConstIterator begin() const { return distribution_.begin(); }

			inline ConstIterator end() const { return distribution_.end(); }
			//@}

		protected:

			/// convolves the distributions @p left and @p right and stores the result in @p result
			void convolve_(ContainerType& result, const ContainerType& left, const ContainerType& right) const;

			/// convolves the distribution @p input @p factor times and stores the result in @p result
			void convolvePow_(ContainerType& result, const ContainerType& input, Size factor) const;

			/// convolves the distribution @p input with itself and stores the result in @p result
			void convolveSquare_(ContainerType& result, const ContainerType& input) const;

			/// maximal isotopes which is used to calculate the distribution
			Size max_isotope_;

			/// stores the isotope distribution
			ContainerType distribution_;
			
			/// cut of for the isotope probabilities, can be used to trim the distribution
			double trim_right_cutoff_;
  };

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_H
