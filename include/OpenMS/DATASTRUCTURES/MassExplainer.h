// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_MASSEXPLAINER_H
#define OPENMS_DATASTRUCTURES_MASSEXPLAINER_H

#include <iostream>
#include <algorithm>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>

namespace OpenMS
{

  /**
    @brief computes empirical formulae for given mass differences using a set of allowed elements
    
  	
  	@ingroup Datastructures
  */
  class OPENMS_DLLAPI MassExplainer
  {

	 public:
			 
		typedef Adduct::AdductsType AdductsType;
		typedef std::vector<Compomer>::const_iterator CompomerIterator;
		
    ///@name Constructors and destructor
    //@{
		/// Default constructor
		MassExplainer()
			: explanations_(),
				adduct_base_(),
				q_min_(1),
				q_max_(5),
				max_span_(3)
		{
				init_(true);
		}
		
		/// Constructor
		MassExplainer(AdductsType adduct_base)
			:	explanations_(),
				adduct_base_(adduct_base),
				q_min_(1),
				q_max_(5),
				max_span_(3)
		{
				init_(true);
		}

		/// Constructor
		MassExplainer(Int q_min, Int q_max, Int max_span, DoubleReal thresh_logp)
			:	explanations_(),
				adduct_base_(),
				q_min_(q_min),
				q_max_(q_max),
				max_span_(max_span),
				thresh_p_(thresh_logp)
		{
				init_(false);
		}

		/// Constructor
		MassExplainer(AdductsType adduct_base, Int q_min, Int q_max, Int max_span, DoubleReal thresh_logp)
			:	explanations_(),
				adduct_base_(adduct_base),
				q_min_(q_min),
				q_max_(q_max),
				max_span_(max_span),
				thresh_p_(thresh_logp)
		{
				init_(false);
		}

		
	private:
		/// check consistency of input
		/// @param init_thresh_p set default threshold (set to "false" to keep current value)
		void init_(bool init_thresh_p)
		{
				if (init_thresh_p)
				{
					// every compound with log_p_<thresh_p_ will be discarded
					// we allow at most two Na+ 
					thresh_p_ = log(0.15)*2 + log(0.7)*(q_max_-2);
				}
				
				// check consistency of members
				if (q_max_ < q_min_)
				{
						Int tmp = q_max_;
						q_max_ = q_min_;
						q_min_ = tmp;
						std::cerr << "__FILE__: Warning! \"q_max < q_min\" needed fixing!\n";
				}

				if (max_span_ > (q_max_ - q_min_ + 1))
				{
						max_span_ = q_max_ - q_min_ + 1;
						std::cerr << "__FILE__: Warning! \"max_span_ > (q_max - q_min + 1)\" needed fixing!\n";
				}
				
				if (adduct_base_.size() == 0)
				{
					//default adducts are: H+, Na+, K+, NH4+
					// do NOT use "+" in empirical formula, as every + will add a proton weight!
					
					adduct_base_.push_back(createAdduct_("H",1, 0.7));
					adduct_base_.push_back(createAdduct_("Na",1, 0.1));
					adduct_base_.push_back(createAdduct_("NH4",1, 0.1));
					adduct_base_.push_back(createAdduct_("K",1, 0.1));
					
				}
		}
		
	public:
		/// Assignment operator
		MassExplainer& operator = (const MassExplainer& rhs)
		{
				if (this == &rhs) return *this;

				explanations_ = rhs.explanations_;
				adduct_base_ = rhs.adduct_base_;
				q_min_ = rhs.q_min_;
				q_max_ = rhs.q_max_;
				max_span_ = rhs.max_span_;
				thresh_p_ = rhs.thresh_p_;

				return *this;
		}
		
		/// Destructor
		virtual ~MassExplainer()
		{
		}
    //@}
    
		
		/// fill map with possible mass-differences along with their explanation
		void compute()
		{
				// calculate some initial boundaries that can be used to shorten the enumeration process
				//Int q_comp_min = -max_span_; //minimal expected charge of compomer
				//Int q_comp_max = max_span_;  //maximal expected charge of compomer
				Int max_pq = q_max_;				 //maximal number of positve adduct-charges for a compomer
				//Int max_nq = q_max_;				 //maximal number of negative adduct-charges for a compomer
				
				for (AdductsType::const_iterator it=adduct_base_.begin(); it!=adduct_base_.end(); ++it)
				{
					std::vector <Adduct> new_adducts; 
					//create new compomers
					Int i=1;
					//warning: the following code assumes that max_nq == max_pq!!
					while (abs(i*it->getCharge()) <= max_pq)
					{
							Adduct a(*it);
							// positive amount
							a.setAmount(i);
							// this might not be a valid compomer (e.g. due to net_charge excess)
							// ... but when combined with other adducts it might become feasible again
							new_adducts.push_back(a);
							++i;
					}
					
					// combine all new compomers with existing compomers
					std::vector <Adduct>::const_iterator new_it;
					std::vector <Adduct>::const_iterator new_begin = new_adducts.begin();
					std::vector <Adduct>::const_iterator new_end   = new_adducts.end();
					
					std::size_t idx_last=explanations_.size();
					for (size_t ci=0; ci < idx_last; ++ci)
					{
							for (new_it = new_begin; new_it!=new_end; ++new_it)
							{
								Compomer cmpl(explanations_[ci]);
								cmpl.add(*new_it,Compomer::LEFT);
								explanations_.push_back(cmpl);
								
								Compomer cmpr(explanations_[ci]);
								cmpr.add(*new_it,Compomer::RIGHT);
								explanations_.push_back(cmpr);
							}
					}
					// finally add new compomers to the list itself
					for (new_it = new_begin; new_it!=new_end; ++new_it)
					{
						Compomer cmpl;
						cmpl.add(*new_it,Compomer::LEFT);
						explanations_.push_back(cmpl);
						
						Compomer cmpr;
						cmpr.add(*new_it,Compomer::RIGHT);
						explanations_.push_back(cmpr);
					}
					
					//std::cout << "valid explanations: " << explanations_.size() << " after " << it->getFormula() << std::endl;
					
				} // END adduct add
				
				
				std::vector< Compomer > valids_only;
				for (size_t ci=0; ci < explanations_.size(); ++ci)
				{
					if (compomerValid_(explanations_[ci])) valids_only.push_back(explanations_[ci]);
				}
				explanations_.swap(valids_only);
				
				// sort according to (in-order) net-charge, mass, probability
				std::sort(explanations_.begin(), explanations_.end());
				
				// set Id´s of compomers, which allows to uniquely identify them (for later lookup)
				for (size_t i=0;i<explanations_.size(); ++i) explanations_[i].setID(i);
				
				//#if DEBUG_FD
				for (size_t ci=0; ci < explanations_.size(); ++ci)
				{
				//		std::cout << explanations_[ci] << " ";
				}
				//#endif
		}
		
		
    //@name Accessors
    //@{
		
		/// Sets the set of possible adducts
		void setAdductBase(AdductsType adduct_base)
		{
			adduct_base_ = adduct_base;
		}
		/// Returns the set of adducts
		AdductsType getAdductBase() const
		{
			return adduct_base_;
		}
		
		/// return a compomer by its Id (useful after a query() ).
		const Compomer& getCompomerById(Size id) const
		{
				return explanations_[id];
		}
		//@}
				
		
		/// search the mass database for explanations
		/// @param thresh_log_p minimal log probability required
		/// @return iterator range with candidates according to net_charge and mass
		SignedSize query(const Int net_charge, 
						  const float mass_to_explain, 
							const float mass_delta,
							const float thresh_log_p,
							std::vector< Compomer >::const_iterator & firstExplanation, 
							std::vector< Compomer >::const_iterator & lastExplanation) const
		{
			#ifdef DEBUG_FD
				if (fabs(mass_to_explain) < 120.0)
				{
					std::cout << "query: qnet=" << net_charge << "; explain_mass=" << mass_to_explain << "; delta=+-" << mass_delta << "\n";
				}
			#endif
	
			Compomer cmp_low(net_charge, mass_to_explain-fabs(mass_delta), 1);
			firstExplanation = lower_bound (explanations_.begin(), explanations_.end(), cmp_low);
			
			Compomer cmp_high(net_charge, mass_to_explain+fabs(mass_delta), thresh_log_p);
			lastExplanation  = lower_bound (explanations_.begin(), explanations_.end(), cmp_high);

			return std::distance(firstExplanation, lastExplanation);
		}

	 protected:
    
	  ///check if the generated compomer is valid jugded by its probability, charges etc
		bool compomerValid_(const Compomer& cmp)
		{
				// probability ok?
				if (cmp.getLogP() < thresh_p_) return false;
				
				// limit the net charge by the maximal span of co-features
				if (abs(cmp.getNetCharge()) >= max_span_) return false;
				
				if (cmp.getNegativeCharges() > q_max_) return false;
				if (cmp.getPositiveCharges() > q_max_) return false;
				
				//TODO include mass?
				//if (abs(cmp.mass_) > mass_max_) return false;
				
				//std::cout << "valid:: " << cmp <<"\n";
				return true;
		}

		/// create a proper adduct from formula and charge and probability
		Adduct createAdduct_(const String & formula, const Int charge, const DoubleReal p) const
		{

			EmpiricalFormula ef(formula);
			//effectively substract charge electron masses: (-H plus one Proton)*charge
			ef -= ("H" + String(charge)); // substracts x hydrogen
			ef.setCharge(charge); // adds x protons
			
			Adduct a(charge, 1, ef.getMonoWeight(), formula, log(p));
			
			return a;
		}
			 
		/// store possible explanations (as formula) for a certain ChargeDifference and MassDifference
		std::vector< Compomer > explanations_;	 
		/// all allowed adducts, whose combination explains the mass difference
		AdductsType adduct_base_;
		/// minimal expected charge
		Int q_min_;
		/// maximal expected charge
		Int q_max_;
		/// maximal span (in terms of charge) for co-features, e.g. a cluster with q={3,6} has span=4
		Int max_span_;
		/// minimum required probability of a compound (all other compounds are discarded)
		DoubleReal thresh_p_;
  };



			
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_MASSEXPLAINER_H

 
