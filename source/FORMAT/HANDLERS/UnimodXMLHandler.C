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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/HANDLERS/UnimodXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{
  
  UnimodXMLHandler::UnimodXMLHandler(vector<ResidueModification*>& mods, const String& filename)
		: XMLHandler(filename, "2.0"),
			avge_mass_(0.0),
			mono_mass_(0.0),
			modification_(0),
			modifications_(mods)
			
  {
  }
   
  UnimodXMLHandler::~UnimodXMLHandler()
  {
    
  }
 
  void UnimodXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	{

		tag_ = String(sm_.convert(qname));
		
		// new modification?
		if (tag_ == "umod:mod" || tag_ == "mod")
		{
			sites_.clear();
			modification_ = new ResidueModification();
			String title(attributeAsString_(attributes, "title"));
			modification_->setId(title);

			String full_name(attributeAsString_(attributes, "full_name"));
			modification_->setFullName(full_name);

			String record_id(attributeAsString_(attributes, "record_id"));
			modification_->setUniModAccession("UniMod:" + record_id);
			return;
		}

		// which residues are allowed?
		if (tag_ == "umod:specificity" || tag_ == "specificity")
		{
			// classification of mod
			String classification(attributeAsString_(attributes, "classification"));
			modification_->setSourceClassification(classification);

			// allowed site
			String site(attributeAsString_(attributes, "site"));
			sites_.push_back(site);

			// allowed positions
			ResidueModification::Term_Specificity position = ResidueModification::ANYWHERE;
			String pos(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("position")))));
			if (pos == "Anywhere")
			{
				position = ResidueModification::ANYWHERE;
			}
			else
			{
				if (pos == "Protein N-term")
				{
					position = ResidueModification::PROTEIN_N_TERM;
				}
				else
				{
					if (pos == "Protein C-term")
					{
						position = ResidueModification::PROTEIN_C_TERM;
					}
					else
					{
						if (pos == "Any C-term")
						{
							position = ResidueModification::C_TERM;
						}
						else
						{
							if (pos == "Any N-term")
							{
								position = ResidueModification::N_TERM;
							}
							else
							{
								warning(LOAD, String("Don't know allowed position called: '") + pos  + "' - setting to anywhere");
							}
						}
					}
				}
			}
			modification_->setTermSpecificity(position);

			return;
		}
	

		if (tag_ == "umod:NeutralLoss" || tag_ == "NeutralLoss")
		{
			// mono_mass="97.976896" avge_mass="97.9952" flag="false"
			//                               composition="H(3) O(4) P">

		}
		
		// delta mass defintions?
		if (tag_ == "umod:delta" || tag_ == "delta")
		{
			// avge_mass="-0.9848" mono_mass="-0.984016" composition="H N O(-1)" >
			avge_mass_ = String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("avge_mass"))))).toDouble();
			mono_mass_ = String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("mono_mass"))))).toDouble();
			return;
		}

		// <umod:element symbol="H" number="1"/>
		if (tag_ == "umod:element")
		{
			String symbol = sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("symbol"))));
			String num = sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("number"))));
			String isotope, tmp_symbol;
      for (Size i = 0; i != symbol.size(); ++i)
      {
      	if (isdigit(symbol[i]))
        {
        	isotope += symbol[i];
        }
        else
        {
          tmp_symbol += symbol[i];
      	}
			}
      
      String formula;
      if (isotope != "")
      {
        formula = '(' + isotope + ')' + tmp_symbol + String(num);
      }
      else
      {
        formula = tmp_symbol + num;
      }
      diff_formula_ += formula;			
			
		}
		
	}
	  
  void UnimodXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
		tag_ = String(sm_.convert(qname));

		// write the modifications to vector
		if (tag_ == "umod:mod" || tag_ == "mod")
		{
			modification_->setAverageMass(avge_mass_);
			modification_->setMonoMass(mono_mass_);
			modification_->setDiffFormula(diff_formula_);

			for (vector<String>::const_iterator it = sites_.begin(); it != sites_.end(); ++it)
			{
				ResidueModification* new_mod = new ResidueModification(*modification_);
				new_mod->setOrigin(*it);
				modifications_.push_back(new_mod);
			}
			
			avge_mass_ = 0.0;
			mono_mass_ = 0.0;
			diff_formula_ = EmpiricalFormula();

			delete modification_;
			return;
		}

		if (tag_ == "umod:NeutralLoss" || tag_ == "NeutralLoss")
		{
			// now diff_formula_ contains the neutral loss diff formula
			modification_->setNeutralLossDiffFormula(diff_formula_);
			modification_->setNeutralLossMonoMass(mono_mass_);
			modification_->setNeutralLossAverageMass(avge_mass_);
			avge_mass_ = 0.0;
			mono_mass_ = 0.0;
			diff_formula_ = EmpiricalFormula();
		}
 	} 

  void UnimodXMLHandler::characters(const XMLCh* const /*chars*/, const XMLSize_t /*length*/)
  {
		// nothing to do here
	}

	} // namespace Internal
} // namespace OpenMS
