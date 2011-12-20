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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>


using namespace xercesc;
using namespace std;

namespace OpenMS 
{

	MascotXMLFile::MascotXMLFile()
		: Internal::XMLFile()
	{
	}

  void MascotXMLFile::load(const String& filename, 
						      					ProteinIdentification& protein_identification, 
						      					vector<PeptideIdentification>& id_data)
  {
  	map<String, vector<AASequence> > peptides;
  	
  	load(filename, protein_identification, id_data, peptides);      
  }  					 
  					 
  void MascotXMLFile::load(const String& filename, 
						      					ProteinIdentification& protein_identification, 
						      					vector<PeptideIdentification>& id_data,
						      					map<String, vector<AASequence> >& peptides)
  {
  	//clear
		protein_identification = ProteinIdentification();
		id_data.clear();

		Internal::MascotXMLHandler handler(protein_identification, id_data, filename, peptides);
		parse_(filename, &handler);
				
		// Since the mascot xml can contain "peptides" without sequences the identifications 
		// without any real peptide hit are removed
  	vector<PeptideHit> peptide_hits;
  	vector<PeptideIdentification>::iterator id_it = id_data.begin();

		while(id_it != id_data.end())
		{
			peptide_hits = id_it->getHits();
			if (peptide_hits.empty() || (peptide_hits.size() == 1 && peptide_hits[0].getSequence() == ""))
			{
				id_it = id_data.erase(id_it);
			}
			else
			{
				++id_it;
			}
		}

		// argh!
		// since Mascot xml 2.2 tends to repeat the first hit (yes it apears twice, we delete one of them)
		for (vector<PeptideIdentification>::iterator it = id_data.begin(); it != id_data.end(); ++it)
		{
			peptide_hits = it->getHits();
			// check if equal, except for rank
			if (peptide_hits.size() > 1 &&
					peptide_hits[0].getScore() == peptide_hits[1].getScore() &&
					peptide_hits[0].getSequence() == peptide_hits[1].getSequence() &&
					peptide_hits[0].getCharge() == peptide_hits[1].getCharge() /* &&
				peptide_hits[0].getProteinAccessions() == peptide_hits[1].getProteinAccessions() &&
					peptide_hits[0].getAABefore() == peptide_hits[1].getAABefore() &&
					peptide_hits[0].getAAAfter() == peptide_hits[1].getAAAfter()*/)
			{
				// erase first hit
				peptide_hits.erase(peptide_hits.begin() + 1);
				it->setHits(peptide_hits);
			}
		}
  }  					 
  					 
} // namespace OpenMS
