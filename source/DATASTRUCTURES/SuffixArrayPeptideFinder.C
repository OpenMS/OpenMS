// -*- mode: C++; tab-width: 2; -*- vi: set ts=2:
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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/SuffixArrayPeptideFinder.h>
#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/CHEMISTRY/ModifierRep.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/SuffixArraySeqan.h>
#include <OpenMS/DATASTRUCTURES/SuffixArrayTrypticSeqan.h>
#include <OpenMS/DATASTRUCTURES/SuffixArrayTrypticCompressed.h>
#include <OpenMS/DATASTRUCTURES/SuffixArray.h>
#include <fstream>

using namespace std;

namespace OpenMS
{
				
SuffixArrayPeptideFinder::SuffixArrayPeptideFinder(const String& f_file, const String& method, const WeightWrapper::WEIGHTMODE weight_mode)
{
	if (!(method=="trypticCompressed" || method=="seqan" || method=="trypticSeqan"))
	{
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"method has to be trypticCompressed,seqan,trypticSeqan",method);
	}

	PepIterator& it = *Factory<PepIterator>::create("FastaIterator");
	try 
	{
		it.setFastaFile(f_file);
	} 
	catch (Exception::FileNotFound&)
	{
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,f_file);
	}
	it.begin();

	while (!it.isAtEnd())
	{
		big_string_.add(*it);
		++it;
	}
	modification_output_method_="mass";

	String saFileName = "";
	/*
	fstream fs;
	String saFileName = (f_file.substr(0, f_file.length() - 6)); /// @todo dangerous remove only suffix! (Chris, Andreas)

	
	const String saFileNameCopy = saFileName;
	
	fs.open((saFileName + ".sa2").c_str());
	if (fs.is_open())
	{
		cout << "sa will be loaded" << endl;
		fs.close();
	} 
	else 
	{
		cout << "sa has to be created! it will be saved afterwards for reuse" << endl;
		
		saFileName = "";
	}
	*/
	
	//cout << "file:" << saFileName << endl;
	if (method == "trypticCompressed")
	{
		sa_ = new SuffixArrayTrypticCompressed(big_string_.getBigString(), saFileName, weight_mode);
	} 
	else if (method == "seqan") 
	{
		sa_ = new SuffixArraySeqan(big_string_.getBigString(), saFileName, weight_mode);
	} 
	else if (method == "trypticSeqan") 
	{
		sa_ = new SuffixArrayTrypticSeqan(big_string_.getBigString(), saFileName, weight_mode);
	}
	
	/*
	cout << "done" << endl;

	if (saFileName=="")
	{
		cout << "saving" << endl;
		sa_->save(saFileNameCopy);
	}*/
}

SuffixArrayPeptideFinder::SuffixArrayPeptideFinder(const SuffixArrayPeptideFinder & source)
	: WeightWrapper(source)
{
	sa_=source.sa_;
	big_string_=source.big_string_;
}

SuffixArrayPeptideFinder::~SuffixArrayPeptideFinder()
{
	delete sa_; // TODO assignment, copy ctor
}

void SuffixArrayPeptideFinder::setTolerance(const DoubleReal t){
	sa_->setTolerance(t);
}

DoubleReal SuffixArrayPeptideFinder::getTolerance() const
{
	return (sa_->getTolerance());
}

void SuffixArrayPeptideFinder::setNumberOfModifications(Size number_of_mods) const
{
	sa_->setNumberOfModifications(number_of_mods);
}

Size SuffixArrayPeptideFinder::getNumberOfModifications() const
{
	return (sa_->getNumberOfModifications());
}

void SuffixArrayPeptideFinder::setTags(const vector<OpenMS::String> & tags)
{
	sa_->setTags(tags);
}

const vector<OpenMS::String> & SuffixArrayPeptideFinder::getTags()
{
	return (sa_->getTags());
}

void SuffixArrayPeptideFinder::setUseTags(bool use_tags)
{
	sa_->setUseTags(use_tags);
}

bool SuffixArrayPeptideFinder::getUseTags()
{
	return (sa_->getUseTags());
}

void SuffixArrayPeptideFinder::setModificationOutputMethod (const String & s)
{
	if (!(s == "mass" || s == "stringUnchecked" || s == "stringChecked"))
	{
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"modification output has to be mass,stringUnchecked,stringChecked",s);
	}
	modification_output_method_ = s;
}

String SuffixArrayPeptideFinder::getModificationOutputMethod ()
{
	return (modification_output_method_);
}

String SuffixArrayPeptideFinder::vToString_ (vector<String> v)
{
	if (v.size()==0) return ("");
	String res = "[";
	for (Size i = 0; i < v.size();++i)
	{
		res += v[i];
		if (i < v.size() - 1)
		{
			res += ",";
		}
	}
	res += "]";
	return (res);
}

void SuffixArrayPeptideFinder::getCandidates(vector<vector<pair<FASTAEntry, String> > >& candidates, const vector<DoubleReal>& spec)
{
	vector<vector<pair<pair<SignedSize, SignedSize>, DoubleReal> > > ca;
	sa_->findSpec(ca, spec);
	
	ModifierRep mod;
	mod.setNumberOfModifications(sa_->getNumberOfModifications());
	for (Size i = 0; i < ca.size(); i++)
	{
		set<String> already_used; // TODO
		// TODO
		if (i % 1000 == 0)
		{
			cerr << i << "/" << ca.size() << endl;
		}

		vector<pair<FASTAEntry, String> > temp;
		for (Size j = 0; j < ca[i].size(); j++)
		{
			FASTAEntry fe;
			big_string_.getPeptide(fe, ca[i][j].first.first, ca[i][j].first.second);

			String mod_str;
			if (ca[i][j].second != 0)
			{
				String mod_str;
				if (modification_output_method_ == "mass")
				{
					//stringstream ss;
					//ss << ca[i][j].second;
					//mm = ss.str();
					mod_str = String(ca[i][j].second);
				} 
				else 
				{
					DoubleReal ma = (DoubleReal)ca[i][j].second;
					if (modification_output_method_ == "stringUnchecked")
					{
						mod_str = vToString_(mod.getModificationsForMass(ma));
					} 
					else
					{
						if (modification_output_method_ == "stringChecked")
						{
							mod_str = vToString_(mod.getModificationsForMass(ma,fe.second));
						}
					}
				}
			}
			if (already_used.find(fe.second) == already_used.end())
			{
				temp.push_back(pair<FASTAEntry, String>(fe, mod_str));
				already_used.insert(fe.second);
			}
		}

		candidates.push_back(temp);
	}
	return;
}

void SuffixArrayPeptideFinder::getCandidates (vector<vector<pair<SuffixArrayPeptideFinder::FASTAEntry, String > > >& candidates, const String& DTA_file)
{
	DTAFile dta_file;
	PeakSpectrum s;
	dta_file.load(DTA_file, s);
	s.sortByPosition();
	PeakSpectrum::ConstIterator it(s.begin());
	vector<DoubleReal> spec;
	for (;it!=s.end();++it)
	{
		spec.push_back(it->getPosition()[0]);
	}
	const vector<DoubleReal> specc(spec);
	getCandidates(candidates, specc);
	return;
}

} // namespace OpenMS
