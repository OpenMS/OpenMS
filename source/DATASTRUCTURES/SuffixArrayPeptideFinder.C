// -*- Mode: C++; tab-width: 2; -*- vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bauer$
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/SuffixArrayPeptideFinder.h>
#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/CHEMISTRY/ModifierRep.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/DSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/DATASTRUCTURES/SuffixArraySeqan.h>
#include <OpenMS/DATASTRUCTURES/SuffixArrayTrypticSeqan.h>
#include <OpenMS/DATASTRUCTURES/SuffixArrayTrypticCompressed.h>
#include <OpenMS/DATASTRUCTURES/SuffixArray.h>
#include <fstream>

using namespace std;

namespace OpenMS
{
				
SuffixArrayPeptideFinder::SuffixArrayPeptideFinder(const String& f_file, const String& method) throw (Exception::FileNotFound,Exception::ParseError,Exception::InvalidValue)
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

	fstream fs;
	String saFileName = (f_file.substr(0, f_file.length() - 6)); // @todo dangerous remove only suffix! (Chris, Andreas)
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
	
	cout << "file:" << saFileName << endl;
	if (method == "trypticCompressed")
	{
		sa_ = new SuffixArrayTrypticCompressed(big_string_.getBigString(),saFileName);
	} 
	else
	{
		if (method == "seqan") 
		{
			sa_ = new SuffixArraySeqan(big_string_.getBigString(),saFileName);
		} 
		else 
		{
			if (method == "trypticSeqan") 
			{
				sa_ = new SuffixArrayTrypticSeqan(big_string_.getBigString(),saFileName);
			}
		}
	}
	
	cout << "done" << endl;

	if (saFileName=="")
	{
		cout << "saving" << endl;
		sa_->save(saFileNameCopy);
	}
}

SuffixArrayPeptideFinder::SuffixArrayPeptideFinder(const SuffixArrayPeptideFinder & source)
{
	sa_=source.sa_;
	big_string_=source.big_string_;
}

SuffixArrayPeptideFinder::~SuffixArrayPeptideFinder()
{

}

void SuffixArrayPeptideFinder::setTolerance(const float t){
	sa_->setTolerance(t);
}

float SuffixArrayPeptideFinder::getTolerance() const
{
	return (sa_->getTolerance());
}

void SuffixArrayPeptideFinder::setNumberOfModifications(UInt number_of_mods) const
{
	sa_->setNumberOfModifications(number_of_mods);
}

UInt SuffixArrayPeptideFinder::getNumberOfModifications() const
{
	return (sa_->getNumberOfModifications());
}

void SuffixArrayPeptideFinder::setTags(const vector<OpenMS::String> & tags) throw (OpenMS::Exception::InvalidValue)
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

void SuffixArrayPeptideFinder::setModificationOutputMethod (const String & s) throw (OpenMS::Exception::InvalidValue)
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
	for (UInt i = 0; i < v.size();++i)
	{
		res += v.at(i);
		if (i < v.size() - 1)
		{
			res += ",";
		}
	}
	res += "]";
	return (res);
}

vector<vector<pair<SuffixArrayPeptideFinder::FASTAEntry , String> > > SuffixArrayPeptideFinder::getCandidates(const vector<double> & spec)
{
	vector<vector<pair<pair<int,int>,float> > > ca = sa_->findSpec(spec);
	ModifierRep* mod = new ModifierRep ();
	mod->setNumberOfModifications(sa_->getNumberOfModifications());
	vector<vector<pair<FASTAEntry,String> > > res;
	for (UInt i = 0; i < ca.size(); i++)
	{
		vector<pair<FASTAEntry,String> > temp;
		for (UInt j = 0; j < ca.at(i).size(); j++)
		{
			FASTAEntry fe (big_string_.getPeptide(ca.at(i).at(j).first.first,ca.at(i).at(j).first.second));
			String mm = "";
			if (modification_output_method_=="mass")
			{
				stringstream ss;
				ss << ca.at(i).at(j).second;
				mm = ss.str();
			} 
			else 
			{
				double ma = (double)ca.at(i).at(j).second;
				if (modification_output_method_ == "stringUnchecked")
				{
					mm = vToString_(mod->getModificationsForMass (ma));
				} 
				else
				{
					if (modification_output_method_ == "stringChecked")
					{
						mm = vToString_(mod->getModificationsForMass (ma,fe.second));
					}
				}
			}
			String mod = (mm=="0") ? "" : mm;
			temp.push_back(pair<FASTAEntry,String>(fe, mod));
		}

		res.push_back (temp);
	}
	return res;
}

vector<vector<pair<SuffixArrayPeptideFinder::FASTAEntry, String > > > SuffixArrayPeptideFinder::getCandidates (const String & DTA_file) throw (Exception::FileNotFound,Exception::ParseError)
{
	DTAFile dta_file;
	DSpectrum<> s;
	dta_file.load(DTA_file,s);
	s.getContainer().sortByPosition();
	DSpectrum<>::ConstIterator it(s.begin());
	vector<double> spec;
	for (;it!=s.end();++it)
	{
		spec.push_back(it->getPosition()[0]);
	}
	const vector<double> specc(spec);
	return (getCandidates(specc));
}

} // namespace OpenMS
