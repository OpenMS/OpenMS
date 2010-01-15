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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------


#include <OpenMS/CHEMISTRY/ModifierRep.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <set>
#include <list>
#include <sstream>

using namespace std;

namespace OpenMS
{

ModifierRep::ModifierRep()
{
	//ResidueDB* rdb = ResidueDB::getInstance();

	//char aa[] = "ARNDCEQGHILKMFPSTWYV";

	for (Size i = 0; i<256;++i)
	{
		modification_table_.push_back(vector<double>());
	}

	number_of_modifications_= 0;

	// @todo reactivate modifications for suffix trees (andreas)
	/*
	for (Size i = 0; i<strlen(aa);++i)
	{
		const Residue * r = rdb.getResidue(aa[i]);
		set< const ResidueModification * > mods (rdb.getModifications(r));
		set< const ResidueModification * >::iterator it (mods.begin());
		for(;it!=mods.end();++it)
		{
			double add_mass = (*it)->getAddAverageWeight();
			double del_mass = (*it)->getDelAverageWeight();
			if (add_mass>0) modification_table_.at((int)aa[i]).push_back(add_mass);
			if (del_mass>0) modification_table_.at((int)aa[i]).push_back(del_mass*-1);
			//cout<<aa[i]<<" : "<<add_mass<<endl;
			//cout<<aa[i]<<" : "<<del_mass<<endl;
		}
	}*/
}

ModifierRep::ModifierRep (const ModifierRep & source):
	modification_table_(source.modification_table_),
	number_of_modifications_(source.number_of_modifications_)
{

}

ModifierRep::~ModifierRep ()
{

}

void ModifierRep::setNumberOfModifications(Size i)
{
	number_of_modifications_ = i;
}

Size ModifierRep::getNumberOfModifications() const
{
	return (number_of_modifications_);
}

const vector<vector<double> > & ModifierRep::getModificationTable ()
{
	return (modification_table_);
}

Size ModifierRep::getMaxModificationMasses ()
{
	if (number_of_modifications_==0)
	{
		return(0);
	}
	map<double, SignedSize> mod_masses;

	for (Size i = 0; i < modification_table_.size();++i)
	{
		for (Size j = 0; j < modification_table_.at(i).size();++j)
		{
			mod_masses[modification_table_.at(i).at(j)]=1;
		}
	}
	vector<double> all_single_mods;

	map<double, SignedSize>::iterator it;
	for (it = mod_masses.begin();it!=mod_masses.end();++it)
	{
		all_single_mods.push_back(it->first);
	}

	for (Size k = 1; k < number_of_modifications_; ++k)
	{
		vector<double> to_add;
		for (Size i = 0; i < all_single_mods.size();++i)
		{
			map<double, SignedSize>::iterator it2;
			for (it2 = mod_masses.begin();it2!=mod_masses.end();++it2)
			{
				to_add.push_back((it2->first)+all_single_mods.at(i));
			}


		}
		for (Size j = 0; j < to_add.size();++j)
		{
			mod_masses[to_add.at(j)]=1;
		}
	}

	return(mod_masses.size());
}

void ModifierRep::refreshModificationList (map<double, SignedSize> & mod_map,const char & c)
{
	if (modification_table_.at((int)c).size()==0) {
		return;
	} else
	{
		for (Size i = 0;i<modification_table_.at(int(c)).size();++i)
		{
			double mod_mass = modification_table_.at(int(c)).at(i);
			map<double, SignedSize>::iterator it;
			vector<pair<double, SignedSize> > to_add;
			for (it=mod_map.begin();it!=mod_map.end();++it){
				if(it->second < (SignedSize)number_of_modifications_)
				{
					to_add.push_back(pair<double, SignedSize>(it->first+mod_mass,it->second+1));
				}
			}
			for (Size j = 0; j<to_add.size();++j){
				mod_map[to_add.at(j).first] = to_add.at(j).second;
			}
			mod_map[mod_mass] = 1;
		}

	}
}

vector<String> ModifierRep::getModificationsForMass (double & m) {
	if (number_of_modifications_==0)
	{
		return(vector<String>());
	}
	// converting double to string
	stringstream ss;
	ss<<m;
	String mm = ss.str();
	/* if it does not exist the mass mapping hat to be calculated. This consits of several steps:
		1. getting all aminoacids with modifications (A,B)
		2. getting all possible combinations of this modifications with maximal number of modifications (A,B,AB,AA,BB - NOT BA because BA == AB)
		3. getting all masses for every combination
		4. joining the results by the mass and create a map <double,vector<string> >
		5. convertig double to string and deleting doubled entrys
		6. saving the mapping
	*/
	if (mass_mapping_.size()==0&&number_of_modifications_>0){
		set<String> all_mods;

		for (Size i = 0; i < modification_table_.size();++i)
		{
			for (Size j = 0; j < modification_table_.at(i).size();++j)
			{
				all_mods.insert(String ((char)i));
			}
		}

		set<String>::iterator it = all_mods.begin();
		vector<String> res;
		for (;it!=all_mods.end();++it)
		{
			res.push_back(*it);
		}

		for (Size k = 1; k < number_of_modifications_; ++k)
		{
			vector<String> to_add;
			set<String>::iterator it2 = all_mods.begin();
			SignedSize c = 0;
			for (;it2!=all_mods.end();++it2)
			{
				for (Size i = c; i < res.size();++i)
				{
					to_add.push_back(res.at(i)+(*it2));
				}
				c++;
			}
			res.insert(res.begin(),to_add.begin(),to_add.end());
		}

		map<double,vector<String> > mass_mapping;
		for (Size i = 0; i < res.size();++i)
		{
			vector<double> masses;
			for (Size j = 0 ; j < res.at(i).length();++j)
			{
				if (masses.size()==0)
				{
					masses.insert(masses.end(), modification_table_.at((int)res.at(i)[j]).begin(),modification_table_.at((int)res.at(i)[j]).end());
				} else
				{
					vector<double> masses_help;
					for (Size k = 0; k < modification_table_.at((int)res.at(i)[j]).size();k++)
					{
						for (Size l = 0 ; l < masses.size(); ++l)
						{
							masses_help.push_back(masses.at(l)+modification_table_.at((int)res.at(i)[j]).at(k));
						}
					}
					masses = masses_help;
				}
			}

			for (Size j = 0; j < masses.size();j++)
			{
				if (mass_mapping.find(masses.at(j))==mass_mapping.end())
				{
					vector<String> help;
					help.push_back(res.at(i));
					mass_mapping[masses.at(j)] = help;
				} else
				{
					mass_mapping[masses.at(j)].push_back(res.at(i));
				}
			}
		}
		map<double,vector<String> >::iterator it3 = mass_mapping.begin();
		vector<double> key_set;
		for (;it3!=mass_mapping.end();++it3)
		{
			key_set.push_back(it3->first);
		}
		map<String,vector<String> > mass_mapping_string;
		for (Size i = 0; i < key_set.size();i++)
		{
			list<String> list_help;
			list_help.insert(list_help.end(),mass_mapping[key_set.at(i)].begin(),mass_mapping[key_set.at(i)].end());
			list_help.unique();
			vector<String> unique_strings;
			unique_strings.insert(unique_strings.end(),list_help.begin(),list_help.end());
			stringstream sst ;
			sst<<key_set.at(i);
			mass_mapping_string[sst.str()]=unique_strings;
		}

		mass_mapping_ =mass_mapping_string;

		/*map<String,vector<String> >::iterator it4 = mass_mapping_.begin();
		for (;it4!=mass_mapping_.end();++it4){
			cout<<endl<<it4->first<<endl;
			for (Size i = 0; i < it4->second.size();++i){
				cout<<it4->second.at(i)<<",";
			}
		}*/
	}
	if (mass_mapping_.find(mm)==mass_mapping_.end())
	{
		return (vector<String>());
	} else
	{
		return (mass_mapping_[mm]);
	}

}


vector<String> ModifierRep::getModificationsForMass (double & m, const String & seq)
{
	vector<String> all_mods = getModificationsForMass(m);
	if (all_mods.size()==0) return (all_mods);

	vector<int> seq_hist;
	for (int i = 0; i < 256;i++)
	{
		seq_hist.push_back(0);
	}
	for (Size i = 0; i < seq.length();++i)
	{
		seq_hist[(int)seq[i]]++;
	}
	vector<String> res;

	for (Size j = 0; j < all_mods.size();j++)
	{
		vector<int> mod_hist;
		for (int i = 0; i < 256;i++)
		{
			mod_hist.push_back(0);
		}
		for (Size i = 0; i < all_mods.at(j).length();++i)
		{
			mod_hist[(int)all_mods.at(j)[i]]++;
		}
		bool to_add = true;
		for (Size i = 0; i < mod_hist.size();i++)
		{
			if (mod_hist[i]>seq_hist[i])
			{
				to_add = false;
			}
		}
		if (to_add) res.push_back(all_mods.at(j));
	}
	return res;

}

} // namespace OpenMS
