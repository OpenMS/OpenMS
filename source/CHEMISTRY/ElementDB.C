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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/SYSTEM/File.h>

#include <cmath>
#include <iostream>

using namespace std;

namespace OpenMS
{
	ElementDB::ElementDB()
	{
		readFromFile_("CHEMISTRY/Elements.xml");
	}

	ElementDB::~ElementDB()
	{
		clear_();
	}

	const Map<String, const Element*>& ElementDB::getNames() const 
	{
		return names_;
	}
	
	const Map<String, const Element*>& ElementDB::getSymbols() const
	{
		return symbols_;
	}

	const Map<UInt, const Element*>& ElementDB::getAtomicNumbers() const
	{
		return atomic_numbers_;
	}

	const Element* ElementDB::getElement(const String& name) const
	{
		if (names_.has(name))
		{
			return names_[name];
		}
		else
		{
			if (symbols_.has(name))
			{
				return symbols_[name];
			}
		}
		return 0;
	}

	const Element* ElementDB::getElement(UInt atomic_number) const
	{
		if (atomic_numbers_.has(atomic_number))
		{
			return atomic_numbers_[atomic_number];
		}
		return 0;
	}
	
	bool ElementDB::hasElement(const String& name) const
	{
		return (names_.has(name) || symbols_.has(name));
	}

	bool ElementDB::hasElement(UInt atomic_number) const
	{
		return (atomic_numbers_.has(atomic_number));
	}
	
	void ElementDB::readFromFile_(const String& file_name) 
	{
		String file = File::find(file_name);

		Param param;
		param.load(file);

		UInt an(0);
		double avg_weight(0), mono_weight(0);
		String name, symbol;
		
		// build prefix
		vector<String> split;
		param.begin().getName().split(':',split);
		String prefix("");
		for (Size i=0;i<split.size()-1;++i)
		{
			prefix += split[i]+":";
		}
		
		Map<UInt, double> distribution;
		
		for (Param::ParamIterator it=param.begin(); it!=param.end(); ++it)
		{
			if (!it.getName().hasPrefix(prefix))
			{
				// update prefix
				it.getName().split(':',split);
				prefix = "";
				for (Size i=0;i<split.size()-1;++i)
				{
					prefix += split[i]+":";
				}

				IsotopeDistribution isotopes = parseIsotopeDistribution_(distribution);

				distribution.clear();

				// new element to be build
				Element* e = new Element(name, symbol, an, avg_weight, mono_weight, isotopes);
				names_[name] = e;
				symbols_[symbol] = e;
				atomic_numbers_[an] = e;

				// add all the individual isotopes as separat elements
				for (IsotopeDistribution::ConstIterator iit = isotopes.begin(); iit != isotopes.end(); ++iit)
				{
					String iso_name = "(" + String(iit->first) + ")" + name;
					String iso_symbol = "(" + String(iit->first) + ")" + symbol;
					DoubleReal iso_avg_weight = mono_weight + (iit->first - isotopes.begin()->first); // TODO exact masses
					DoubleReal iso_mono_weight = iso_avg_weight;
					IsotopeDistribution iso_isotopes;
					vector<pair<Size, double> > iso_container;
					iso_container.push_back(make_pair(iit->first, 1.0));
					iso_isotopes.set(iso_container);
					Element* iso_e = new Element(iso_name, iso_symbol, an, iso_avg_weight, iso_mono_weight, iso_isotopes);
					names_[iso_name] = iso_e;
					names_[iso_symbol] = iso_e;
				}
			}

			// read the contents of the element section
			it.getName().split(':',split);
			String key = split[2];
			String value = it->value;
			value.trim();

			if (key == "AtomicNumber")
			{
				an = (UInt)value.toInt();
			}
			else
			{
				if (key == "AverageWeight")
				{
					avg_weight = value.toDouble();
				}
				else
				{
					if (key == "Isotopes")
					{
						distribution[UInt(split[3].toInt())] = double(value.toDouble()/100);
					}
					else
					{
						if (key == "Name")
						{
							name = value;
						}
						else
						{
							if (key == "Symbol")
							{
								symbol = value;
							}
							else
							{
								if (key == "MonoWeight")
								{
									mono_weight = value.toDouble();
								}
								else
								{
									cerr << "read unknown tag: " << key << endl;
								}
							}
						}
					}
				}
			}
		}

		// build last element
		IsotopeDistribution isotopes = parseIsotopeDistribution_(distribution);
		Element * e = new Element(name, symbol, an, avg_weight, mono_weight, isotopes);
		names_[name] = e;
		symbols_[symbol] = e;
		atomic_numbers_[an] = e;
	}

	IsotopeDistribution ElementDB::parseIsotopeDistribution_(const Map<UInt, double>& distribution) 
	{
		IsotopeDistribution::ContainerType dist;
		for (Map<UInt, double>::ConstIterator it=distribution.begin(); it!=distribution.end(); ++it)
		{
			dist.push_back(make_pair(it->first, it->second));
		}
	
		IsotopeDistribution iso_dist;
		iso_dist.set(dist);
		iso_dist.setMaxIsotope(100);
	
		return iso_dist;
	}
	
	void ElementDB::clear_()
	{
		Map<String, const Element*>::Iterator it=names_.begin();
		for (;it!=names_.end();++it)
		{
			delete it->second;
		}
		names_.clear();
		symbols_.clear();
		atomic_numbers_.clear();
	}
} // namespace OpenMS

