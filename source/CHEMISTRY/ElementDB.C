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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <cmath>
#include <iostream>

using namespace std;

namespace OpenMS
{
	ElementDB::ElementDB()
	{
		readFromFile_(OPENMS_PATH "/data/" OPENMS_CHEMISTRY_ELEMENTDB_DEFAULT_FILE);
	}

	ElementDB::~ElementDB()
	{
		clear_();
	}

	const HashMap<String, const Element*>& ElementDB::getNames() const 
	{
		return names_;
	}
	
	const HashMap<String, const Element*>& ElementDB::getSymbols() const
	{
		return symbols_;
	}

	const HashMap<Size, const Element*>& ElementDB::getAtomicNumbers() const
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

	const Element* ElementDB::getElement(Size atomic_number) const
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

	bool ElementDB::hasElement(Size atomic_number) const
	{
		return (atomic_numbers_.has(atomic_number));
	}
	
	void ElementDB::readFromFile_(const String& file_name) throw(Exception::FileNotFound, Exception::ParseError)
	{
		Param param;
		param.load(file_name);

		Size an(0);
		float avg_weight(0), mono_weight(0);
		String name, symbol;
		
		// build prefix
		vector<String> split;
		String(param.begin()->first).split(':',split);
		String prefix("");
		for (Size i=0;i<split.size()-1;++i)
		{
			prefix += split[i]+":";
		}
		
		HashMap<UnsignedInt, double> distribution;
		
		for (Param::ConstIterator it=param.begin(); it!=param.end(); ++it)
		{
			if (!String(it->first).hasPrefix(prefix))
			{
				// update prefix
				String(it->first).split(':',split);
				prefix = "";
				for (Size i=0;i<split.size()-1;++i)
				{
					prefix += split[i]+":";
				}

				IsotopeDistribution isotopes = parseIsotopeDistribution_(distribution);

				distribution.clear();

				// new element to be build
				Element * e = new Element(name, symbol, an, avg_weight, mono_weight, isotopes);
				names_[name] = e;
				symbols_[symbol] = e;
				atomic_numbers_[an] = e;
			}

			// read the contents of the element section
			String(it->first).split(':',split);
			String key = split[2];
			String value = String(it->second.toString()).trim();

			if (key == "AtomicNumber")
			{
				an = (Size)value.toInt();
			}
			else
			{
				if (key == "AverageWeight")
				{
					avg_weight = value.toFloat();
				}
				else
				{
					if (key == "Isotopes")
					{
						distribution[UnsignedInt(split[3].toInt())] = double(value.toFloat()/100);
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
									mono_weight = value.toFloat();
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

	IsotopeDistribution ElementDB::parseIsotopeDistribution_(const HashMap<UnsignedInt, double>& distribution) 
		throw(Exception::ParseError)
	{
		IsotopeDistribution::ContainerType dist;
		for (HashMap<UnsignedInt, double>::ConstIterator it=distribution.begin(); it!=distribution.end(); ++it)
		{
			dist.push_back(make_pair<UnsignedInt, double>(it->first, it->second));
		}
	
		// TODO sorting!
		
		IsotopeDistribution iso_dist;
		iso_dist.set(dist);
		iso_dist.setMaxIsotope(100);
	
		return iso_dist;
	}
	
	void ElementDB::clear_()
	{
		HashMap<String, const Element*>::Iterator it=names_.begin();
		for (;it!=names_.end();++it)
		{
			delete it->second;
		}
		names_.clear();
		symbols_.clear();
		atomic_numbers_.clear();
	}
} // namespace OpenMS

