// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_CHEMISTRY_MODIFIERREP_H
#define OPENMS_CHEMISTRY_MODIFIERREP_H

#include <map>
#include <vector>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS {

/**
@brief Implements modification for suffix arrays
*/

class OPENMS_DLLAPI ModifierRep {

public:
/**
@brief Constructor
*/
ModifierRep ();
/**
@brief copy constructor
*/
ModifierRep (const ModifierRep & source);
/**
@brief destructor
*/
virtual ~ModifierRep ();

/**
@brief setter for number of modifications
@param i number of modifications
*/
void setNumberOfModifications(Size i);

/**
@brief getter for number of modifications
@return number of modifications
*/
Size getNumberOfModifications() const;

/**
@brief getter for the modification table
@return const reference to the modofication table
*/
const std::vector<std::vector<double> > & getModificationTable ();
/**
@brief updates the modifications list if with the modifications of a new amino acid
@param mod_map reference to the map holding the possible modifications
@param c const character for the amino acid

all modification for the given amino acid are added to the map. The key is the mass, the value is the number of modifications for that mass. So the advantage of using a map is that a value will occure only once. the mass of the modification is added to all elements in the map whose number of modifications is smaller than the maximal number of modifications.
*/
void refreshModificationList (std::map<double, SignedSize> & mod_map,const char & c);

/**
@brief calculates the maximal number unique masses of combinations of modifications (maximal possible size of the modification map)
@return Size
*/
Size getMaxModificationMasses ();

/**
@brief gets all modification possibilities for a given mass
@param m masse
@return vector of strings
therefor at first the massmapping is calculated. This massmapping will be saved as private member mass_mapping_ so that it has to be calculated only once.
*/
std::vector<String> getModificationsForMass (double & m);

/**
@brief gets all modification possibilities for a given mass and for the given peptide
@param m masse
@param seq peptide sequence
@return vector of strings
the getModificationsForMass (double & m) will be used. then a chacater histogramm of the sequence is created as well as for every possible combination of modifications. Then only modifications that are possible are returned.
*/
std::vector<String> getModificationsForMass (double & m, const String & seq);


protected:

std::vector<std::vector<double> > modification_table_; ///< all possible modifications

Size number_of_modifications_; ///< number of maximal modifications

std::map<String,std::vector<String> > mass_mapping_; ///< maps a mass to the combination of modifications

};

}

#endif //OPENMS_CHEMISTRY_MODIFIERREP_H
