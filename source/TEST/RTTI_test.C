// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/RTTI.h>

///////////////////////////

START_TEST(RTTI, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::RTTI;

CHECK(RTTI::isKindOf<>())
//TODO 
//	Protein p;
//	Protein*			p_ptr(&p);
//	Molecule*			m_ptr(&p);
//	AtomContainer*	b_ptr(&p);
//	TEST_EQUAL(isKindOf<Molecule>(*p_ptr), true)
//	TEST_EQUAL(isKindOf<Molecule>(*m_ptr), true)
//	TEST_EQUAL(isKindOf<Molecule>(*b_ptr), true)
//	TEST_EQUAL(isKindOf<Protein>(*p_ptr), true)
//	TEST_EQUAL(isKindOf<Protein>(*m_ptr), true)
//	TEST_EQUAL(isKindOf<Protein>(*b_ptr), true)
//	TEST_EQUAL(isKindOf<AtomContainer>(*p_ptr), true)
//	TEST_EQUAL(isKindOf<AtomContainer>(*m_ptr), true)
//	TEST_EQUAL(isKindOf<AtomContainer>(*b_ptr), true)
//	TEST_EQUAL(isKindOf<Residue>(*p_ptr), false)
//	TEST_EQUAL(isKindOf<Residue>(*m_ptr), false)
//	TEST_EQUAL(isKindOf<Residue>(*b_ptr), false)
//	Molecule m;
//	m_ptr = &m;
//	b_ptr = &m;
//	TEST_EQUAL(isKindOf<Molecule>(*m_ptr), true)
//	TEST_EQUAL(isKindOf<Molecule>(*b_ptr), true)
//	TEST_EQUAL(isKindOf<Protein>(*m_ptr), false)
//	TEST_EQUAL(isKindOf<Protein>(*b_ptr), false)
//	TEST_EQUAL(isKindOf<AtomContainer>(*m_ptr), true)
//	TEST_EQUAL(isKindOf<AtomContainer>(*b_ptr), true)
//	TEST_EQUAL(isKindOf<Residue>(*m_ptr), false)
//	TEST_EQUAL(isKindOf<Residue>(*b_ptr), false)
//	AtomContainer b;
//	b_ptr = &b;
//	TEST_EQUAL(isKindOf<Molecule>(*b_ptr), false)
//	TEST_EQUAL(isKindOf<Protein>(*b_ptr), false)
//	TEST_EQUAL(isKindOf<AtomContainer>(*b_ptr), true)
//	TEST_EQUAL(isKindOf<Residue>(*b_ptr), false)
//
//	System*	s_ptr = new System;
//	TEST_EQUAL(isKindOf<System>(*s_ptr), true)
//	delete s_ptr;
//
//	s_ptr = 0;
//	TEST_EQUAL(isKindOf<System>(*s_ptr), false)
RESULT											

CHECK(RTTI::isInstanceOf<>())
//TODO 
//	Protein p;
//	Protein*			p_ptr(&p);
//	Molecule*			m_ptr(&p);
//	AtomContainer*	b_ptr(&p);
//	TEST_EQUAL(isInstanceOf<Molecule>(*p_ptr), false)
//	TEST_EQUAL(isInstanceOf<Molecule>(*m_ptr), false)
//	TEST_EQUAL(isInstanceOf<Molecule>(*b_ptr), false)
//	TEST_EQUAL(isInstanceOf<Protein>(*p_ptr), true)
//	TEST_EQUAL(isInstanceOf<Protein>(*m_ptr), true)
//	TEST_EQUAL(isInstanceOf<Protein>(*b_ptr), true)
//	TEST_EQUAL(isInstanceOf<AtomContainer>(*p_ptr), false)
//	TEST_EQUAL(isInstanceOf<AtomContainer>(*m_ptr), false)
//	TEST_EQUAL(isInstanceOf<AtomContainer>(*b_ptr), false)
//	TEST_EQUAL(isInstanceOf<Residue>(*p_ptr), false)
//	TEST_EQUAL(isInstanceOf<Residue>(*m_ptr), false)
//	TEST_EQUAL(isInstanceOf<Residue>(*b_ptr), false)
//	Molecule m;
//	m_ptr = &m;
//	b_ptr = &m;
//	TEST_EQUAL(isInstanceOf<Molecule>(*m_ptr), true)
//	TEST_EQUAL(isInstanceOf<Molecule>(*b_ptr), true)
//	TEST_EQUAL(isInstanceOf<Protein>(*m_ptr), false)
//	TEST_EQUAL(isInstanceOf<Protein>(*b_ptr), false)
//	TEST_EQUAL(isInstanceOf<AtomContainer>(*m_ptr), false)
//	TEST_EQUAL(isInstanceOf<AtomContainer>(*b_ptr), false)
//	TEST_EQUAL(isInstanceOf<Residue>(*m_ptr), false)
//	TEST_EQUAL(isInstanceOf<Residue>(*b_ptr), false)
//	AtomContainer b;
//	b_ptr = &b;
//	TEST_EQUAL(isInstanceOf<Molecule>(*b_ptr), false)
//	TEST_EQUAL(isInstanceOf<Protein>(*b_ptr), false)
//	TEST_EQUAL(isInstanceOf<AtomContainer>(*b_ptr), true)
//	TEST_EQUAL(isInstanceOf<Residue>(*b_ptr), false)
//
//	Atom*					a_ptr = new Atom;
//	System*				s_ptr = new System;
//	PDBAtom*		 pa_ptr = new PDBAtom;
//	Bond*				 bo_ptr = new Bond;
//	NucleicAcid* na_ptr = new NucleicAcid;
//	Nucleotide*  nu_ptr = new Nucleotide;
//	Chain*				c_ptr = new Chain;
//	SecondaryStructure* ss_ptr = new SecondaryStructure;
//	TEST_EQUAL(isInstanceOf<Atom>(*a_ptr), true)
//	TEST_EQUAL(isInstanceOf<System>(*s_ptr), true)
//	TEST_EQUAL(isInstanceOf<PDBAtom>(*pa_ptr), true)
//	TEST_EQUAL(isInstanceOf<Bond>(*bo_ptr), true)
//	TEST_EQUAL(isInstanceOf<NucleicAcid>(*na_ptr), true)
//	TEST_EQUAL(isInstanceOf<Nucleotide>(*nu_ptr), true)
//	TEST_EQUAL(isInstanceOf<Chain>(*c_ptr), true)
//	TEST_EQUAL(isInstanceOf<SecondaryStructure>(*ss_ptr), true)
//	delete a_ptr;
//	delete s_ptr;
//	delete pa_ptr;
//	delete bo_ptr;
//	delete na_ptr;
//	delete nu_ptr;
//	delete c_ptr;
//	delete ss_ptr;
RESULT											

CHECK(getDefault<>())
//TODO 
//	const Protein& p = getDefault<Protein>();
//	TEST_EQUAL(p.isValid(), true)
//	const Atom& a = getDefault<Atom>();
//	TEST_EQUAL(a.isValid(), true)
//	const System& s = getDefault<System>();
//	TEST_EQUAL(s.isValid(), true)
//	PDBAtom pa;
//	pa = getDefault<PDBAtom>();
//	Bond b;
//	b = getDefault<Bond>();
//	NucleicAcid na;
//	na = getDefault<NucleicAcid>();
//	Nucleotide n;
//	n = getDefault<Nucleotide>();
//	SecondaryStructure ss;
//	ss = getDefault<SecondaryStructure>();
//	Residue r;
//	r = getDefault<Residue>();
//	Chain c;
//	c = getDefault<Chain>();
//	Molecule m;
//	m = getDefault<Molecule>();
//	AtomContainer bf;
//	bf = getDefault<AtomContainer>();
//	Fragment f;
//	f = getDefault<Fragment>();
RESULT

CHECK(getNew<>())
//TODO 
//	Protein* p = (Protein*)getNew<Protein>();
//	TEST_NOT_EQUAL(p, 0)
//	delete p;
RESULT

CHECK(getName<>())
// there is not much to check - each damned compiler 
// tries his own demangling!
//TODO 
//TEST_EQUAL(String(getName<Protein>()).hasSubstring("Protein"), true)
RESULT

CHECK(getStreamName<>())
	// there is not much to check - each damned compiler 
	// tries his own demangling!
	TEST_EQUAL(getStreamName<float>(), String("float"))
	TEST_EQUAL(getStreamName<Index>(), String("OpenMS::Index"))
	TEST_EQUAL(getStreamName<bool>(), String("bool"))
	TEST_EQUAL(getStreamName<char>(), String("char"))
	TEST_EQUAL(getStreamName<double>(), String("double"))
RESULT

CHECK(castTo<>())
	//TODO
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
