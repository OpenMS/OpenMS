// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>
#include <set>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModifierRep.h>
//////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModifierRep, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModifierRep* ptr = nullptr;
ModifierRep* nullPointer = nullptr;

START_SECTION(ModifierRep())
	ptr = new ModifierRep();
  TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION(~ModifierRep())
	delete ptr;
END_SECTION

START_SECTION(ModifierRep(const ModifierRep &source))
	ptr = new ModifierRep();
	ptr->setNumberOfModifications(2);
	ModifierRep * new_ptr = new ModifierRep(*ptr);
	TEST_EQUAL (ptr->getNumberOfModifications(),new_ptr->getNumberOfModifications());
	TEST_EQUAL (ptr->getModificationTable().size(),new_ptr->getModificationTable().size());
END_SECTION

START_SECTION(void setNumberOfModifications(Size i))
	ptr = new ModifierRep();
	TEST_EQUAL (0,ptr->getNumberOfModifications());
	ptr->setNumberOfModifications(1);
	TEST_EQUAL (1,ptr->getNumberOfModifications());
END_SECTION

START_SECTION(Size getNumberOfModifications() const )
	ptr = new ModifierRep();
	ptr->setNumberOfModifications(1);
	TEST_EQUAL (1,ptr->getNumberOfModifications());
	ptr->setNumberOfModifications(2);
	TEST_EQUAL (2,ptr->getNumberOfModifications());

END_SECTION



START_SECTION(const std::vector<std::vector<double> >& getModificationTable())
	/// missing functionality in ResidueDB? NOT_TESTABLE for now... (ek)
	#if 0
	ptr = new ModifierRep();
	std::vector<std::vector<double> > mod_table = ptr->getModificationTable();
	ResidueDB* rdb = ResidueDB::getInstance();
		
	char aa[] = "ARNDCEQGHILKMFPSTWYV";
	for (Size i = 0; i<strlen(aa);++i)
	{
		const Residue* r = rdb->getResidue(aa[i]);
		std::set< const ResidueModification * > mods (rdb->getModifications(r));
		std::set< const ResidueModification * >::iterator it (mods.begin());
		for(;it!=mods.end();++it)
		{
			double add_mass = (*it)->getDiffAverageMass();
			bool found_add = false;
			if (add_mass>0) 
			{
				for (Size j = 0 ; j <mod_table.at((int)aa[i]).size();++j)
				{
					if (mod_table.at((int)aa[i]).at(j)==add_mass) found_add = true;
				}
				TEST_EQUAL (found_add,1)
			}
			bool found_del = false;
			if (del_mass>0) 
			{
				for (Size j = 0 ; j <mod_table.at((int)aa[i]).size();++j)
				{
					if (mod_table.at((int)aa[i]).at(j)==-del_mass) found_del = true;
				}
				TEST_EQUAL (found_del,1)
			}
			std::cerr<<aa[i]<<" : "<<add_mass<<std::endl;
			std::cerr<<aa[i]<<" : "<<del_mass<<std::endl;
		}
	}
	#endif
	
	TEST_EQUAL(0,0)
END_SECTION

START_SECTION(Size getMaxModificationMasses())
	ptr = new ModifierRep();
	TEST_EQUAL(0,ptr->getMaxModificationMasses());
	std::vector<std::vector<double> > mod_table = ptr->getModificationTable();
	std::set<double> mod_masses_set;
	for (Size i = 0; i < mod_table.size();++i)
	{
		for (Size j = 0; j < mod_table.at(i).size();++j)
		{
			mod_masses_set.insert(mod_table.at(i).at(j));
		}
	}
	ptr->setNumberOfModifications(1);
	TEST_EQUAL(mod_masses_set.size(),ptr->getMaxModificationMasses());
END_SECTION

START_SECTION(void refreshModificationList(std::map< double, SignedSize > &mod_map, const char &c))
	
	#if 0
	ptr = new ModifierRep();
	std::map<double,int> mods;
	std::vector<std::vector<double> > mod_table = ptr->getModificationTable();
	const char aa = 'C';
	ptr->setNumberOfModifications(1);
	ptr->refreshModificationList(mods,aa);
	TEST_EQUAL (mods.size(),2);
	for (Size i = 0;i < mod_table.at((int)'C').size();++i)
	{
		TEST_EQUAL(mods[mod_table.at((int)'C').at(i)],1);
	}
	ptr->refreshModificationList(mods,aa);
	TEST_EQUAL (mods.size(),2);
	for (Size i = 0;i < mod_table.at((int)'C').size();++i)
	{
		TEST_EQUAL(mods[mod_table.at((int)'C').at(i)],1);
	}
	ptr->setNumberOfModifications(2);
	ptr->refreshModificationList(mods,aa);
	TEST_EQUAL (mods.size(),5);
	const char aa2 = 'X';
	ptr->refreshModificationList(mods,aa2);
	TEST_EQUAL (mods.size(),5);
	#endif
	
	TEST_EQUAL(0,0)
END_SECTION

START_SECTION(std::vector<String> getModificationsForMass(double &m))
	#if 0
	ptr = new ModifierRep();
	std::vector<std::vector<double> > mod_table = ptr->getModificationTable();
	ptr->setNumberOfModifications(2);
	char aa[] = "ARNDCEQGHILKMFPSTWYV";
	std::map<double,int> mods;
	for (Size i = 0; i < strlen(aa);++i){
		for (Size j = 0; j < strlen(aa);++j){
			if ( mod_table.at((int)aa[i]).size()>0 && mod_table.at((int)aa[j]).size()>0){
				double m = mod_table.at((int)aa[i]).at(0)+mod_table.at((int)aa[j]).at(0);
				vector <String> res = ptr->getModificationsForMass(m);
				bool found = false;
				String tupel = (String)aa[i]+(String)aa[j];
				String tupel_changed = (String)aa[j]+(String)aa[i];
				//std::cout<<std::endl<<tupel<< ":" << m<<std::endl;
				
				for (Size k = 0; k < res.size();++k)
				{
					//std::cout<<" f:"<<res.at(k)<<",";
					if (res.at(k)==tupel||res.at(k)==tupel_changed) found = true;
				}
				TEST_EQUAL (found, true);
			}
		}
	}
	#endif
	
	TEST_EQUAL(0,0)
END_SECTION


START_SECTION(std::vector<String> getModificationsForMass(double &m, const String &seq))
	#if 0
	ptr = new ModifierRep();
	std::vector<std::vector<double> > mod_table = ptr->getModificationTable();
	ptr->setNumberOfModifications(2);
	char aa[] = "ARNDCEQGHILKMFPSTWYV";
	std::map<double,int> mods;
	for (Size i = 0; i < strlen(aa);++i){
		for (Size j = 0; j < strlen(aa);++j){
			if ( mod_table.at((int)aa[i]).size()>0 && mod_table.at((int)aa[j]).size()>0){
				double m = mod_table.at((int)aa[i]).at(0)+mod_table.at((int)aa[j]).at(0);
				String tupel = (String)aa[i]+(String)aa[j];
				String tupel_changed = (String)aa[j]+(String)aa[i];
				vector <String> res = ptr->getModificationsForMass(m,tupel);
				TEST_EQUAL (res.at(0)==tupel||res.at(0)==tupel_changed,true)
				TEST_EQUAL (res.size(),1);
			}
		}
	}
	#endif
	
	TEST_EQUAL(0,0)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
