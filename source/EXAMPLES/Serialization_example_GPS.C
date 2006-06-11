// $Id: Serialization_example_GPS.C,v 1.1 2006/06/01 14:35:25 groepl Exp $
// $Author: groepl $
// $Maintainer: Clemens Groepl$


#include <fstream>
#include <string>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// include headers that implement a archive in simple xml format
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include <boost/serialization/nvp.hpp>

/////////////////////////////////////////////////////////////

struct BaseA
{

  BaseA() { name = "BaseA default";}
  BaseA(int rhs) { id[0] = rhs; id[1] = 2*rhs; name = "BaseA from int";}


 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int /* version */ )
  {
#if 1 // here we serialize the array by hand
    ar & boost::serialization::make_nvp("id0",id[0]);
    ar & boost::serialization::make_nvp("id1",id[1]);
#else // boost_serialization can iterate over the array using its own style
    ar & BOOST_SERIALIZATION_NVP(id);
#endif
    ar & BOOST_SERIALIZATION_NVP(name);
  }

 public:
  int id[2];
  std::string name;
};

/////////////////////////////////////////////////////////////
// gps coordinate
//
// illustrates serialization for a simple type
//
class gps_position : public BaseA
{
 private:
  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int /* version */ )
  {
    ar & boost::serialization::make_nvp("BaseA",boost::serialization::base_object<BaseA>(*this));
    ar & BOOST_SERIALIZATION_NVP(degrees);
    ar & BOOST_SERIALIZATION_NVP(minutes);
    ar & BOOST_SERIALIZATION_NVP(seconds);
  }
  int degrees;
  int minutes;
  float seconds;
 public:
  gps_position(){};
  gps_position(int d, int m, float s)
    : BaseA ( int(d+m+s) ),
      degrees(d), minutes(m), seconds(s)
  {}
};

int main() {
  // create and open a character archive for output
  std::ofstream text_ofs("archived.txt");

  // create and open a character archive for output
  std::ofstream xml_ofs("archived.xml");

  // create class instance
  const gps_position g(35, 59, 24.567f), h;

  // save data to archive
  {
    boost::archive::text_oarchive text_oa(text_ofs);
    // write class instance to archive
    text_oa << g << h;
    // archive and stream closed when destructors are called
  }
  {
    boost::archive::xml_oarchive xml_oa(xml_ofs);
    // write class instance to archive
    xml_oa << boost::serialization::make_nvp("position",g) << boost::serialization::make_nvp("another_position",h);
    // archive and stream closed when destructors are called
  }

  // ... some time later restore the class instance to its orginal state
  gps_position newg_from_text, newh_from_text;
  {
    // create and open an archive for input
    std::ifstream text_ifs("archived.txt", std::ios::binary);
    boost::archive::text_iarchive text_ia(text_ifs);
    // read class state from archive
    text_ia >> boost::serialization::make_nvp("position",newg_from_text) >> boost::serialization::make_nvp("another_position",newh_from_text);
    // archive and stream closed when destructors are called was the same
  }
  
  gps_position newg_from_xml, newh_from_xml;
  {
    // create and open an archive for input
    std::ifstream xml_ifs("archived.xml", std::ios::binary);
    boost::archive::xml_iarchive xml_ia(xml_ifs);
    // read class state from archive
    xml_ia >> boost::serialization::make_nvp("position",newg_from_xml) >> boost::serialization::make_nvp("another_position",newh_from_xml);
    // archive and stream closed when destructors are called was the same
  }

  // create and open a character archive for output
  std::ofstream text_ofs2("archived2.txt");

  // create and open a character archive for output
  std::ofstream xml_ofs2("archived2.xml");

  // save data to archive
  {
    boost::archive::text_oarchive text_oa2(text_ofs2);
    // write class instance to archive
    text_oa2 << boost::serialization::make_nvp("position",newg_from_text) << boost::serialization::make_nvp("another_position",newh_from_text );
    // archive and stream closed when destructors are called
  }
  {
    boost::archive::xml_oarchive xml_oa2(xml_ofs2);
    // write class instance to archive
    xml_oa2 << boost::serialization::make_nvp("position",newg_from_xml) << boost::serialization::make_nvp("another_position",newh_from_xml );
    // archive and stream closed when destructors are called
  }
  
  

  return 0;
}

