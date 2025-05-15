
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <iostream>
#include <string>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

using namespace OpenMS::Internal;





START_TEST(StringManager, "$Id$")


const XMLCh russianHello[] = {
    0x041F, 0x0440, 0x0438, 0x0432, 0x0435, 0x0442, 0x043C, 
    0x0438, 0x0440, // "Привет мир" (Hello World in Russian)
};
XMLSize_t r_length = xercesc::XMLString::stringLen(russianHello);

const XMLCh ascii[] = { 
    0x0048,0x0065,0x006C,0x006C,0x006F,0x002C,0x0057,0x006F,
    0x0072,0x006C,0x0064,0x0021, 0x0000};
XMLSize_t a_length = xercesc::XMLString::stringLen(ascii); 

const XMLCh mixed[] = { 
    0x0048, 0x0065,0x0432, 0x0435, 0x0442, 0x043C, 0x006F,
    0x0072,0x006C,0x0064, 0x0021, 0x0000 };
XMLSize_t m_length = xercesc::XMLString::stringLen(mixed);

const XMLCh empty[] = {0};
XMLSize_t e_length = xercesc::XMLString::stringLen(empty);

const XMLCh upperBoundary [] = {0x00FF,0x00FF,0x0000};
XMLSize_t u_length = xercesc::XMLString::stringLen(upperBoundary);

bool isAscii = false;

START_SECTION(isASCII(const XMLCh * chars, const XMLSize_t length))
  isAscii = StringManager::isASCII(ascii,a_length);
  std::cout << "1 \n";
  TEST_TRUE(isAscii)
  isAscii = StringManager::isASCII(russianHello,r_length);
  std::cout << "2 \n";
  TEST_FALSE(isAscii)
  isAscii  = StringManager::isASCII(mixed,m_length);
  std::cout << "3 \n";
  TEST_FALSE(isAscii)
  isAscii = StringManager::isASCII(empty,e_length);
  std::cout << "4 \n";
  TEST_FALSE(isAscii)
  isAscii = StringManager::isASCII(upperBoundary,u_length);
  std::cout << "5 \n";
  TEST_TRUE(isAscii)
END_SECTION

const XMLCh eight_block_negative[] = {0x0148,0x0165,0x016C,0x016C,0x016F,0x012C,0x0157,0x016F};

const XMLCh eight_block[] = {0x0048,0x0065,0x006C,0x006C,0x006F,0x002C,0x0057,0x006F};

const XMLCh eight_block_mixed[] ={0x0042,0x0045,0x004C,0x0041,0x0142,0x0145,0x014C,0x0141};

const XMLCh eight_block_kadabra[] = {
    0x004B, // K
    0x0041, // A
    0x0044, // D
    0x0041, // A
    0x0042, // B
    0x0052, // R
    0x0041, // A
    0x0021  // !
};

START_SECTION(compress64 (const XMLCh* input_it, char* output_it))
    std::string o1_str(8,'\0');
    StringManager::compress64(eight_block,o1_str.data());
    std::string res1_str = "Hello,Wo";
    TEST_STRING_EQUAL(o1_str,res1_str);
    
   
    std::string o2_str(8,'\0'); 
    StringManager::compress64(eight_block_negative,o2_str.data());
    std::string res2_str = res1_str;
    TEST_STRING_EQUAL(o2_str, res2_str);

    
    std::string o3_str(8,'\0');
    // char res3 [9] = {0x42,0x45,0x4C,0x41,0x42,0x45,0x4C,0x41};
    // res3[8] = '\0';
    StringManager::compress64(eight_block_mixed,o3_str.data());
    std::string res3_str = {0x42,0x45,0x4C,0x41,0x42,0x45,0x4C,0x41};
    TEST_STRING_EQUAL(o3_str, res3_str);

    std::string o4_str(12,'\0');
    o4_str [0]  ='A';
    o4_str [1]  ='B';
    o4_str [2]  ='R';
    o4_str [3]  ='A';
    
    StringManager::compress64(eight_block_kadabra,((o4_str.data())+4));
    std::string res4_str = "ABRAKADABRA!";
    TEST_STRING_EQUAL(o4_str, res4_str);

END_SECTION

//Tests Number of Chars not Dividable by 8
OpenMS::String o5_str;
std::string res5_str = "Hello,World!";

//Checks how the Function handles Data thats already stored in Output string
OpenMS::String o6_str = "Gruess Gott und ";
std::string res6_str = "Gruess Gott und Hello,World!";

OpenMS::String o7_str;
std::string res7_str = "";


START_SECTION(appendASCII(const XMLCh * chars, const XMLSize_t length, String & result))

    StringManager::appendASCII(ascii,a_length,o5_str);
    TEST_STRING_EQUAL(o5_str, res5_str);

    StringManager::appendASCII(ascii,a_length,o6_str);
    TEST_STRING_EQUAL(o6_str, res6_str);

    StringManager::appendASCII(empty,e_length,o7_str);
    TEST_STRING_EQUAL(o7_str, res7_str);
    std::cout << o7_str.size() << std::endl;

END_SECTION
XMLCh* nullPointer = nullptr;
START_SECTION(appendASCII(const XMLCh * chars, const XMLSize_t length, String & result))
    int o_length = StringManager::strLength(ascii);
    TEST_EQUAL(o_length, a_length);
    o_length = StringManager::strLength(empty);
    TEST_EQUAL(o_length, e_length);
    o_length = StringManager::strLength(upperBoundary);
    TEST_EQUAL(o_length, u_length);
    o_length = StringManager::strLength(nullPointer);
END_SECTION

END_TEST


    

