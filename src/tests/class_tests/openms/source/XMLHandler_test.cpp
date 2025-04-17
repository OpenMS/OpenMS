
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
    0x0072,0x006C,0x0064,0x0021};
XMLSize_t a_length = xercesc::XMLString::stringLen(ascii); 

const XMLCh mixed[] = { 
    0x0048, 0x0065,0x0432, 0x0435, 0x0442, 0x043C, 0x006F,
    0x0072,0x006C,0x0064, 0x0021 };
XMLSize_t m_length = xercesc::XMLString::stringLen(mixed);

const XMLCh empty[] = {0};
XMLSize_t e_length = xercesc::XMLString::stringLen(empty);
std::cout << e_length << std::endl;

const XMLCh upperBoundary [] = {0x00FF,0x00FF};
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

const XMLCh eight_block_negative[] = {0xFFFF,0xFFFE,0xFFFB,0xFFF6,0xFFEC,0xFFCE,0xFF9C,0xFE00};

const XMLCh eight_block[] = {0x0048,0x0065,0x006C,0x006C,0x006F,0x002C,0x0057,0x006F};

const XMLCh eight_block_mixed[] ={0x0042,0x0045,0x004C,0x0041,0xFFFF,0xFFFE,0xFFFB,0xFFF6};

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
    char* output1 = new char [9];
    output1[8] = '\0';
    StringManager::compress64(eight_block,output1);
    std::string res1_str = "Hello,Wo";
    std::string o1_str (output1);
    TEST_STRING_EQUAL(o1_str,res1_str);
    delete[] output1;
    
   
    char* output2 = new char [9];
    output2 [8] = '\0';
    char res2 [9] = {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};
    res2[8] = '\0';
    StringManager::compress64(eight_block_negative,output2);
    std::string res2_str (res2);
    std::string o2_str(output2);
    TEST_STRING_EQUAL(o2_str, res2_str);
    delete[] output2;

    char* output3 = new char [9];
    output3 [8] = '\0';
    char res3 [9] = {0x42,0x45,0x4C,0x41,0x00,0x00,0x00,0x00};
    res3[8] = '\0';
    StringManager::compress64(eight_block_mixed,output3);
    std::string res3_str (res3);
    std::string o3_str(output3);
    TEST_STRING_EQUAL(o3_str, res3_str);
    delete[] output3;

    char* output4 = new char [13];
    output4 [0]  ='A';
    output4 [1]  ='B';
    output4 [2]  ='R';
    output4 [3]  ='A';
    output3 [12] = '\0';
    
    StringManager::compress64(eight_block_kadabra,(output4+4));
    std::string res4_str = "ABRAKADABRA!";
    std::string o4_str(output4);
    TEST_STRING_EQUAL(o4_str, res4_str);
    delete[] output4;

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

END_TEST


    

