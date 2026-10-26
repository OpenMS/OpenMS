
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <string>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CONCEPT/Init.h>
#include <random>

class StringManager_test : public OpenMS::Internal::StringManager
{
public:
  StringManager_test() = default;
  ~StringManager_test() = default;

  static void compress64(const XMLCh* input_it, char* output_it)
  {
    StringManager::compress64_(input_it, output_it);
  }
};

using namespace OpenMS::Internal;





START_TEST(StringManager, "$Id$")


const XMLCh russianHello[] = {
    0x041F, 0x0440, 0x0438, 0x0432, 0x0435, 0x0442, 0x043C,
    0x0438, 0x0440,0x0000 // "Привет мир" (Hello World in Russian)
};
size_t r_length = StringManager::strLength(russianHello);

const XMLCh ascii[] = {
    0x0048,0x0065,0x006C,0x006C,0x006F,0x002C,0x0057,0x006F,
    0x0072,0x006C,0x0064,0x0021, 0x0000};
size_t a_length = StringManager::strLength(ascii);

const XMLCh mixed[] = {
    0x0048, 0x0065,0x0432, 0x0435, 0x0442, 0x043C, 0x006F,
    0x0072,0x006C,0x0064, 0x0021, 0x0000 };
size_t m_length = StringManager::strLength(mixed);

const XMLCh empty[] = {0};
size_t e_length = StringManager::strLength(empty);

const XMLCh upperBoundary [] = {0x00FF,0x00FF,0x0000};
size_t u_length = StringManager::strLength(upperBoundary);

bool isAscii = false;

START_SECTION(isASCII(const XMLCh * chars, const size_t length))
  isAscii = StringManager::isASCII(ascii,a_length);
  TEST_TRUE(isAscii)

  isAscii = StringManager::isASCII(russianHello,r_length);
  TEST_FALSE(isAscii)

  isAscii  = StringManager::isASCII(mixed,m_length);
  TEST_FALSE(isAscii)

  isAscii = StringManager::isASCII(empty,e_length);
  TEST_TRUE(isAscii)

  // U+00FF is Latin-1, not ASCII: it needs the UTF-8 transcoding path
  isAscii = StringManager::isASCII(upperBoundary,u_length);
  TEST_FALSE(isAscii)
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
    StringManager_test::compress64(eight_block,o1_str.data());
    std::string res1_str = "Hello,Wo";
    TEST_STRING_EQUAL(o1_str,res1_str);


    std::string o2_str(8,'\0');
    StringManager_test::compress64(eight_block_negative,o2_str.data());
    std::string res2_str = res1_str;
    TEST_STRING_EQUAL(o2_str, res2_str);


    std::string o3_str(8,'\0');
    StringManager_test::compress64(eight_block_mixed,o3_str.data());
    std::string res3_str = {0x42,0x45,0x4C,0x41,0x42,0x45,0x4C,0x41};
    TEST_STRING_EQUAL(o3_str, res3_str);

    std::string o4_str(12,'\0');
    o4_str [0]  ='A';
    o4_str [1]  ='B';
    o4_str [2]  ='R';
    o4_str [3]  ='A';

    StringManager_test::compress64(eight_block_kadabra,((o4_str.data())+4));
    std::string res4_str = "ABRAKADABRA!";
    TEST_STRING_EQUAL(o4_str, res4_str);

END_SECTION

//Tests Number of Chars not Dividable by 8
std::string o5_str;
std::string res5_str = "Hello,World!";

//Checks how the Function handles Data thats already stored in Output string
std::string o6_str = "Gruess Gott und ";
std::string res6_str = "Gruess Gott und Hello,World!";

std::string o7_str;
std::string res7_str;


START_SECTION(appendASCII(const XMLCh * chars, const size_t length, String & result))

    StringManager::appendASCII(ascii,a_length,o5_str);
    TEST_STRING_EQUAL(o5_str, res5_str);

    StringManager::appendASCII(ascii,a_length,o6_str);
    TEST_STRING_EQUAL(o6_str, res6_str);

    StringManager::appendASCII(empty,e_length,o7_str);
    TEST_STRING_EQUAL(o7_str, res7_str);


END_SECTION
START_SECTION(static bool widensASCII())
  OpenMS::Internal::xercesInitialize();
  // the tests run with an ASCII-compatible locale (C or UTF-8)
  TEST_TRUE(StringManager::widensASCII())
END_SECTION

START_SECTION(static XercesString convert(const char * str))
  // ASCII text is widened directly: the same as the transcoder's result (convertPtr())
  std::string all_ascii;
  for (int c = 1; c < 128; ++c) { all_ascii += static_cast<char>(c); }
  for (const std::string& text : {std::string(""), std::string("mono_mass"), all_ascii, std::string("caf\xC3\xA9")})
  {
    const std::basic_string<char16_t> widened = StringManager::convert(text);
    const unique_xerces_ptr<char16_t> transcoded = StringManager::convertPtr(text);
    TEST_TRUE(widened == std::basic_string<char16_t>(transcoded.get()))
  }
END_SECTION

START_SECTION(static Int parseInt(const char16_t * str))
  // as XMLString::parseInt() of Xerces 3.3 (checked against it): XML white space trimmed, then all of it a base-10 long
  // (strtol, which also skips leading \\v and \\f) cut to an int, else it throws
  // (the strings that are not ASCII or longer than 31 characters go to XMLString::parseInt() itself)
  const auto outcome = [](const std::u16string& text) -> std::pair<bool, long long>
  {
    try { return {true, StringManager::parseInt(text.c_str())}; }
    catch (...) { return {false, 0}; }
  };
  const std::vector<std::pair<std::u16string, long long>> taken = {{u"0", 0}, {u"-0", 0}, {u"+5", 5}, {u" 12", 12}, {u"12 ", 12},
    {u"\t-42\n", -42}, {u"\r\n7\r\n", 7}, {u"2147483647", 2147483647LL}, {u"-2147483648", -2147483648LL}, {u" +3 ", 3},
    {u"007", 7}, {u"000000000000000000000000000000000001", 1}, {u"-000000000000000000000000000000001", -1},
    {u"2147483648", -2147483648LL}, {u"-2147483649", 2147483647LL}, {u"4294967297", 1}, {u"9223372036854775807", -1}, {u"\v12", 12},
    {u"\f12", 12}};
  for (const auto& [text, value] : taken)
  {
    const auto o = outcome(text);
    TEST_TRUE(o.first)
    TEST_EQUAL(o.second, value)
  }
  const std::vector<std::u16string> refused = {u"99999999999999999999", u"9223372036854775808", u"12\v",
    u"", u"   ", u"1a", u"0x10", u"1 2", u"+", u"-", u"1e3", u"12.0", u"\u0661\u0662", u"\uFF11\uFF12", u"12\u00A0", u"\u00A012"};
  for (const std::u16string& text : refused) { TEST_FALSE(outcome(text).first) }
  std::mt19937_64 rng(42);
  for (int i = 0; i < 1000; ++i)
  {
    const long long v = static_cast<long long>(rng() >> (rng() % 64)) * ((rng() & 1) ? -1 : 1);
    const std::string t = std::to_string(v);
    const auto o = outcome(std::u16string(t.begin(), t.end()));
    TEST_TRUE(o.first)
    TEST_EQUAL(o.second, static_cast<long long>(static_cast<OpenMS::Int>(v)))
  }
END_SECTION

XMLCh* nullPointer = nullptr;
START_SECTION(strLength(const XMLCh* input_ptr))
    int o_length = StringManager::strLength(ascii);
    TEST_EQUAL(o_length, a_length);
    o_length = StringManager::strLength(empty);
    TEST_EQUAL(o_length, e_length);
    o_length = StringManager::strLength(upperBoundary);
    TEST_EQUAL(o_length, u_length);
    o_length = StringManager::strLength(nullPointer);
    TEST_EQUAL(o_length, 0);
END_SECTION

END_TEST




