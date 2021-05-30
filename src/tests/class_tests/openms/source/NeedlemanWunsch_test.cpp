#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/SEQUENCE/NeedlemanWunsch.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(NeedlemanWunsch, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

NeedlemanWunsch* ptr = nullptr;
START_SECTION(NeedlemanWunsch(ScoringMatrix matrix, int penalty)())
{
  ptr = new NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::PAM30MS, 5);
  TEST_EQUAL(ptr == nullptr, false)
}
END_SECTION

START_SECTION(~NeedlemanWunsch())
{
  delete (ptr);
}
END_SECTION

String seq1 = "IGGATLIGQLAIQQAHVHL";
String seq2 = "IGGATLIGALDQVVAQQAHVHL";

START_SECTION(double align(const String& seq1, const String& seq2))
{
  NeedlemanWunsch alignment = NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::identity, 5);
  TEST_EQUAL(alignment.align(seq1, seq2), 1);
  TEST_EQUAL(alignment.align(seq1, seq1), 19);
  TEST_EQUAL(alignment.align(seq2, seq2), 22);
}
END_SECTION

START_SECTION(void setMatrix(const ScoringMatrix& matrix))
{
  NeedlemanWunsch alignment = NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::identity, 5);
  alignment.setMatrix(NeedlemanWunsch::ScoringMatrix::PAM30MS);
  TEST_EQUAL(alignment.align(seq1, seq2), 93);
  TEST_EQUAL(alignment.align(seq1, seq1), 131);
  TEST_EQUAL(alignment.align(seq2, seq2), 151);
}
END_SECTION

START_SECTION(void setMatrix(const std::string& matrix))
{
  NeedlemanWunsch alignment = NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::PAM30MS, 5);

  TEST_EXCEPTION(Exception::IllegalArgument, alignment.setMatrix("Identity"))

  alignment.setMatrix("identity");
  TEST_EQUAL(alignment.align(seq1, seq2), 1);
  TEST_EQUAL(alignment.align(seq1, seq1), 19);
  TEST_EQUAL(alignment.align(seq2, seq2), 22);
}
END_SECTION

START_SECTION(void setPenalty(const int penalty))
{
  NeedlemanWunsch alignment = NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::PAM30MS, 5);
  alignment.setPenalty(1);
  TEST_EQUAL(alignment.align(seq1, seq2), 113);
  TEST_EQUAL(alignment.getPenalty(), 1);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

