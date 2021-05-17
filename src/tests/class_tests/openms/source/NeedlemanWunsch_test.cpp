#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <iostream>
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
  ptr = new NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::PAM30MSMatrix, -5);
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

START_SECTION(double align_(const String& seq1, const String& seq2))
{
  NeedlemanWunsch object = NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::identityMatrix, -5);
  TEST_EQUAL(object.align_(seq1, seq2), 1);
  TEST_EQUAL(object.align_(seq1, seq1), 19);
  TEST_EQUAL(object.align_(seq2, seq2), 22);
}
END_SECTION

START_SECTION(void setMatrix_(const ScoringMatrix& matrix))
{
  NeedlemanWunsch object = NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::identityMatrix, -5);
  object.setMatrix_(NeedlemanWunsch::ScoringMatrix::PAM30MSMatrix);
  TEST_EQUAL(object.align_(seq1, seq2), 93);
  TEST_EQUAL(object.align_(seq1, seq1), 131);
  TEST_EQUAL(object.align_(seq2, seq2), 151);
  //TEST_EQUAL(object.getMatrix_(), NeedlemanWunsch::ScoringMatrix::PAM30MSMatrix); kein == operator definiert für ScoringMatrix?
}
END_SECTION

START_SECTION(void setPenalty_(const ScoringMatrix& matrix))
{
  NeedlemanWunsch object = NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix::PAM30MSMatrix, -5);
  TEST_EXCEPTION(Exception::IllegalArgument, object.setPenalty_(5))
  object.setPenalty_(-1);
  TEST_EQUAL(object.align_(seq1, seq2), 113);
  TEST_EQUAL(object.getPenalty_(), -1);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

