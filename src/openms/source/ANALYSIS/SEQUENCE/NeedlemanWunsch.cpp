#include <OpenMS/ANALYSIS/SEQUENCE/NeedlemanWunsch.h>
#include <iostream>
#include <OpenMS/CONCEPT/Exception.h>
#include  <utility> //swap

using namespace std;
namespace OpenMS
{

std::vector<int> adaptedIdentity
    {
        //      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
        /* A */ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* R */ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* N */ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* D */ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* C */ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* Q */ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* E */ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* G */ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* H */ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* I */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* L */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* K */ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* M */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* F */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* P */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* S */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -17,
        /* T */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -17,
        /* W */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -17,
        /* Y */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -17,
        /* V */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -17,
        /* B */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -17,
        /* Z */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -17,
        /* X */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* * */ -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, 1
    };

std::vector<int> PAM30MS
    {
        //        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
        /* A */   6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5, -8, -2,  0, -1,-13, -8, -2, -7, -6,  0,-17,
        /* R */  -7,  8, -6,-10, -8, -2, -9, -9, -2, -5, -7,  0, -4, -9, -4, -3, -6, -2,-10, -8,  5, -1,  0,-17,
        /* N */  -4, -6,  8,  2,-11, -3, -2, -3,  0, -5, -6, -1, -9, -9, -6,  0, -2, -8, -4, -8, -4, -2,  0,-17,
        /* D */  -3,-10,  2,  8,-14, -2,  2, -3, -4, -7,-10, -4,-11,-15, -8, -4, -5,-15,-11, -8, -7, -3,  0,-17,
        /* C */  -6, -8,-11,-14, 10,-14,-14, -9, -7, -6,-11,-14,-13,-13, -8, -3, -8,-15, -4, -6,-11,-14,  0,-17,
        /* Q */  -4, -2, -3, -2,-14,  8,  1, -7,  1, -8, -7, -3, -4,-13, -3, -5, -5,-13,-12, -7, -3,  4,  0,-17,
        /* E */  -2, -9, -2,  2,-14,  1,  8, -4, -5, -5, -7, -4, -7,-14, -5, -4, -6,-17, -8, -6, -7, -2,  0,-17,
        /* G */  -2, -9, -3, -3, -9, -7, -4,  6, -9,-11,-11, -7, -8, -9, -6, -2, -6,-15,-14, -5, -8, -7,  0,-17,
        /* H */  -7, -2,  0, -4, -7,  1, -5, -9,  9, -9, -8, -6,-10, -6, -4, -6, -7, -7, -3, -6, -4, -3,  0,-17,
        /* I */  -5, -5, -5, -7, -6, -8, -5,-11, -9,  8,  5, -6, -1, -2, -8, -7, -2,-14, -6,  2, -6, -7,  0,-17,
        /* L */  -6, -7, -6,-10,-11, -7, -7,-11, -8,  5,  5, -7,  0, -3, -8, -8, -5,-10, -7,  0, -7, -7,  0,-17,
        /* K */  -7,  0, -1, -4,-14, -3, -4, -7, -6, -6, -7,  7, -2,-14, -6, -4, -3,-12, -9, -9,  5,  4,  0,-17,
        /* M */  -5, -4, -9,-11,-13, -4, -7, -8,-10, -1,  0, -2, 11, -4, -8, -5, -4,-13,-11, -1, -3, -3,  0,-17,
        /* F */  -8, -9, -9,-15,-13,-13,-14, -9, -6, -2, -3,-14, -4,  9,-10, -6, -9, -4,  2, -8,-12,-14,  0,-17,
        /* P */  -2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -8, -6, -8,-10,  8, -2, -4,-14,-13, -6, -5, -5,  0,-17,
        /* S */   0, -3,  0, -4, -3, -5, -4, -2, -6, -7, -8, -4, -5, -6, -2,  6,  0, -5, -7, -6, -4, -5,  0,-17,
        /* T */  -1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -5, -3, -4, -9, -4,  0,  7,-13, -6, -3, -5, -4,  0,-17,
        /* W */ -13, -2, -8,-15,-15,-13,-17,-15, -7,-14,-10,-12,-13, -4,-14, -5,-13, 13, -5,-15, -7,-13,  0,-17,
        /* Y */  -8,-10, -4,-11, -4,-12, -8,-14, -3, -6, -7, -9,-11,  2,-13, -7, -6, -5, 10, -7,-10,-11,  0,-17,
        /* V */  -2, -8, -8, -8, -6, -7, -6, -5, -6,  2,  0, -9, -1, -8, -6, -6, -3,-15, -7,  7, -9, -8,  0,-17,
        /* B */  -7,  5, -4, -7,-11, -3, -7, -8, -4, -6, -7,  5, -3,-12, -5, -4, -5, -7,-10, -9,  5,  1,  0,-17,
        /* Z */  -6, -1, -2, -3,-14,  4, -2, -7, -3, -7, -7,  4, -3,-14, -5, -5, -4,-13,-11, -8,  1,  4,  0,-17,
        /* X */   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,-17,
        /* * */ -17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,  1
    };

  NeedlemanWunsch::NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix matrix, int penalty)
{
 if (penalty >= 0)
 {
   String msg = "Gap penalty should be negative";
   throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    msg);
 }

  gapPenalty_ = penalty;

  if (matrix == ScoringMatrix::identityMatrix)
  {
    matrixPtr_ = &adaptedIdentity;
  }

  else if (matrix == ScoringMatrix::PAM30MSMatrix)
  {
    matrixPtr_ = &PAM30MS;
  }
  else
  {
    String msg = "Matrix is not known! Valid choices are: "
                                       "'identityMatrix', 'PAM30MSMatrix'.";
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     msg);
  }
}

void NeedlemanWunsch::setMatrix_(const NeedlemanWunsch::ScoringMatrix& matrix)
{
  if (matrix == ScoringMatrix::identityMatrix)
  {
    matrixPtr_ = &adaptedIdentity;
  }

  else if (matrix == ScoringMatrix::PAM30MSMatrix)
  {
    matrixPtr_ = &PAM30MS;
  }

  else
  {
    String msg = "Matrix is not known! Valid choices are: "
                                       "'identityMatrix', 'PAM30MSMatrix'.";
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     msg);
  }
}

void NeedlemanWunsch::setPenalty_(const int& penalty)
{

    if (penalty >= 0)
    {
      String msg = "Gap penalty should be negative";
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       msg);
    }
      gapPenalty_ = penalty;
}

NeedlemanWunsch::ScoringMatrix NeedlemanWunsch::getMatrix_() const
{
    if (*matrixPtr_ == adaptedIdentity)
    {
      return ScoringMatrix::identityMatrix;
    }
    else
    {
      return ScoringMatrix::PAM30MSMatrix;
    }
}

int NeedlemanWunsch::getPenalty_() const
{
    return gapPenalty_;
}

int NeedlemanWunsch::getIndex_(const char& a, const char& b) const //noch optimieren (Tina)
{
vector<pair<char,int>> vec =
    {
        {'A', 0}, {'R', 1}, {'N', 2},
        {'D', 3}, {'C', 4}, {'Q', 5},
        {'E', 6}, {'G', 7}, {'H', 8},
        {'I', 9}, {'L', 10}, {'K', 11},
        {'M', 12}, {'F', 13}, {'P', 14},
        {'S', 15}, {'T', 16}, {'W', 17},
        {'Y', 18}, {'V', 19}, {'B', 20},
        {'Z', 21}, {'X', 22}, {'*', 23}
    };
int x = -1;
int y = -1;
for (int i = 0; i < vec.size(); ++i)
{
if (vec[i].first == a)
x = vec[i].second;
if (vec[i].first == b)
y = vec[i].second;
}
if (x == -1)
x = 23;
if (y == -1)
y = 23;
return x + y*vec.size();
}
/*
double NeedlemanWunsch::align_(const String& seq1, const String& seq2) //vollständige matrix
{
  seq1len_ = seq1.length();
  seq2len_ = seq2.length();

  vector<int> matrix((seq1len_+1)*(seq2len_+1), 0);//matrix mit 0en initialisieren
  for (unsigned i = 1; i <= seq1len_; ++i) //vertikale mit gapkkosten initialisieren
    matrix[i*(seq2len_+1)]=i*gapPenalty_;
  for (unsigned i =0; i<=seq2len_;++i)//horizontale mit gapkosten initialieren
    matrix[i]=i*gapPenalty_;
  for (unsigned i=1;i<=seq1len_;++i)
  {
    for (unsigned j=1;j<=seq2len_;++j)
    {
      matrix[i*(seq2len_ +1)+j]=max(max((matrix[i*(seq2len_+1)+j-1]+gapPenalty_), (matrix[(i-1)*(seq2len_+1)+j]+gapPenalty_)),
                                    (matrix[(i-1)*(seq2len_+1)+j-1])+ (*matrixPtr_)[getIndex_(seq1[i-1], seq2[j-1])]);
    }
  }
  cout<<matrix[(seq1len_+1)*(seq2len_+1)-1]<<endl;
  return matrix[(seq1len_+1)*(seq2len_+1)-1];
}
 */

//linear space (2 Zeilen)
 double NeedlemanWunsch::align_(const String& seq1, const String& seq2)
 {
   seq1len_ = seq1.length();
   seq2len_ = seq2.length();

   vector<int> firstRow{};
   vector<int> secondRow(seq2len_+1,0);
   vector<int>* firstRowPtr = &firstRow;
   vector<int>* secondRowPtr = &secondRow;


   for (unsigned i = 0; i <= seq2len_; ++i)//horizontale mit gapkosten initialieren
   {
     firstRow.push_back(i * gapPenalty_);
   }

   for (unsigned i = 1;i <= seq1len_; ++i) //second row berechnen und swappen
   {
     (*secondRowPtr)[0] = i * gapPenalty_; //erster wert in der zeile mit gapkosten
     for (unsigned j = 1; j <= seq2len_; ++j) //secondRow berechnen
     {
       (*secondRowPtr)[j] = (max(max(((*secondRowPtr)[j-1] + gapPenalty_), ((*firstRowPtr)[j] + gapPenalty_)),
                                    ((*firstRowPtr)[j-1]) + (*matrixPtr_)[getIndex_(seq1[i-1], seq2[j-1])]));//statt getIndex: [seq1[i-1] - 'A'] [seq2[j-1] - 'A'] und matrix entsprechend aufbauen
                                    //cout<<(*matrixPtr_)[getIndex_(seq1[i-1], seq2[j-1])]<<endl;
     }
     swap(firstRowPtr, secondRowPtr);
   }
   cout<<(*firstRowPtr)[seq2len_]<<endl;
   return (*firstRowPtr)[seq2len_];
 }

}