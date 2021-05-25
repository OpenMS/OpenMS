#include <OpenMS/ANALYSIS/SEQUENCE/NeedlemanWunsch.h>
#include <iostream>
#include <OpenMS/CONCEPT/Exception.h>
#include  <utility> //swap


using namespace std;
namespace OpenMS
{
    static int adaptedIdentity[26][26]
            {

                    //         A  B  C  D  E  F  G  H  I      J     K  L  M  N      O     P  Q  R  S  T      U     V  W  X  Y  Z
                    /* A */   {1, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* B */   {0, 1, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* C */   {0, 0, 1, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* D */   {0, 0, 0, 1, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* E */   {0, 0, 0, 0, 1, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* F */   {0, 0, 0, 0, 0, 1, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* G */   {0, 0, 0, 0, 0, 0, 1, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* H */   {0, 0, 0, 0, 0, 0, 0, 1, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* I */   {0, 0, 0, 0, 0, 0, 0, 0, 1, INT8_MAX, 0, 1, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* J */   {INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX,},
                    /* K */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 1, 0, 0, 0, INT8_MAX, 0, 1, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* L */   {0, 0, 0, 0, 0, 0, 0, 0, 1, INT8_MAX, 0, 1, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* M */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 1, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* N */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 1, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* O */   {INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX,},
                    /* P */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 1, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* Q */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 1, 0, 0, 0, INT8_MAX, 0, 1, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* R */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 1, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* S */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 1, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* T */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 1, INT8_MAX, 0, 0, 0, 0, 0},
                    /* U */   {INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX,},
                    /* V */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 1, 0, 0, 0, 0},
                    /* W */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 1, 0, 0, 0},
                    /* X */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0},
                    /* Y */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 1, 0},
                    /* Z */   {0, 0, 0, 0, 0, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, INT8_MAX, 0, 0, 0, 0, 0, INT8_MAX,0, 0, 0, 0, 1}
            };


    static int PAM30MS[26][26]
            {

        //             A    B   C   D   E   F   G   H   I      J      K   L   M   N     O       P   Q   R   S   T      U      V   W  X   Y   Z
        /* A */        {6, -7, -6, -3, -2, -8, -2, -7, -5, INT8_MAX, -7, -6, -5, -4, INT8_MAX, -2, -4, -7,  0, -1, INT8_MAX, -2,-13, 0, -8, -6},
        /* B */       {-7,  5,-11, -7, -7,-12, -8, -4, -6, INT8_MAX,  5, -7, -3, -4, INT8_MAX, -5, -3,  5, -4, -5, INT8_MAX, -9, -7, 0,-10,  1},
        /* C */       {-6,-11, 10,-14,-14,-13, -9, -7, -6, INT8_MAX,-14,-11,-13,-11, INT8_MAX, -8,-14, -8, -3, -8, INT8_MAX, -6,-15, 0, -4,-14},
        /* D */       {-3, -7,-14,  8,  2,-15, -3, -4, -7, INT8_MAX, -4,-10,-11,  2, INT8_MAX, -8, -2,-10, -4, -5, INT8_MAX, -8,-15, 0,-11, -3},
        /* E */       {-2, -7,-14,  2,  8,-14, -4, -5, -5, INT8_MAX, -4, -7, -7, -2, INT8_MAX, -5,  1, -9, -4, -6, INT8_MAX, -6,-17, 0, -8, -2},
        /* F */       {-8,-12,-13,-15,-14,  9, -9, -6, -2, INT8_MAX,-14, -3, -4, -9, INT8_MAX,-10,-13, -9, -6, -9, INT8_MAX, -8, -4, 0,  2,-14},
        /* G */       {-2, -8, -9, -3, -4, -9,  6, -9,-11, INT8_MAX, -7,-11, -8, -3, INT8_MAX, -6, -7, -9, -2, -6, INT8_MAX, -5,-15, 0,-14, -7},
        /* H */       {-7, -4, -7, -4, -5, -6, -9,  9, -9, INT8_MAX, -6, -8,-10,  0, INT8_MAX, -4,  1, -2, -6, -7, INT8_MAX, -6, -7, 0, -3, -3},
        /* I */       {-5, -6, -6, -7, -5, -2,-11, -9,  8, INT8_MAX, -6,  5, -1, -5, INT8_MAX, -8, -8, -5, -7, -2, INT8_MAX,  2,-14, 0, -6, -7},
        /* J */       {INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX,},
        /* K */       {-7,  5,-14, -4, -4,-14, -7, -6, -6, INT8_MAX,  7, -7, -2, -1, INT8_MAX, -6, -3,  0, -4, -3, INT8_MAX, -9,-12, 0, -9,  4},
        /* L */       {-6, -7,-11,-10, -7, -3,-11, -8,  5, INT8_MAX, -7,  5,  0, -6, INT8_MAX, -8, -7, -7, -8, -5, INT8_MAX,  0,-10, 0, -7, -7},
        /* M */       {-5, -3,-13,-11, -7, -4, -8,-10, -1, INT8_MAX, -2,  0, 11, -9, INT8_MAX, -8, -4, -4, -5, -4, INT8_MAX, -1,-13, 0,-11, -3},
        /* N */       {-4, -4,-11,  2, -2, -9, -3,  0, -5, INT8_MAX, -1, -6, -9,  8, INT8_MAX, -6, -3, -6,  0, -2, INT8_MAX, -8, -8, 0, -4, -2},
        /* O */       {INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX,},
        /* P */       {-2, -5, -8, -8, -5,-10, -6, -4, -8, INT8_MAX, -6, -8, -8, -6, INT8_MAX,  8, -3, -4, -2, -4, INT8_MAX, -6,-14, 0,-13, -5},
        /* Q */       {-4, -3,-14, -2,  1,-13, -7,  1, -8, INT8_MAX, -3, -7, -4, -3, INT8_MAX, -3,  8, -2, -5, -5, INT8_MAX, -7,-13, 0,-12,  4},
        /* R */       {-7,  5, -8,-10, -9, -9, -9, -2, -5, INT8_MAX,  0, -7, -4, -6, INT8_MAX, -4, -2,  8, -3, -6, INT8_MAX, -8, -2, 0, 10, -1},
        /* S */       {0,  -4, -3, -4, -4, -6, -2, -6, -7, INT8_MAX, -4, -8, -5,  0, INT8_MAX, -2, -5, -3,  6,  0, INT8_MAX, -6, -5, 0, -7, -5},
        /* T */       {-1, -5, -8, -5, -6, -9, -6, -7, -2, INT8_MAX, -3, -5, -4, -2, INT8_MAX, -4, -5, -6,  0,  7, INT8_MAX, -3,-13, 0, -6, -4},
        /* I */       {INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX,},
        /* V */       {-2, -9, -6, -8, -6, -8, -5, -6,  2, INT8_MAX, -9,  0, -1, -8, INT8_MAX, -6, -7, -8, -6, -3, INT8_MAX,  7,-15, 0, -7, -8},
        /* W */       {-13,-7,-15,-15,-17, -4,-15, -7,-14, INT8_MAX,-12,-10,-13, -8, INT8_MAX,-14,-13, -2, -5,-13, INT8_MAX,-15, 13, 0, -5,-13},
        /* X */       {0,   0,  0,  0,  0,  0,  0,  0,  0, INT8_MAX,  0,  0,  0,  0, INT8_MAX,  0,  0,  0,  0,  0, INT8_MAX,  0,  0, 0,  0,  0},
        /* Y */       {-8,-10, -4,-11, -8,  2,-14, -3, -6, INT8_MAX, -9, -7,-11, -4, INT8_MAX,-13,-12,-10, -7, -6, INT8_MAX, -7, -5, 0, 10,-11},
        /* Z */       {-6,  1,-14, -3, -2,-14, -7, -3, -7, INT8_MAX,  4, -7, -3, -2, INT8_MAX, -5,  4, -1, -5, -4, INT8_MAX, -8,-13, 0,-11,  4}

            };



NeedlemanWunsch::NeedlemanWunsch(NeedlemanWunsch::ScoringMatrix matrix, int penalty)
{
  setMatrix(matrix);
  setPenalty(penalty);

}

void NeedlemanWunsch::setMatrix(const NeedlemanWunsch::ScoringMatrix& matrix)
{
  if (matrix == ScoringMatrix::identity)
  {
    matrixPtr_ = &adaptedIdentity;
  }

  else if (matrix == ScoringMatrix::PAM30MS)
  {
    matrixPtr_ = &PAM30MS;
  }
}

void NeedlemanWunsch::setMatrix(const std::string& matrix)
{
  auto first = &validMatrices_[0];
  auto last = &validMatrices_[2];
  const auto it = std::find(first, last, matrix);
  if (it == last)
  {
    String msg = "Matrix is not known! Valid choices are: "
                 "'identity', 'PAM30MS'.";
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     msg);
  }
  setMatrix(static_cast<NeedlemanWunsch::ScoringMatrix>(it - first));
}



void NeedlemanWunsch::setPenalty(const int& penalty)
{
  gapPenalty_ = penalty;
}

NeedlemanWunsch::ScoringMatrix NeedlemanWunsch::getMatrix() const
{
    if (*matrixPtr_ == adaptedIdentity)
    {
      return ScoringMatrix::identity;
    }
    else
    {
      return ScoringMatrix::PAM30MS;
    }
}

int NeedlemanWunsch::getPenalty() const
{
    return gapPenalty_;
}
/*
int NeedlemanWunsch::align(const String& seq1, const String& seq2) //vollständige matrix
{
  seq1len_ = seq1.length();
  seq2len_ = seq2.length();

  vector<int> matrix((seq1len_+1)*(seq2len_+1), 0);//matrix mit 0en initialisieren
  for (unsigned i = 1; i <= seq1len_; ++i) //vertikale mit gapkkosten initialisieren
    matrix[i*(seq2len_+1)]=i*(-gapPenalty_);
  for (unsigned i =0; i<=seq2len_;++i)//horizontale mit gapkosten initialieren
    matrix[i]=i*(-gapPenalty_);
  for (unsigned i=1;i<=seq1len_;++i)
  {
    for (unsigned j=1;j<=seq2len_;++j)
    {
      matrix[i*(seq2len_ +1)+j]=max(max((matrix[i*(seq2len_+1)+j-1]-gapPenalty_), (matrix[(i-1)*(seq2len_+1)+j]-gapPenalty_)),
                                    (matrix[(i-1)*(seq2len_+1)+j-1])+ (*matrixPtr_)[seq1[i-1] - 'A'] [seq2[j-1] - 'A']);
    }
  }
  return matrix[(seq1len_+1)*(seq2len_+1)-1];
}
*/
/*
//linear space (2 Zeilen) //seit vectoren member sind: munmap_chunk(): invalid pointer
 int NeedlemanWunsch::align(const String& seq1, const String& seq2)
 {
   seq1len_ = seq1.length();
   seq2len_ = seq2.length();

   firstRow_.resize(seq1len_);
   secondRow_.resize(seq2len_+1);
   vector<int>* firstRowPtr = &firstRow_;
   vector<int>* secondRowPtr = &secondRow_;


   for (unsigned i = 0; i <= seq2len_; ++i)//horizontale mit gapkosten initialieren
   {
     firstRow_[i] = i * ((-1 * gapPenalty_));
   }

   for (unsigned i = 1;i <= seq1len_; ++i) //second row berechnen und swappen
   {
     (*secondRowPtr)[0] = i * ((-1 * gapPenalty_)); //erster wert in der zeile mit gapkosten
     for (unsigned j = 1; j <= seq2len_; ++j) //secondRow berechnen
     {
       (*secondRowPtr)[j] = (max(max(((*secondRowPtr)[j-1] - gapPenalty_), ((*firstRowPtr)[j] - gapPenalty_)),
                                    ((*firstRowPtr)[j-1]) + (*matrixPtr_)[seq1[i-1] - 'A'] [seq2[j-1] - 'A']));//[getIndex_(seq1[i-1], seq2[j-1])]));//statt getIndex: [seq1[i-1] - 'A'] [seq2[j-1] - 'A'] und matrix entsprechend aufbauen
     }
     swap(firstRowPtr, secondRowPtr);
   }
   return (*firstRowPtr)[seq2len_];
 }
 */


  int NeedlemanWunsch::align(const String& seq1, const String& seq2)
  {
    seq1len_ = seq1.length();
    seq2len_ = seq2.length();

    firstRow_.resize(seq2len_+1); // both rows have the same length
    secondRow_.resize(seq2len_+1);

    int* firstRowPtr = &(firstRow_[0]);
    int* secondRowPtr = &(secondRow_[0]);


    for (unsigned i = 0; i <= seq2len_; ++i)//horizontale mit gapkosten initialieren
    {
      firstRow_[i] = i * (-gapPenalty_);
    }

    for (unsigned i = 1;i <= seq1len_; ++i) //second row berechnen und swappen
    {
      (*secondRowPtr) = i * (-gapPenalty_); //erster wert in der zeile mit gapkosten //second row pointer muss auf die erste stelle zeigen
      for (unsigned j = 1; j <= seq2len_; ++j) //secondRow berechnen
      {
        (*(secondRowPtr+j)) = max(max(((*(secondRowPtr+j-1)) - gapPenalty_), ((*(firstRowPtr+j)) - gapPenalty_)),
                                  ((*(firstRowPtr+j-1)) + (*matrixPtr_)[seq1[i-1] - 'A'] [seq2[j-1] - 'A']));
      }
      swap(firstRowPtr, secondRowPtr);
    }
    return (*(firstRowPtr+seq2len_));
  }

}