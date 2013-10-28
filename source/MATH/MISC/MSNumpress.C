/*
        MSNumpress.cpp
        johan.teleman@immun.lth.se
        
        This distribution goes under the BSD 3-clause license. If you prefer to use Apache
        version 2.0, that is also available at https://github.com/fickludd/ms-numpress
        Copyright (c) 2013, Johan Teleman
        All rights reserved.

        Redistribution and use in source and binary forms, with or without modification,
        are permitted provided that the following conditions are met:

*         Redistributions of source code must retain the above copyright notice, this list
        of conditions and the following disclaimer.
*        Redistributions in binary form must reproduce the above copyright notice, this
        list of conditions and the following disclaimer in the documentation and/or other
        materials provided with the distribution.
*        Neither the name of the Lund University nor the names of its contributors may be
        used to endorse or promote products derived from this software without specific
        prior written permission.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
        EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
        OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
        SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
        OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
        HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <iostream>
#include <cmath>
#include <OpenMS/MATH/MISC/MSNumpress.h>

namespace ms {
namespace numpress {
namespace MSNumpress {

using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;


const int ONE = 1;
bool is_big_endian() {
	return *((char*)&(ONE)) == 1;
}
bool IS_BIG_ENDIAN = is_big_endian();



/////////////////////////////////////////////////////////////

void encodeFixedPoint(
		double fixedPoint, 
		unsigned char *result
) {
	int i;
	unsigned char *fp = (unsigned char*)&fixedPoint;
	for (i=0; i<8; i++) {
		result[i] = fp[IS_BIG_ENDIAN ? (7-i) : i];
	}
}



double decodeFixedPoint(
		const unsigned char *data
) {
	int i;
	double fixedPoint;
	unsigned char *fp = (unsigned char*)&fixedPoint;
		
	for (i=0; i<8; i++) {
		fp[i] = data[IS_BIG_ENDIAN ? (7-i) : i];
	}
	
	return fixedPoint;
}

/////////////////////////////////////////////////////////////

/**
 * Encodes the int x as a number of halfbytes in res. 
 * res_length is incremented by the number of halfbytes, 
 * which will be 1 <= n <= 9
 */
void encodeInt(
		const int x,
		unsigned char* res,
		size_t *res_length	
) {
	int i, l, m;
	int mask = 0xf0000000;
	int init = x & mask;

	if (init == 0) {
		l = 8;
		for (i=0; i<8; i++) {
			m = mask >> (4*i);
			if ((x & m) != 0) {
				l = i;
				break;
			}
		}
		res[0] = l;
		for (i=l; i<8; i++) {
			res[1+i-l] = x >> (4*(i-l));
		}
		*res_length += 1+8-l;

	} else if (init == mask) {
		l = 7;
		for (i=0; i<8; i++) {
			m = mask >> (4*i);
			if ((x & m) != m) {
				l = i;
				break;
			}
		}
		res[0] = l + 8;
		for (i=l; i<8; i++) {
			res[1+i-l] = x >> (4*(i-l));
		}
		*res_length += 1+8-l;

	} else {
		res[0] = 0;
		for (i=0; i<8; i++) {
			res[1+i] = x >> (4*i);
		}
		*res_length += 9;

	}
}



/**
 * Decodes an int from the half bytes in bp. Lossless reverse of encodeInt 
 */
void decodeInt(
		const unsigned char *data,
		size_t *di,
		int *half,
		int *res
) {
	size_t n;
	size_t i;
	int mask, m;
	unsigned char head;
	unsigned char hb;

	if (*half == 0) {
		head = data[*di] >> 4;
	} else {
		head = data[*di] & 0xf;
		(*di)++;
	}

	*half = 1-(*half);
	*res = 0;
	
	if (head <= 8) {
		n = head;
	} else { // leading ones, fill n half bytes in res
		n = head - 8;
		mask = 0xf0000000;
		for (i=0; i<n; i++) {
			m = mask >> (4*i);
			*res = *res | m;
		}
	}
	
	if (n == 8) {
		return;
	}
	
	for (i=n; i<8; i++) {
		if (*half == 0) {
			hb = data[*di] >> 4;
		} else {
			hb = data[*di] & 0xf;
			(*di)++;
		}
		*res = *res | (hb << ((i-n)*4));
		*half = 1 - (*half);
	}
}




/////////////////////////////////////////////////////////////


double optimalLinearFixedPoint(
		const double *data, 
		size_t dataSize
) {
	/*
	 * safer impl - apparently not needed though
	 *
	if (dataSize == 0) return 0;
	
	double maxDouble = 0;
	double x;

	for (size_t i=0; i<dataSize; i++) {
		x = data[i];
		maxDouble = max(maxDouble, x);
	}

	return floor(0xFFFFFFFF / maxDouble);
	*/
	if (dataSize == 0) return 0;
	if (dataSize == 1) return floor(0x7FFFFFFFl / data[0]);
	double maxDouble = max(data[0], data[1]);
	double extrapol;
	double diff;

	for (size_t i=2; i<dataSize; i++) {
		extrapol = data[i-1] + (data[i-1] - data[i-2]);
		diff = data[i] - extrapol;
		maxDouble = max(maxDouble, ceil(abs(diff)+1));
	}

	return floor(0x7FFFFFFFl / maxDouble);
}

size_t encodeLinear(
		const double *data, 
		size_t dataSize, 
		unsigned char *result,
		double fixedPoint
) {
	unsigned long long ints[3];
	size_t i, ri;
	unsigned char halfBytes[10];
	size_t halfByteCount;
	size_t hbi;
	long long extrapol;
	int diff;

	//printf("Encoding %d doubles with fixed point %f\n", (int)dataSize, fixedPoint);
	encodeFixedPoint(fixedPoint, result);


	if (dataSize == 0) return 8;

	ints[1] = data[0] * fixedPoint + 0.5;
	for (i=0; i<4; i++) {
		result[8+i] = (ints[1] >> (i*8)) & 0xff;
	}

	if (dataSize == 1) return 12;

	ints[2] = data[1] * fixedPoint + 0.5;
	for (i=0; i<4; i++) {
		result[12+i] = (ints[2] >> (i*8)) & 0xff;
	}

	halfByteCount = 0;
	ri = 16;

	for (i=2; i<dataSize; i++) {
		ints[0] = ints[1];
		ints[1] = ints[2];
		ints[2] = data[i] * fixedPoint + 0.5;
		extrapol = ints[1] + (ints[1] - ints[0]);
		diff = ints[2] - extrapol;
		//printf("%lu %lu %lu,   extrapol: %ld    diff: %d \n", ints[0], ints[1], ints[2], extrapol, diff);
		encodeInt(diff, &halfBytes[halfByteCount], &halfByteCount);
		/*
		printf("%d (%d):  ", diff, (int)halfByteCount);
		for (size_t j=0; j<halfByteCount; j++) {
			printf("%x ", halfBytes[j] & 0xf);
		}
		printf("\n");
		*/
		
		
		for (hbi=1; hbi < halfByteCount; hbi+=2) {
			result[ri] = (halfBytes[hbi-1] << 4) | (halfBytes[hbi] & 0xf);
			//printf("%x \n", result[ri]);
			ri++;
		}
		if (halfByteCount % 2 != 0) {
			halfBytes[0] = halfBytes[halfByteCount-1];
			halfByteCount = 1;
		} else {
			halfByteCount = 0;
		}
	}
	if (halfByteCount == 1) {
		result[ri] = halfBytes[0] << 4;
		ri++;
	}
	return ri;
}


size_t decodeLinear(
		const unsigned char *data,
		const size_t dataSize,
		double *result
) {
	size_t i;
	size_t ri = 0;
	unsigned int init;
	int diff;
	long long ints[3];
	//double d;
	size_t di;
	int half;
	long long extrapol;
	long long y;
	double fixedPoint;
	
	//printf("Decoding %d bytes with fixed point %f\n", (int)dataSize, fixedPoint);

	if (dataSize < 8) return -1;
	
	fixedPoint = decodeFixedPoint(data);


	if (dataSize < 12) return -1;
	
	try {
		ints[1] = 0;
		for (i=0; i<4; i++) {
			ints[1] = ints[1] | ((0xff & (init = data[8+i])) << (i*8));
		}
		result[0] = ints[1] / fixedPoint;

		if (dataSize == 12) return 1;
		if (dataSize < 16) return -1;

		ints[2] = 0;
		for (i=0; i<4; i++) {
			ints[2] = ints[2] | ((0xff & (init = data[12+i])) << (i*8));
		}
		result[1] = ints[2] / fixedPoint;
			
		half = 0;
		ri = 2;
		di = 16;
		
		while (di < dataSize) {
			ints[0] = ints[1];
			ints[1] = ints[2];
			if (di == (dataSize - 1) && half == 1) {
				if ((data[di] & 0xf) != 0x8) {
					break;
				}
			}
			decodeInt(data, &di, &half, &diff);
			
			extrapol = ints[1] + (ints[1] - ints[0]);
			y = extrapol + diff;
			//printf("%lu %lu,   extrapol: %ld    diff: %d \n", ints[0], ints[1], extrapol, diff);
			result[ri++] 	= y / fixedPoint;
			ints[2] 		= y;
		}
	} catch (...) {
		cerr << "DECODE ERROR" << endl;
		cerr << "i: " << i << endl;
		cerr << "ri: " << ri << endl;
		cerr << "di: " << di << endl;
		cerr << "half: " << half << endl;
		cerr << "dataSize: " << dataSize << endl;
		cerr << "ints[]: " << ints[0] << ", " << ints[1] << ", " << ints[2] << endl;
		cerr << "extrapol: " << extrapol << endl;
		cerr << "y: " << y << endl;

		for (i = di - 3; i < min(di + 3, dataSize); i++) {
			cerr << "data[" << i << "] = " << data[i];
		}
		cerr << endl;
	}
	
	return ri;
}

void encodeLinear(
		const std::vector<double> &data, 
		std::vector<unsigned char> &result,
		double fixedPoint
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 5);
	size_t encodedLength = encodeLinear(&data[0], dataSize, &result[0], fixedPoint);
	result.resize(encodedLength);
}

void decodeLinear(
		const std::vector<unsigned char> &data,
		std::vector<double> &result
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 2);
	size_t decodedLength = decodeLinear(&data[0], dataSize, &result[0]);
	result.resize(decodedLength);
}

/////////////////////////////////////////////////////////////


size_t encodePic(
		const double *data, 
		size_t dataSize, 
		unsigned char *result
) {
	size_t i, ri, count;
	unsigned char halfBytes[10];
	size_t halfByteCount;
	size_t hbi;

	//printf("Encoding %d doubles\n", (int)dataSize);

	halfByteCount = 0;
	ri = 0;

	for (i=0; i<dataSize; i++) {
		count = data[i] + 0.5;
		//printf("%d %d %d,   extrapol: %d    diff: %d \n", ints[0], ints[1], ints[2], extrapol, diff);
		encodeInt(count, &halfBytes[halfByteCount], &halfByteCount);
		/*
		printf("%d (%d):  ", count, (int)halfByteCount);
		for (j=0; j<halfByteCount; j++) {
			printf("%x ", halfBytes[j] & 0xf);
		}
		printf("\n");
		*/
		
		for (hbi=1; hbi < halfByteCount; hbi+=2) {
			result[ri] = (halfBytes[hbi-1] << 4) | (halfBytes[hbi] & 0xf);
			//printf("%x \n", result[ri]);
			ri++;
		}
		if (halfByteCount % 2 != 0) {
			halfBytes[0] = halfBytes[halfByteCount-1];
			halfByteCount = 1;
		} else {
			halfByteCount = 0;
		}
	}
	if (halfByteCount == 1) {
		result[ri] = halfBytes[0] << 4;
		ri++;
	}
	return ri;
}





size_t decodePic(
		const unsigned char *data,
		const size_t dataSize,
		double *result
) {
	size_t i, ri;
	int count;
	//double d;
	size_t di;
	int half;

	try {
		half = 0;
		ri = 0;
		di = 0;
		
		while (di < dataSize) {
			if (di == (dataSize - 1) && half == 1) {
				if ((data[di] & 0xf) != 0x8) {
					break;
				}
			}
			decodeInt(&data[0], &di, &half, &count);
			
			//printf("count: %d \n", count);
			result[ri++] 	= count;
		}
	} catch (...) {
		cerr << "DECODE ERROR" << endl;
		cerr << "ri: " << ri << endl;
		cerr << "di: " << di << endl;
		cerr << "half: " << half << endl;
		cerr << "dataSize: " << dataSize << endl;
		cerr << "count: " << count << endl;

		for (i = di - 3; i < min(di + 3, dataSize); i++) {
			cerr << "data[" << i << "] = " << data[i];
		}
		cerr << endl;
	}
	return ri;
}


void encodePic(
		const std::vector<double> &data,  
		std::vector<unsigned char> &result
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 5);
	size_t encodedLength = encodePic(&data[0], dataSize, &result[0]);
	result.resize(encodedLength);
}



void decodePic(
		const std::vector<unsigned char> &data,  
		std::vector<double> &result
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 2);
	size_t decodedLength = decodePic(&data[0], dataSize, &result[0]);
	result.resize(decodedLength);
}

/////////////////////////////////////////////////////////////

double optimalSlofFixedPoint(
		const double *data, 
		size_t dataSize
) {
	if (dataSize == 0) return 0;
	
	double maxDouble = 1;
	double x;
	double fp;

	for (size_t i=0; i<dataSize; i++) {
		x = log(data[i]+1);
		maxDouble = max(maxDouble, x);
	}

	fp = floor(0xFFFF / maxDouble);

	//cout << "    max val: " << maxDouble << endl;
	//cout << "fixed point: " << fp << endl;

	return fp;
}

size_t encodeSlof(
		const double *data, 
		size_t dataSize, 
		unsigned char *result,
		double fixedPoint
) {
	size_t i, ri;
	unsigned short x;
	encodeFixedPoint(fixedPoint, result);

	ri = 8;
	for (i=0; i<dataSize; i++) {
		x = log(data[i]+1) * fixedPoint + 0.5;
		result[ri++] = x & 0xff;
		result[ri++] = x >> 8;
	}
	return ri;
}



size_t decodeSlof(
		const unsigned char *data, 
		const size_t dataSize, 
		double *result
) {
	size_t i, ri;
	unsigned short x;
	ri = 0;
	double fixedPoint;

	if (dataSize < 8) return -1;
	fixedPoint = decodeFixedPoint(data);

	for (i=8; i<dataSize; i+=2) {
		x = data[i] | (data[i+1] << 8);
		result[ri++] = exp(x / fixedPoint) - 1;
	}
	return ri;
}

void encodeSlof(
		const std::vector<double> &data,  
		std::vector<unsigned char> &result,
		double fixedPoint
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 2);
	size_t encodedLength = encodeSlof(&data[0], dataSize, &result[0], fixedPoint);
	result.resize(encodedLength);
}



void decodeSlof(
		const std::vector<unsigned char> &data,  
		std::vector<double> &result
) {
	size_t dataSize = data.size();
	result.resize(dataSize / 2);
	size_t decodedLength = decodeSlof(&data[0], dataSize, &result[0]);
	result.resize(decodedLength);
}

}
} // namespace numpress
} // namespace ms
