/*
FPMSComp: mass spectrometry data compression

Ruochen Yang
email: rcyang624@126.com
*/

#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <cmath>
#include <bitset>
#include <stdint.h>
#include <string.h>
#include <string>
#include <iostream>
#include <malloc.h>
#include "base64.h"
#include "tinyxml2.h"
#include <cstring>

//using namespace std;
using namespace tinyxml2;


namespace Comp
{
	void pairsComp(FILE * fpW,XMLElement * peaks);

	int match_type;
	// mz
	int * this_mz_hex;
	int * last_mz_hex;
	int * zero_len;
	int * diff_reshape;
	uint8_t * diff_reshape_byte;
	int diff_reshape_len;
	int diff[8];
	int cnt;
	bool flag_cnt_zero;
	// intensity
	int * this_intensity_hex;
	int * find_intensity_hex;
	int * pointer;
	unsigned char * residual;
	unsigned char * unmatch_data;
	bool flag_find_match;
	// arithmetic compression
	uint8_t * output_bitstream;
	int * output_len;
	int * cum_cnt;

	void arithmeticEncoder(int * ori,int ori_len,int range,int * cum_cnt,
		uint8_t * output_bitstream,int * output_len);
	void sendBit(int type,uint8_t * buffer,int * buffer_len,uint8_t * output_bitstream,int * output_len);

}

namespace DeComp
{
	void pairsDecomp(FILE ** fp, XMLElement * scan, XMLDocument * doc);
	
	int pairs_len;
	int range;
	int * cum_cnt;
	int input_len;
	uint8_t * input_bitstream;
	int input_len2;
	uint8_t * input_bitstream2;
	int diff_reshape_len;
	int match_len;
	uint8_t * residual;
	int unmatch_len;
	uint8_t * unmatch_data;
	int * zero_len;
	int * diff_reshape_hex;
	int * last_mz_hex;
	int * this_mz_hex;
	int * diff_hex;
	unsigned char *  bin;
	int pos;
	int match_type;
	int * pointer ;
	int unmatch_cnt;
	int match_cnt;
	int * find_int_hex;
	int * this_int_hex;
	int flen;
	char * base64char;

	void arithmeticDecoder(uint8_t * input_bitstream,int input_len,
		int range,int * cum_cnt,int * ori,int ori_len);
	void getTag(uint8_t * input_bitstream,int * byte_pos,int * bit_pos,uint64_t * tag);
}
