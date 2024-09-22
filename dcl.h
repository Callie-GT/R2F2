#ifndef __DCL_H__
#define __DCL_H__

using namespace std;

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <random>
#include <iostream>
#include <fstream>

#include "ap_fixed.h"
#include "hls_half.h"

// standard 32-bit floating point bitwidth
#define std_FP32_TOTAL_BIT      32
#define std_FP32_SIGN_BIT       1
#define std_FP32_EXPN_BIT       8
#define std_FP32_MANT_BIT       23
#define std_FP32_MANT_MUL_BIT   48

#define std_FP32_BIAS           127

// used in templated floating point multiplication
// just to store values for simplicity; the computations use the actual bitwidth
#define myFP_MAX_TOTAL_BIT      16
#define myFP_MAX_SIGN_BIT       1
#define myFP_MAX_EXPN_BIT       6 // in testing, exp can go up to 6, although standard is 5-bit
#define myFP_MAX_MANT_BIT       13 // in testing, mant can go up to 13, although standard is 10-bit
#define myFP_MAX_MANT_MUL_BIT   28
// note: use myFP_MAX_EXPN_BIT = 5, myFP_MAX_MANT_BIT = 10, and myFP_MAX_MANT_MUL_BIT = 22 for synthesis

// bitwidth of the flexible run-time reconfigurable floating point
#define TOTAL_FLEX  16  // including the sign bit
#define FIX_EXPN    3
#define FIX_MANT    9
#define FLEX_BIT    3
#define MANT_APPX   3   // how many bits will be kept (approximation) after FIX_MANT


typedef ap_uint<myFP_MAX_TOTAL_BIT> myFP;
typedef ap_uint<myFP_MAX_SIGN_BIT> myFP_sign;
typedef ap_uint<myFP_MAX_EXPN_BIT + 2> myFP_expn; // account for exp overflow
typedef ap_uint<myFP_MAX_MANT_BIT + 1> myFP_mant; // account for the implicit 1 before mantissa

#define DATA_PAIRS  10000

// multiplication functions
float templated_float_multiply(float a, float b, bool *overflow, int EB, int MB);
float flexible_float_multiply(float a, float b, bool *overflow, int EB, int MB, bool adjust);
float flexible_float_multiply_retry(float a, float b, bool adjust);

void float2myFP(float *f, myFP *h, int EB, int MB, bool *overflow);
float myFP2float(const myFP h, int EB, int MB);

void printFloatBits(float f);
void printHalfBits(half f);
bool compareFloatsBitwise(float a, float b);
float float_multiply_32std(float a, float b, bool *overflow);
float float_multiply_16std(float a, float b, bool *overflow);
std::string floatToBinary(half num);


float randomFloat32();
float randomFloat16();
float randomFloatRange(float min, float max);
void read_operands(string file_name, float* values);

#endif