#ifndef __HDR__
#define __HDR__ 1

#if defined(__cplusplus)
extern "C"{
#endif   //(__cplusplus)

#include "src/DFT.h"

//for test
//#define HDR_DEFAULT_CROSSOVER_MIN    -15;
//#define HDR_DEFAULT_CROSSOVER_MAX    -3;

//for bach
#define HDR_DEFAULT_CROSSOVER_MIN    -9;
#define HDR_DEFAULT_CROSSOVER_MAX    -3;

#define HDR_DEFAULT_RELEASE_COEFF    0.999;

/*--------------------------------------------------------------------*/
typedef struct Opaque_HDR_Struct HDR;

typedef void (*hdr_onprocess_t)(void* onprocess_self, dft_sample_t* real, int N);

HDR* hdr_new(int window_size);
HDR* hdr_destroy(HDR* self);
void hdr_process(HDR* self, dft_sample_t* high_gain_input, dft_sample_t* low_gain_input, int len, double sample_rate, hdr_onprocess_t onprocess, void* onprocess_self);

/* get the last little bit out of the output buffer at the end of the stream */
void hdr_flush(HDR* self, hdr_onprocess_t onprocess, void* onprocess_self);

void   hdr_set_crossover_min(HDR* self, double coeff_dB);
double hdr_get_crossover_min(HDR* self);

void   hdr_set_crossover_max(HDR* self, double coeff_dB);
double hdr_get_crossover_max(HDR* self);

void   hdr_set_release_coeff(HDR* self, double coeff /*0 to 1*/);
double hdr_get_release_coeff(HDR* self);

#if defined(__cplusplus)
}
#endif   //(__cplusplus)

#endif   // __HDR__
