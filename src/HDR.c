#include "../HDR.h"

#include <stdlib.h> //calloc
#include <math.h>

/*--------------------------------------------------------------------*/
struct Opaque_HDR_Struct
{
  int           window_size;
  int           fft_N;
  int           fft_N_over_2;
  int           overlap;
  int           hop_size;
  int           return_size;
  int           sample_counter;
  int           input_index;
  dft_sample_t* window;
  dft_sample_t* running_input_h;
  dft_sample_t* running_input_l;
  dft_sample_t* running_output;
  dft_sample_t* real_h;
  dft_sample_t* imag_h;
  dft_sample_t* real_l;
  dft_sample_t* imag_l;
  
  dft_sample_t* filtered_crossover_gains;
  dft_sample_t  filtered_crossover_gain;

  double release_coeff;
  double crossover_min;
  double crossover_max;
};
  
/*--------------------------------------------------------------------*/
HDR* hdr_new(int window_size)
{
  HDR* self = calloc(1, sizeof(*self));
  if(self != NULL)
    {
      int i;
      self->window_size          = window_size;
      self->fft_N                = 2 * window_size;
      
      self->fft_N_over_2         = self->fft_N >> 1;
      self->overlap              = 2;
      self->hop_size             = window_size / self->overlap;

      self->return_size        = self->window_size;
      
      self->sample_counter       = 0;
      self->input_index          = 0;

      self->window               = calloc(self->window_size, sizeof(*(self->window)));
      self->running_input_h      = calloc(self->window_size, sizeof(*(self->running_input_h)));
      self->running_input_l      = calloc(self->window_size, sizeof(*(self->running_input_l)));
      
      self->running_output       = calloc(self->fft_N       , sizeof(*(self->running_input_h)));
      self->real_h               = calloc(self->fft_N       , sizeof(*(self->real_h)));
      self->imag_h               = calloc(self->fft_N       , sizeof(*(self->imag_h)));
      self->real_l               = calloc(self->fft_N       , sizeof(*(self->real_l)));
      self->imag_l               = calloc(self->fft_N       , sizeof(*(self->imag_l)));
      self->filtered_crossover_gains               = calloc(self->fft_N_over_2, sizeof(*(self->filtered_crossover_gains)));
      
      if(self->window           == NULL) return hdr_destroy(self);
      if(self->running_input_h  == NULL) return hdr_destroy(self);
      if(self->running_input_h  == NULL) return hdr_destroy(self);
      if(self->running_output   == NULL) return hdr_destroy(self);
      if(self->real_h           == NULL) return hdr_destroy(self);
      if(self->imag_h           == NULL) return hdr_destroy(self);
      if(self->real_l           == NULL) return hdr_destroy(self);
      if(self->imag_l           == NULL) return hdr_destroy(self);
      if(self->filtered_crossover_gains   == NULL) return hdr_destroy(self);

      dft_init_half_sine_window(self->window, self->window_size);
      
      self->crossover_min = HDR_DEFAULT_CROSSOVER_MIN;
      self->crossover_max = HDR_DEFAULT_CROSSOVER_MAX;
      self->release_coeff = HDR_DEFAULT_RELEASE_COEFF;
    }
    
  return self;
}

/*--------------------------------------------------------------------*/
HDR* hdr_destroy(HDR* self)
{
  if(self != NULL)
    {
      if(self->window != NULL)
        free(self->window);
      if(self->running_input_h != NULL)
        free(self->running_input_h);
      if(self->running_input_l != NULL)
        free(self->running_input_l);
      if(self->running_output != NULL)
        free(self->running_output);
      if(self->real_h != NULL)
        free(self->real_h);
      if(self->imag_h != NULL)
        free(self->imag_h);
      if(self->real_l != NULL)
        free(self->real_l);
      if(self->imag_l != NULL)
        free(self->imag_l);
      if(self->filtered_crossover_gains != NULL)
        free(self->filtered_crossover_gains);
      free(self);
    }
  return (HDR*) NULL;
}

/*--------------------------------------------------------------------*/
void hdr_flush(HDR* self, hdr_onprocess_t onprocess, void* onprocess_self)
{
  onprocess(onprocess_self, self->running_output+self->hop_size, self->window_size-self->hop_size);
}

/*--------------------------------------------------------------------*/
void   hdr_set_crossover_min(HDR* self, double coeff_dB)
{
  self->crossover_min = coeff_dB;
}

/*--------------------------------------------------------------------*/
double hdr_get_crossover_min(HDR* self)
{
  return self->crossover_min;
}

/*--------------------------------------------------------------------*/
void   hdr_set_crossover_max(HDR* self, double coeff_dB)
{
  self->crossover_max = coeff_dB;
}

/*--------------------------------------------------------------------*/
double hdr_get_crossover_max(HDR* self)
{
  return self->crossover_max;
}

/*--------------------------------------------------------------------*/
void   hdr_set_release_coeff(HDR* self, double coeff)
{
  coeff = (coeff < 0) ? 0 : (coeff > 1) ? 1 : coeff;
  self->release_coeff = coeff;
}

/*--------------------------------------------------------------------*/
double hdr_get_release_coeff(HDR* self)
{
  return self->release_coeff;
}

/*--------------------------------------------------------------------*/
double hdr_calculate_rms(HDR* self, dft_sample_t* signal, int n)
{
  if(n<=0) return 0;
  
  int i;
  double result = 0;
  for(i=0; i<n; i++)
    {
      result += *signal * *signal;
      ++signal;
    }
  result /= n;
  result = sqrt(result);
  
  return result;
}

/*--------------------------------------------------------------------*/
double hdr_calculate_peak_amp(HDR* self, dft_sample_t* signal, int n)
{
  int i;
  double temp, peak = 0;
  for(i=0; i<n; i++)
    {
      temp = fabs(*signal++);
      if(temp > peak)
        peak = temp;
    }
  
  return peak;
}

/*--------------------------------------------------------------------*/
void hdr_process(HDR* self, dft_sample_t* high_gain_input, dft_sample_t* low_gain_input, int len, double sample_rate, hdr_onprocess_t onprocess, void* onprocess_self)
{
  int i, j;
  double amplitude;

  for(i=0; i<len; i++)
    {
      self->running_input_h[self->input_index] = high_gain_input[i];
      self->running_input_l[self->input_index] = low_gain_input[i];

      ++self->input_index;
      self->input_index %= self->window_size;
    
      if(++self->sample_counter == self->hop_size)
        {
          self->sample_counter = 0;
        
          for(j=0; j<self->window_size; j++)
            {
              self->real_h[j] = self->running_input_h[(self->input_index+j) % self->window_size];
              self->real_l[j] = self->running_input_l[(self->input_index+j) % self->window_size];
            }
            
          for(j=self->window_size; j<self->fft_N; j++)
            self->real_h[j] = self->real_l[j] = 0;

          dft_apply_window(self->real_l, self->window, self->window_size);
          dft_apply_window(self->real_h, self->window, self->window_size);
          dft_2_real_forward_dfts(self->real_h, self->real_l, self->imag_h, self->imag_l, self->fft_N);

          if(self->crossover_max > self->crossover_min)
            {
              double magnitude, min_crossover, max_crossover, crossover_gain;
              for(j=0; j<self->fft_N_over_2; j++)
              {
                  magnitude = hypot(self->real_h[j], self->imag_h[j]);
                  magnitude = 20 * log10(magnitude);

                  //todo: could preecompute some of  this
                  crossover_gain = 44.723 + -20 * log10(dft_frequency_of_bin(j+1, sample_rate, self->fft_N));
                  min_crossover = crossover_gain + self->crossover_min;
                  max_crossover = crossover_gain + self->crossover_max;
                  crossover_gain = (magnitude - min_crossover) / (max_crossover - min_crossover);
                  if(crossover_gain > 1) crossover_gain = 1;
                  if(crossover_gain < 0) crossover_gain = 0;
                  
                  if(crossover_gain > self->filtered_crossover_gains[j])
                    self->filtered_crossover_gains[j] = crossover_gain;
                  else
                    self->filtered_crossover_gains[j] = (self->filtered_crossover_gains[j] * self->release_coeff) + (crossover_gain * (1-self->release_coeff));
                  
                  crossover_gain = self->filtered_crossover_gains[j];
                  
                  //fprintf(stderr, "%f\r\n", magnitude);
                  
                  /*
                  if(j==4) //689 Hz
                    fprintf(stderr, "%f\t", crossover_gain);
                  if(j==40) //6891 Hz
                    fprintf(stderr, "%f\r\n", crossover_gain);
                  */
                  
                  self->real_h[j] = (self->real_h[j] * (1-crossover_gain)) + (self->real_l[j] * crossover_gain);
                  self->imag_h[j] = (self->imag_h[j] * (1-crossover_gain)) + (self->imag_l[j] * crossover_gain);
                }
            }

          dft_real_inverse_dft(self->real_h, self->imag_h, self->fft_N);
          dft_apply_window(self->real_h, self->window, self->window_size);
              
          for(j=0; j<self->window_size-self->hop_size; j++)
            {
              self->real_h[j] += self->running_output[j+self->hop_size];
              self->running_output[j] = self->real_h[j];
            }
          for(; j<self->window_size; j++)
            self->running_output[j] = self->real_h[j];

          onprocess(onprocess_self, self->real_h, self->hop_size);
        }
    }
}
