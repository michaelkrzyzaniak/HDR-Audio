//gcc *.c ../src/*.c

#include "../HDR.h"
#include "MKAiff.h"

#define WINDOW_SIZE 128

void hdr_onprocess_callback(void* SELF, dft_sample_t* real, int N);

/*--------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  if(argc < 4)
    {fprintf(stderr, "This command should be called with three arguments. 1) A high-gain wav or aiff file; 2) A low-gain wav or aiff file; 3) A filename where you would like the output to be saved as a wav file.\r\n"); exit(-1);}
  
  /* open audio files and check conformity */
  MKAiff* high_gain = aiffWithContentsOfFile(argv[1]);
  MKAiff* low_gain  = aiffWithContentsOfFile(argv[2]);
  
  if(high_gain == NULL)
    {perror("unable to open High Gain aiff"); exit(-1);}
  if(low_gain == NULL)
    {perror("unable to open Low Gain aiff"); exit(-1);}

  int num_chans_h = aiffNumChannels(high_gain);
  int num_chans_l = aiffNumChannels(low_gain);
  
  if(num_chans_h != 1)
    {fprintf(stderr, "%s has %i channels but must have exactly one\r\n", argv[1], num_chans_h); exit(-1);}
  if(num_chans_l != 1)
    {fprintf(stderr, "%s has %i channels but must have exactly one\r\n", argv[2], num_chans_l); exit(-1);}

  double sample_rate_h = aiffSampleRate(high_gain);
  double sample_rate_l = aiffSampleRate(low_gain);
  if(sample_rate_h != sample_rate_l)
    {fprintf(stderr, "%s has a sample rate of %f Hz and %s has %f Hz but they must match\r\n", argv[1], sample_rate_h, argv[2], sample_rate_l); exit(-1);}
  
  MKAiff* output_aiff    = aiffWithDurationInSamples(1, (unsigned long)sample_rate_h, aiffBitsPerSample(high_gain), aiffDurationInSamples(high_gain));
  if(output_aiff == NULL)
    {perror("Unable to create output wave file object."); exit(-1);}
    
  HDR* hdr = hdr_new(WINDOW_SIZE);
  if(hdr == NULL)
    {perror("Unable to create output wave file object."); exit(-1);}
  
  aiffRewindPlayheadToBeginning(high_gain);
  aiffRewindPlayheadToBeginning(low_gain);
  
  fprintf(stderr, "Processing %s (High Gain) and %s (Low Gain) ... \r\n", argv[1], argv[2]);
  
  for(;;)
    {
      /* Here audio will be fed in the the HDR object one sample at a time
      but it can be done in buffers of any size or even all at once, it will
      internally break the stream up into windows of the specified size and call
      the specified callback after processing each window.
      */
      
      int32_t high_gain_sample_i, low_gain_sample_i;
      double high_gain_sample, low_gain_sample;
      
      int num_samples_read_high = aiffReadIntegerSamplesAtPlayhead (high_gain, &high_gain_sample_i, 1, aiffYes);
      int num_samples_read_low  = aiffReadIntegerSamplesAtPlayhead (low_gain , &low_gain_sample_i , 1, aiffYes);
      
      high_gain_sample = high_gain_sample_i / (double)0x7FFFFFFF;
      low_gain_sample  = low_gain_sample_i  / (double)0x7FFFFFFF;
  
      if((num_samples_read_high==1) && (num_samples_read_low==1))
        hdr_process(hdr, &high_gain_sample, &low_gain_sample, 1, sample_rate_l, hdr_onprocess_callback, output_aiff);
      else
        break;
    }
  
  hdr_flush(hdr, hdr_onprocess_callback, output_aiff);

  aiffSaveWaveWithFilename(output_aiff, argv[3]);
  fprintf(stderr, "Done. Saved output wav file as %s\r\n", argv[3]);
  
  /* cleanup memory */
  aiffDestroy(high_gain);
  aiffDestroy(low_gain);
  aiffDestroy(output_aiff);
  hdr_destroy(hdr);
}

/*--------------------------------------------------------------------*/
void hdr_onprocess_callback(void* SELF, dft_sample_t* real, int N)
{
  MKAiff* output_aiff = SELF;
  aiffAppendFloatingPointSamples(output_aiff, real, N, aiffDoubleSampleType);
}
