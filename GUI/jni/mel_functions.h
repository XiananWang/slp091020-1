#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

float* powspec1(float* inputBuffer,float winpts,float steppts,int sr);
float* spec2cep(float* spec);
float* lifter(float* x,int lift);
float* mel2hz(int nfilts,float minmel,float maxmel);
float hz2mel(int f);
float** fft2melmx();
float* audspec(float* pspectrum, int nfreqs, int nfilts);
float max(float a,float b);
float min(float a,float b);
float* hanning(int winpts);
float* melfcc(float* inputBuffer);
