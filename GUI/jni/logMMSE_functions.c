#include "logMMSE_functions.h"

#define PI 3.14159265

void logMMSE(float* inputBuffer,int sr,int winduration,int hoptime)
{
	int len,len1,len2;
	float* WINDOW;
	int NFFT;
	int select;
	float* buffer_window;
	int i;

	winduration = 11;
	hoptime = 5;
	sr = 8000;
	// -- Initial --


	len = (int)(round(winduration*sr));
	len2 = (int)(round(hoptime*sr));
	len1 = len - len2;

	WINDOW = hanning(len);

	for(i = 0;i<len;i++)
	{
		buffer_window[i] = inputBuffer[i]*WINDOW[i];
	}

	NFFT = pow(2,ceil(logf(len)/logf(2)));
	Transform* inTransform1 = newTransform(1,NFFT);

	inTransform1->doTransform(inTransform1,inputBuffer);

	// -- Initialize Parameters --

	int k;
	float aa,eta,mu,c,qk,qkr,ksi_min;
	int RangePriormax,RangePriormin;
	int RangePostmin,RangePostmax;
	int Tkm;
	float alphap,alphad;
	float wmin;
	int NumWinPmin;
	float SafetyNetB,SafetyNetEta;

	k=1;
	aa=0.98;
	eta= 0.15;
	mu=0.98;
	c=sqrt(PI)/2;
	qk=0.3;
	qkr=(1-qk)/qk;
	ksi_min=powf(10,((float)(-25)/(float)(10)));
	RangePriormax = 40;
	RangePriormin = -19;
	RangePostmin = -30;
	RangePostmax = 40;
	Tkm = 4;
	alphap = 0.1;
	alphad = 0.85;

	//SafetyNetPmin = realmax.*ones(1,nFFT);
	wmin = 0.8*sr;
	NumWinPmin = floor(wmin/len1);
	//SafetyNetP = zeros(NumWinPmin,nFFT);
	SafetyNetB = 1.5;
	SafetyNetEta = 0.1;
}
