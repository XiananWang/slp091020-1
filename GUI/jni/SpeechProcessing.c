#include <jni.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Timer.h"
#include "Transforms.h"
#include "mel_functions.h"
#include "logMMSE_functions.h"

#define PI 3.14159265

typedef struct Melfcc {
    float wintime;
    float hoptime;
    int numcep;
    float lifterexp;
    int sumpower;
    float preemph;
    int dither;
    int minfreq;
    int maxfreq;
    int nbands;
    int bwidth;
    int dcttype;
    int usecmp;
    int modelorder;
    int boraden;
    int useenergy;
} Melfcc;

typedef struct Variables {
    Timer* timer;
    float* inputBuffer;
    float* outputBuffer;
    float* testBuffer;
    int frequency;
    int stepSize;
    int windowSize;
    int overlap;
    //---------------TVPhamVAD--------------------
    float D,Pl,Ph,Px;
	int Na,Nf;
	float Dw,Dc,Ds,PreDs;
	float *Dbuf,*Dbufsort;
	float Tqb,epsqb,alpha,aDs;
	float PreTqb;
	int qb;
	int Firstrunflag;
	int DbufInd;
	int interVadDec,NoVoiceCount;
	int VADflag;
	int VadDec;
	//----------------VAD mid variant--------------
	float* XL;
	float* XH;
	//-----------------Xianan Define Monitor-------
	int monitor;
	//-----------------Melfcc----------------------
	Melfcc* inMelfcc;
} Variables;

/*typedef struct PowerSpec{
	float mel_error;
	float* fft_outputBuffer;
	float size;
}PowerSpec;*/

static const float LD[] =
{
	0.015404109327027373,
	0.0034907120842174702,
	-0.11799011114819057,
	-0.048311742585633,
	0.4910559419267466,
	0.787641141030194,
	0.3379294217276218,
	-0.07263752278646252,
	-0.021060292512300564,
	0.04472490177066578,
	0.0017677118642428036,
	-0.007800708325034148
};

static const float HD[] =
{
	0.007800708325034148,
	0.0017677118642428036,
	-0.04472490177066578,
	-0.021060292512300564,
	0.07263752278646252,
	0.3379294217276218,
	-0.787641141030194,
	0.4910559419267466,
	0.048311742585633,
	-0.11799011114819057,
	-0.0034907120842174702,
	0.015404109327027373
};

void DecideVad(Variables *inParam);
//void melfcc(Variables *inParam);
//PowerSpec* powspec(Variables *inParam,float winpts,float steppts);
//void audspec(void);
//float* hanning(int winpts,Variables *inParam);
//void melfcc_Initial(Melfcc* inMelfcc, int freq);
//PowerSpec* newPowerSpec(int size);

/*void melfcc(Variables *inParam)
{
	int i,j,index;
	int stepsize,order,overlap;
	float pre_temp = 0;
	float pre_coeff[2] = {};
	float winpts;
	float steppts;
	int sr;
	int zero_num;
	float wintime,hoptime;
	float* pre_emp_buffer = (float *)calloc(inParam->windowSize,sizeof(float));

	stepsize = inParam->stepSize;
	overlap = inParam->overlap;
	order = 2;
	pre_coeff[0] = 1;
	pre_coeff[1] = 0;

	wintime = inParam->inMelfcc->wintime;
	hoptime = inParam->inMelfcc->hoptime;
	sr = inParam->frequency;

	winpts = roundf(wintime*sr);
	steppts = roundf(hoptime*sr);
	zero_num = winpts - steppts;

	for(i=0;i<inParam->windowSize;i++)
	{
		pre_temp = 0;
		for(j=0;j<order;j++)
		{
			index = overlap + i - j;
			if(index >= 0)
			{
				if(index >= zero_num)
				{
					pre_temp += inParam->inputBuffer[index-zero_num]*pre_coeff[j];
				}
			}
		}
		pre_emp_buffer[i] = pre_temp;
	}

	powspec(inParam,winpts,steppts);
}*/

/*PowerSpec* powspec(Variables *inParam,float winpts,float steppts)
{
	float steptime;
	float wintime;
	int sr;
	int near_temp;
	float NFFT;
	float* WINDOW;
	int NOVERLAP,SAMPRATE;
	int i,j;
	float y_aft_wind;
	int z_pad_size;
	int select;

	inParam->testBuffer[0] = inParam->windowSize;
	inParam->testBuffer[1] = inParam->overlap;
	inParam->testBuffer[2] = inParam->stepSize;

	NFFT = pow(2,(ceil(log(winpts)/log(2))));
	Transform* inTransform = newTransform(1, NFFT);
	WINDOW = hanning(winpts,inParam);
	NOVERLAP = winpts - steppts;
	SAMPRATE = sr;

	float* z_pad_inputBuffer = (float *)calloc(inParam->windowSize,sizeof(float));

	for(i = 0;i<inParam->windowSize;i++)
	{
		z_pad_inputBuffer[i] = inParam->inputBuffer[i]*WINDOW[i];
	}

	inTransform = newTransform(1, NFFT);
	inTransform ->doTransform(inTransform, z_pad_inputBuffer);

	select = NFFT/2;

	//float* ds_real_Buffer = (float *)calloc(select+1,sizeof(float));
	//float* ds_imag_Buffer = (float *)calloc(select+1,sizeof(float));
	float* fft_outputBuffer = (float *)calloc(select+1,sizeof(float));
	PowerSpec* inPowerSpec = newPowerSpec(select+1);
	float ds_real;
	float ds_imag;
	float mel_error;

	mel_error = 0;

	for(i = 0;i<=select;i++)
	{
		ds_real = inTransform->real[i];
		ds_imag = inTransform->imaginary[i];
		inPowerSpec->fft_outputBuffer[i] = ds_real * ds_real + ds_imag * ds_imag;
		inPowerSpec->mel_error += fft_outputBuffer[i];
	}

	inPowerSpec->size = select+1;

	destroyTransform(&inTransform);

	return inPowerSpec;
}*/

/*void audspec(Variables* inParam,PowerSpec* inPowerSpec)
{
	int nfreqs,nfft;
	int nfilts,sr,width;
	float minfreq,maxfreq;
	int htkmel = 0;
	float minmel,maxmel,binfrqs;

	nfreqs = inPowerSpec->size;
	nfft = (nfreqs-1)*2;
	nfilts = inParam->inMelfcc->nbands;
	sr = inParam->frequency;
	width = inParam->inMelfcc->bwidth;
	minfreq = inParam->inMelfcc->minfreq;
	maxfreq = inParam->inMelfcc->maxfreq;

	minmel = hz2mel(minfreq);
	maxmel = hz2mel(maxfreq);

	binfrqs = mel2hz(minmel+[]);
	binfrqs = mel2hz(nfilts,)
}*/

/*PowerSpec* newPowerSpec(int size)
{
	PowerSpec* inPowerSpec = (PowerSpec*)malloc(sizeof(PowerSpec));
	inPowerSpec->fft_outputBuffer = (float *)calloc(size,sizeof(float));
	inPowerSpec->mel_error = 0;
	inPowerSpec->size = 0;
	return inPowerSpec;
}

float* hanning(int winpts,Variables *inParam)
{
	float* hanning_window = (float *)calloc(winpts,sizeof(float));
	int i;

	for(i = 1;i<=winpts;i++)
	{
		hanning_window[i-1] = 0.5*(1-cosf(2*PI*i/(float)(winpts)));
	}

	return hanning_window;
}*/

void DecideVad(Variables *inParam)
{
	float D,Px;
	float Pl,Ph;
	float D_test;
	float Dw,Dc,Ds;
	float Tqb;
	int qb;
	float temp;
	int i,j;
	int Na,winSize,Nf;

	winSize = inParam->stepSize;
	Na = winSize/2;
	Nf = inParam->frequency/winSize/2;

	Px = 0;

	//----------
	int counter = 0;
	int marker = 0;
	//----------

	for(i=0;i<winSize;i++)
	{
		Px += (inParam->inputBuffer[inParam->overlap+i])*(inParam->inputBuffer[inParam->overlap+i]);
	}

	Pl=0;
	Ph=0;

	for(i=0;i<Na;i++)
	{
		Pl += (inParam->XL[i])*(inParam->XL[i]);
		Ph += (inParam->XH[i])*(inParam->XH[i]);
	}

	D = (Pl-Ph)/(float)(Na);
	//D = (inParam->Pl - inParam->Ph)/(1.0*inParam->Na);
	if(D<0) D = -1 * D; // marker
	//inParam->outputBuffer[4] = D;
	Dw = D*(0.5 + (16/log(2.0))*log(1+2*Px)); // marker
	//inParam->outputBuffer[5] = Dw;
	Dw = -2*Dw;

	Dc = (1-exp(Dw))/(1+exp(Dw));
	Ds = Dc + inParam->PreDs*0.65;
	inParam->PreDs = Ds;
	inParam->Dbuf[inParam->DbufInd++] = Ds;

	if((inParam->DbufInd)>=Nf)
	{
		inParam->DbufInd = 0;
	}

	for(i=0;i<Nf;i++)
	{
		inParam->Dbufsort[i] = inParam->Dbuf[i];
	}
	for(i=Nf-1;i>0;i--) //marker
	{
		for(j=0;j<i;j++)
		{
			if((inParam->Dbufsort[j])>(inParam->Dbufsort[j+1]))
			{
				temp = inParam->Dbufsort[j];
				inParam->Dbufsort[j] = inParam->Dbufsort[j+1];
				inParam->Dbufsort[j+1] = temp;
			}
		}
	}
	qb = 4;

	while((((inParam->Dbufsort[qb])-(inParam->Dbufsort[qb-4]))<(inParam->epsqb)) && (qb<Nf-1))
	{
		qb++;
	}

	Tqb = inParam->Dbufsort[qb];
	if(inParam->Firstrunflag)
		inParam->Firstrunflag = 0;
	else
		Tqb = (inParam->alpha)*(inParam->PreTqb) + (1-inParam->alpha)*(Tqb);

	inParam->PreTqb = Tqb;

	if(Ds>Tqb)
		inParam->interVadDec = 1;
	else
		inParam->interVadDec = 0;

	if((inParam->interVadDec)!=0)
	{
		inParam->NoVoiceCount = 0;
		inParam->VADflag = 0;
		inParam->VadDec = 1;
	}
	else
	{
		if(inParam->VADflag)
		{
			inParam->VadDec = 0;
		}
		else
		{
			inParam->VadDec = 1;
			inParam->NoVoiceCount++;
			if(inParam->NoVoiceCount >25)//guard time for 25 windows = 200ms
				inParam->VADflag = 1;
		}
	}
	inParam->monitor++;
}

static void
compute(JNIEnv *env, jobject thiz,  jlong memoryPointer, jshortArray input)
{
    Variables* inParam = (Variables*) memoryPointer;
    //TVPhamVAD* VADRefer = (TVPhamVAD*) malloc(sizeof(TVPhamVAD));

    start(inParam->timer);

    short *_in = (*env)->GetShortArrayElements(env, input, NULL);

    int i,j,overlap, stepsize,order,index;
    order = 12;
    stepsize = inParam->stepSize;
    float* cepstra;

    	float XL_temp = 0;
    	float XH_temp = 0;

    	//---------------------------
    	int counter = 0;
    	//----------------------------

    	overlap = inParam->overlap;
    	stepsize = inParam->stepSize;

    	for(i=0; i<overlap; i++)
    	{
    		inParam->inputBuffer[i] = inParam->inputBuffer[stepsize + i];
    	}

    	for (i=0; i<stepsize; i++)
    	{
    		inParam->inputBuffer[overlap + i] = _in[i]/32768.0;
    		//inParam->outputBuffer[i] = inParam->inputBuffer[overlap + i];
    	}

    	(*env)->ReleaseShortArrayElements(env, input, _in, 0);

    	//Processing logic here. For now this placeholder will copy the input
    	//to the output array for the processing function and scale it by 1/2.
    	//For switching between suppressed and original output a GUI option
    	//can be used to select the output channel.

    	for(i=0;i<stepsize;i++)
    	{
    		XL_temp = 0;
    		XH_temp = 0;
    		for(j=0;j<order;j++)
    		{
    			index = overlap + i - j;
    			if(index>=0)
    			{
    				XL_temp += inParam->inputBuffer[overlap + i - j]*LD[j];
    				XH_temp += inParam->inputBuffer[overlap + i - j]*HD[j];
    			}
    		}

    		if(i%2 == 0)
    		{
    			inParam->XL[i/2] = XL_temp;
    			inParam->XH[i/2] = XH_temp;
    			//inParam->outputBuffer[i/2] = inParam->XH[i/2];
    		}
    	}

    	DecideVad(inParam);
    	cepstra = melfcc(inParam->inputBuffer);
    	//melfcc(inParam);
    /*for(i=0;i<stepsize;i++)
    {
        inParam->outputBuffer[i] = inParam->inputBuffer[overlap+i]*0.5f;
    }*/

    	/*for(i=0;i<stepsize;i++)
		{
			inParam->outputBuffer[i] = XL[i];
		}*/


    stop(inParam->timer);
}

static jlong
initialize(JNIEnv* env, jobject thiz, jint freq, jint stepsize, jint windowsize)
{
    Variables* inParam = (Variables*) malloc(sizeof(Variables));
    inParam->inMelfcc = (Melfcc*) malloc(sizeof(Melfcc));
    inParam->timer = newTimer();
    inParam->frequency = freq;
    inParam->stepSize = stepsize;
    inParam->windowSize = windowsize;
    inParam->overlap = windowsize-stepsize;
    //inParam->overlap = windowsize-stepsize;
    inParam->inputBuffer = (float*) calloc(windowsize,sizeof(float));
    //inParam->inputBuffer = (float*) calloc(windowsize,sizeof(float));
    //inParam->outputBuffer = (float*) malloc(stepsize*sizeof(float));
    inParam->outputBuffer = (float*) calloc(stepsize,sizeof(float));
    inParam->testBuffer = (float*) calloc(windowsize,sizeof(float));
    //-----------------------TVPhamVAD-----------------------------------------
    int winSize = inParam->stepSize;
	int Na = winSize/2;
	int Nf = inParam->frequency/winSize/2;

	inParam->Nf = (inParam->frequency/inParam->windowSize)/2;
	inParam->epsqb = 0.001;
	inParam->alpha = 0.975;
	inParam->aDs = 0.65;
	inParam->Na = (inParam->windowSize)/2;
	inParam->Dbuf = (float*) calloc(Nf,sizeof(float));
	inParam->Dbufsort = (float*) calloc(Nf,sizeof(float));
	inParam->DbufInd = 0;
	inParam->PreDs = 0;
	inParam->PreTqb = 0;
	inParam->Firstrunflag = 1;
	inParam->VADflag = 0;
	inParam->NoVoiceCount = 0;
	inParam->VadDec = 0;
	inParam->interVadDec = 0;
	//-----------------------VAD------------------------------------------------
	inParam->XL = (float*)malloc(sizeof(float)*stepsize);
	inParam->XH = (float*)malloc(sizeof(float)*stepsize);
    //initialTVPhamVAD(inParam);
	//-----------------------Monitor--------------------------------------------
	inParam->monitor = 0;
	//--------------------------------------------------------------------------
	melfcc_Initial(inParam->inMelfcc,inParam->frequency);

    return (jlong)inParam;
}

void melfcc_Initial(Melfcc* inMelfcc,int freq)
{
	inMelfcc->wintime = 0.011;
	inMelfcc->hoptime = 0.005;
	inMelfcc->numcep = 13;
	inMelfcc->lifterexp = 0.6;
	inMelfcc->sumpower = 1;
	inMelfcc->preemph = 0;
	inMelfcc->dither = 0;
	inMelfcc->minfreq = 0;
	inMelfcc->maxfreq = freq/2;
	inMelfcc->nbands = 40;
	inMelfcc->bwidth = 1.0;
	inMelfcc->dcttype = 2;
	inMelfcc->usecmp = 0;
	inMelfcc->modelorder = 0;
	inMelfcc->boraden = 0;
	inMelfcc->useenergy = 0;
}

static void
finish(JNIEnv* env, jobject thiz, jlong memoryPointer)
{
    Variables* inParam = (Variables*) memoryPointer;
    //cleanup memory
    if(inParam != NULL){
        tellTime(inParam->timer);
        destroy(&(inParam->timer));
        if(inParam->inputBuffer != NULL){
            free(inParam->inputBuffer);
            inParam->inputBuffer = NULL;
        }
        if(inParam->outputBuffer != NULL){
            free(inParam->outputBuffer);
            inParam->outputBuffer = NULL;
        }
        free(inParam);
        inParam = NULL;
    }
}

static jfloat
getTime(JNIEnv* env, jobject thiz, jlong memoryPointer)
{
    Variables* inParam = (Variables*) memoryPointer;
    return getMS(inParam->timer);
}

static jfloatArray
getOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelect)
{
    Variables* inParam = (Variables*) memoryPointer;

    jshortArray output = (*env)->NewShortArray(env, inParam->stepSize);
    short *_output = (*env)->GetShortArrayElements(env, output, NULL);

    if(outputSelect == 0) { //Case 1 - Original input signal
        int i;
        for(i=0;i<inParam->stepSize;i++)
        {
            _output[i] = inParam->inputBuffer[inParam->overlap+i];
        	//_output[i] = 6;
        }

    } else {                //Case 2 - Processed output signal
                            //This should be your synthesized output.
        int i;
        for(i=0;i<inParam->stepSize;i++)
        {
            _output[i] = inParam->outputBuffer[i];
        	//_output[i] = 7;
        }
    }

    (*env)->ReleaseShortArrayElements(env, output, _output, 0);
    return output;
}

static jfloatArray
getDebug(JNIEnv* env, jobject thiz, jlong memoryPointer, jint debugSelect)
{
    Variables* inParam = (Variables*) memoryPointer;

    jfloatArray debugOutput = NULL;

    if(debugSelect == 0) {          //Test Case 1 - inputBuffer contents

        debugOutput = (*env)->NewFloatArray(env, inParam->windowSize);
        float *_debugOutput = (*env)->GetFloatArrayElements(env, debugOutput, NULL);

        int i;
        for (i=0; i<inParam->windowSize;i++)
        {
            _debugOutput[i] = inParam->inputBuffer[i];
        }

        (*env)->ReleaseFloatArrayElements(env, debugOutput, _debugOutput, 0);
    } else if (debugSelect == 1) {  //Test Case 2 - outputBuffer contents

        /*debugOutput = (*env)->NewFloatArray(env, inParam->stepSize);
        float *_debugOutput = (*env)->GetFloatArrayElements(env, debugOutput, NULL);

        int i;
        for (i=0; i<inParam->stepSize;i++)
        {
            _debugOutput[i] = inParam->outputBuffer[i];
        }*/

    	debugOutput = (*env)->NewFloatArray(env, inParam->windowSize);
    	float *_debugOutput = (*env)->GetFloatArrayElements(env, debugOutput, NULL);
    	int i;
        for(i = 0;i<inParam->windowSize;i++)
        {
        	_debugOutput[i] = inParam->testBuffer[i];
        }

        (*env)->ReleaseFloatArrayElements(env, debugOutput, _debugOutput, 0);
    }

    //Add additional cases to output other data that may be needed.

    return debugOutput;
}

////////////////////////////////////////////////////////////////////////////////////////////
// JNI Setup - Functions and OnLoad
////////////////////////////////////////////////////////////////////////////////////////////

static JNINativeMethod nativeMethods[] =
    {//     Name                            Signature           Pointer
            {"compute",                     "(J[S)V",           (void *)&compute                },
            {"initialize",                  "(III)J",           (void *)&initialize             },
            {"finish",                      "(J)V",             (void *)&finish                 },
            {"getTime",                     "(J)F",             (void *)&getTime                },
            {"getOutput",                   "(JI)[S",           (void *)&getOutput              },
            {"getDebug",                    "(JI)[F",           (void *)&getDebug               }
    };

jint
JNI_OnLoad(JavaVM* vm, void* reserved)
{
    JNIEnv* env;
    jint result;
    //get a hook to the environment
    result = (*vm)->GetEnv(vm, (void**) &env, JNI_VERSION_1_6);
    if (result == JNI_OK) {
        //find the java class to hook the native methods to
        jclass filters = (*env)->FindClass(env, "com/dsp/speechpipeline/SpeechProcessing");
        if (filters != NULL) {
            result = (*env)->RegisterNatives(env, filters, nativeMethods, sizeof(nativeMethods)/sizeof(nativeMethods[0]));
            (*env)->DeleteLocalRef(env, filters);
            if(result == JNI_OK){
                return JNI_VERSION_1_6;
            } else {
                //something went wrong with the method registration
                return JNI_ERR;
            }
        } else {
            //class wasn't found
            return JNI_ERR;
        }
    } else {
        //could not get environment
        return JNI_ERR;
    }
}
