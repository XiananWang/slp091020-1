#include <jni.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Timer.h"

typedef struct Variables {
    Timer* timer;
    float* inputBuffer;
    float* outputBuffer;
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

} Variables;

/*typedef struct TVPhamVAD {
	float D,Pl,Ph,Px;
	int Na,WinSize,Nf;
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
}TVPhamVAD;*/

/*void initialTVPhamVAD()
{
	Variables* inParam = (Variables*)malloc(sizeof(Variables));
	TVPhamVAD* VADRefer = (TVPhamVAD*) malloc(sizeof(TVPhamVAD));
	VADRefer->WinSize = inParam->windowSize;
	VADRefer->Nf = (inParam->frequency/VADRefer->WinSize)/2;
	VADRefer->epsqb = 0.001;
	VADRefer->alpha = 0.975;
	VADRefer->aDs = 0.65;
	VADRefer->Na = (VADRefer->WinSize)/2;
	VADRefer->Dbuf = (float*) calloc(VADRefer->Nf,sizeof(float));
	VADRefer->Dbufsort = (float*) calloc(VADRefer->Nf,sizeof(float));
	VADRefer->DbufInd = 0;
	VADRefer->PreDs = 0;
	VADRefer->PreTqb = 0;
	VADRefer->Firstrunflag = 1;
	VADRefer->VADflag = 0;
	VADRefer->NoVoiceCount = 0;
	VADRefer->VadDec = 0;
	VADRefer->interVadDec = 0;
	VADRefer->NoVoiceCount = 0;
}*/

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

void DecideVad(Variables *inParam)
{
	//TVPhamVAD* VADRefer = (TVPhamVAD*) malloc(sizeof(TVPhamVAD));
	//Variables* inParam = (Variables*) memoryPointer;
	float D,Pl,Ph,Px;
	float Dw,Dc,Ds;
	float Tqb;
	int i,qb;
	Px = 0;

	for(i=0;i<(inParam->windowSize);i++)
	{
		Px += inParam->inputBuffer[i]*inParam->inputBuffer[i];
	}
	Pl=0;
	Ph=0;
	for(i=0;i<(inParam->Na);i++)
	{
		Pl += (inParam->XL[i])*(inParam->XL[i]);
		//inParam->outputBuffer[i] = xL[i];
		Ph += (inParam->XH[i])*(inParam->XH[i]);
	}
	inParam->outputBuffer[0] = Pl;
	D = abs(Pl-Ph)/inParam->Na;
	Dw = D*(0.5 + (16/logf(2.0))*log(1+2*Px));
	Dw = -2*Dw;
	Dc = (1-exp(Dw))/(1+exp(Dw));
	Ds = Dc + inParam->PreDs*0.65;
	inParam->PreDs = Ds;
	inParam->Dbuf[inParam->DbufInd++] = Ds;
	if((inParam->DbufInd)>=(inParam->Nf))
		inParam->DbufInd = 0;

	//sort Dbuf
	sortDbuf();
	qb = 4;
	while((((inParam->Dbufsort[qb])-(inParam->Dbufsort[qb-4]))<(inParam->epsqb)) && (qb<(inParam->Nf)-1))
	{
		qb++;
	}
	Tqb = inParam->Dbufsort[qb];
	if(inParam->Firstrunflag)
		inParam->Firstrunflag = 0;
	else
		Tqb = (inParam->alpha)*(inParam->PreTqb) + (1-inParam->alpha)*(inParam->Tqb);

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
			inParam->VadDec = 0;
		else
		{
			inParam->VadDec = 1;
			inParam->NoVoiceCount++;
			if(inParam->NoVoiceCount >25)//guard time for 25 windows = 200ms
				inParam->VADflag = 1;
		}
	}
	inParam->outputBuffer[inParam->monitor++] = inParam->VADflag;
}

static void
compute(JNIEnv *env, jobject thiz,  jlong memoryPointer, jshortArray input)
{
    Variables* inParam = (Variables*) memoryPointer;
    //TVPhamVAD* VADRefer = (TVPhamVAD*) malloc(sizeof(TVPhamVAD));

    start(inParam->timer);

    short *_in = (*env)->GetShortArrayElements(env, input, NULL);

    int i,j,overlap, stepsize,order;
    order = 12;
    stepsize = inParam->stepSize;
    	float XL_temp = 0;
    	float XH_temp = 0;
    	overlap = inParam->overlap;
    	stepsize = inParam->stepSize;

    	for(i=0; i<overlap; i++)
    	{
    		inParam->inputBuffer[i] = inParam->inputBuffer[stepsize + i];
    	}

    	for (i=0; i<stepsize; i++)
    	{
    		inParam->inputBuffer[overlap + i] = _in[i];
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
    			XL_temp += inParam->inputBuffer[overlap + i - j]*LD[j];
    			XH_temp += inParam->inputBuffer[overlap + i - j]*HD[j];
    		}
    		if(i%2 == 0)
    		{
    			inParam->XL[i/2] = XL_temp;
    			inParam->XH[i/2] = XH_temp;
    			//inParam->outputBuffer[i/2] = inParam->XL[i/2];
    		}
    	}

    	/*for(i = 0;i<stepsize;i++)
    	{
    		inParam->outputBuffer[i] = inParam->XL[i];
    	}*/

    	DecideVad(inParam);
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
    inParam->timer = newTimer();
    inParam->frequency = freq;
    inParam->stepSize = stepsize;
    inParam->windowSize = windowsize;
    inParam->overlap = windowsize-stepsize;
    inParam->inputBuffer = (float*) calloc(windowsize,sizeof(float));
    //inParam->outputBuffer = (float*) malloc(stepsize*sizeof(float));
    inParam->outputBuffer = (float*) calloc(stepsize,sizeof(float));
    //-----------------------TVPhamVAD-----------------------------------------
	inParam->Nf = (inParam->frequency/inParam->windowSize)/2;
	inParam->epsqb = 0.001;
	inParam->alpha = 0.975;
	inParam->aDs = 0.65;
	inParam->Na = (inParam->windowSize)/2;
	inParam->Dbuf = (float*) calloc(inParam->Nf,sizeof(float));
	inParam->Dbufsort = (float*) calloc(inParam->Nf,sizeof(float));
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
    return (jlong)inParam;
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
        }

    } else {                //Case 2 - Processed output signal
                            //This should be your synthesized output.
        int i;
        for(i=0;i<inParam->stepSize;i++)
        {
            _output[i] = inParam->outputBuffer[i];
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

        debugOutput = (*env)->NewFloatArray(env, inParam->stepSize);
        float *_debugOutput = (*env)->GetFloatArrayElements(env, debugOutput, NULL);

        int i;
        for (i=0; i<inParam->stepSize;i++)
        {
            _debugOutput[i] = inParam->outputBuffer[i];
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
