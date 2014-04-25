#include "mel_functions.h"

float* powspec1(float* inputBuffer,float winpts,float steppts,int sr)
{
	float* WINDOW;
	int i;
	float* z_pad_window;

	/*inParam->testBuffer[0] = inParam->windowSize;
	inParam->testBuffer[1] = inParam->overlap;
	inParam->testBuffer[2] = inParam->stepSize;*/
	WINDOW = hanning(winpts);

	z_pad_window = (float*)calloc(winpts,sizeof(float));

	for(i = 0;i<winpts;i++)
	{
		z_pad_window[i] = inputBuffer[i]*WINDOW[i];
	}

	return z_pad_window;
}

float* spec2cep(float* spec)
{
	int nrow,i,j,ncep;
	nrow = 40;
	ncep = 13;
	float** dctm = (float**)malloc(ncep*sizeof(float*));

	float temp1,temp2;
	float* cep = (float*)calloc(ncep,sizeof(float));

	for(i = 0;i < ncep;i++)
	{
		dctm[i] = (float*)malloc(nrow*sizeof(float));
	}

	temp2 = (float)(sqrt((float)(2)/(float)(nrow)));

	for(i = 0;i < ncep;i++)
	{
		for(j = 0;j < nrow;j++)
		{
			temp1 = 1 + j*2;
			dctm[i][j] = cos(i*temp1/(2*nrow)*PI)*temp2;
		}
	}

	temp2 = (float)(sqrt(2));
	for(j = 0;j < nrow;j++)
	{
		dctm[0][j] = dctm[0][j]/temp2;
	}

	for(i = 0;i < ncep;i++)
	{
		cep[i] = 0;
		for(j = 0;j < nrow;j++)
		{
			temp2 = logf(spec[j]);
			cep[i] += dctm[i][j] * temp2;
		}
	}

	free(dctm);
	return(cep);
}

float* lifter(float* x,int lift)
{
	int ncep;
	int i;

	ncep = 13;
	lift = 0.6;

	float* liftwts = (float*)calloc(ncep,sizeof(float));

	for(i = 0; i < ncep; i++)
	{
		if(i == 0)
		{
			liftwts[i] = 1;
		}
		else
		{
			liftwts[i] = powf(i,lift);
		}
	}

	for(i = 0;i < ncep;i++)
	{
		liftwts[i] = liftwts[i]*x[i];
	}

	return liftwts;
}

float hz2mel(int f)
{
	int f_0;
	float f_sp;
	int brkfrq;
	float brkpt;
	float logstep;
	float z;

	f_0 = 0;
	f_sp = 200.0/3.0;
	brkfrq = 1000;
	brkpt = ((float)(brkfrq - f_0))/f_sp;
	logstep = exp(logf(6.4)/27);

	if(f<brkfrq)
	{
		z = (float)(f - f_0)/f_sp;
	}
	else
	{
		z = brkpt+(logf(f/brkfrq))/(logf(logstep));
	}

	return z;
}

float* mel2hz(int nfilts,float minmel,float maxmel)
{
	int f_0;
	float f_sp;
	int brkfrq;
	float brkpt;
	float logstep;
	int i;
	int num_nfilts;
	float* f = (float*)calloc(nfilts+2,sizeof(float));

	num_nfilts = nfilts + 2;
	f_0 = 0;
	f_sp = 200.0/3.0;
	brkfrq = 1000;
	brkpt = ((float)(brkfrq - f_0))/f_sp;
	logstep = exp(logf(6.4)/27);

	for(i = 0;i < num_nfilts;i ++)
	{
		f[i] = minmel + i * (maxmel - minmel)/(nfilts + 1);
		if(f[i] < brkpt)
		{
			f[i] = f_0 + f_sp * f[i];
		}
		else
		{
			f[i] = brkfrq * exp(logf(logstep)*(f[i] - brkpt));
		}
	}

	return f;
}

float** fft2melmx(void)
{
	float minmel;
	float maxmel;
	int minfrq;
	int maxfrq;
	float* binfrqs;
	int nfilts;
	int i,j;
	float fs[3];
	int width;
	float loslope;
	float hislope;
	int nfft;
	int binfrqs_size;
	float temp;
	int sr;
	float fftfrqs;
	int nfreqs;


	minfrq = 0;
	maxfrq = 4000;
	nfilts = 40;
	width = 1;
	nfft = 128;
	sr = 8000;
	nfreqs = 65;

	binfrqs_size = round(nfft/2)+1;

	float** wts ;
	wts= (float**)malloc(nfilts*sizeof(float*));

	for(i = 0;i < nfilts;i++)
	{
		wts[i] = (float*)malloc(nfft*sizeof(float));
	}

	minmel = hz2mel(minfrq);
	maxmel = hz2mel(maxfrq);
	binfrqs = mel2hz(nfilts,minmel,maxmel);

	for(i = 0;i < nfilts;i++)
	{
		for(j = 0;j < 3;j++)
		{
			fs[j] = binfrqs[i+j];
		}
		for(j = 0;j < 3;j++)
		{
			fs[j] = fs[1] + width * (fs[j] - fs[1]);
		}
		for(j = 0;j < binfrqs_size;j++)
		{
			fftfrqs = (float)(j*sr)/(float)(nfft);
			loslope = (fftfrqs - fs[0])/(fs[1] - fs[0]);
			hislope = (fs[2] - fftfrqs)/(fs[2] - fs[1]);
			temp = min(loslope,hislope);
			wts[i][j] = max(0,temp);
		}
	}

	// -- Matrix Multiply -- Diagnal * Normal Matrix
	for(i = 0;i < nfilts;i++)
	{
		for(j = 0;j < nfft;j++)
		{
			temp = 2/(binfrqs[i+2] - binfrqs[i]);
			wts[i][j] = wts[i][j] * temp;
		}
	}

	for(i = 0;i < nfilts;i++)
	{
		for(j = nfreqs;j < nfft;j++)
		{
			wts[i][j] = 0;
		}
	}

	return wts;
}

float* audspec(float* pspectrum, int nfreqs, int nfilts)
{
	float** wts;
	float binfrqs_size;
	int nfft;
	int i,j;
	float* aspectrum = (float*)calloc(nfilts,sizeof(float));

	nfft = (nfreqs - 1)*2;
	binfrqs_size = round(nfft/2)+1;

	wts = fft2melmx();

	for(i = 0;i < nfilts;i++)
	{
		aspectrum[i] = 0;
		for(j = 0;j < binfrqs_size;j++)
		{
			aspectrum[i] += wts[i][j] * pspectrum[j];
		}
	}

	return aspectrum;
}

float max(float a,float b)
{
	float max_value;

	if(a > b)
	{
		max_value = a;
	}
	else
	{
		max_value = b;
	}

	return max_value;
}

float min(float a,float b)
{
	float min_value;

	if(a > b)
	{
		min_value = b;
	}
	else
	{
		min_value = a;
	}

	return min_value;
}

float* hanning(int winpts)
{
	float* hanning_window = (float *)calloc(winpts,sizeof(float));
	int i;
	int bound;

	bound = winpts/2;

	for(i = 0;i<bound;i++)
	{
		hanning_window[i] = 0.5*(1-cosf(2*PI*(i+1)/(float)(winpts)));
	}

	for(i = bound;i<winpts;i++)
	{
		hanning_window[i] = hanning_window[winpts -1 - i];
	}

	return hanning_window;
}

float* melfcc(float* inputBuffer)
{
	int i;
	float winpts;
	float steppts;
	int sr;
	float wintime,hoptime;
	float* z_pad_inputBuffer;
	float NFFT;
	int select;
	float magnitude;
	float* fft_outputBuffer;
	float ds_real;
	float ds_imag;
	float sum_fft;
	int nfreqs,nfilts;
	float* aspectrum;
	float* cepstra;
	float lifterexp;
	float* cepstra_lift;

	wintime = 0.011;
	hoptime = 0.005;
	sr = 8000;

	winpts = roundf(wintime*sr);
	steppts = roundf(hoptime*sr);

	z_pad_inputBuffer = powspec1(inputBuffer,winpts,steppts,sr);

	NFFT = powf(2,(ceil(log(winpts)/log(2))));

	Transform* inTransform = newTransform(1, NFFT);

	inTransform ->doTransform(inTransform, z_pad_inputBuffer);

	select = (int)(round(NFFT/2));

	fft_outputBuffer = (float *)calloc(select+1,sizeof(float));

	sum_fft = 0;

	for(i = 0;i<=select;i++)
	{
		ds_real = inTransform->real[i];
		ds_imag = inTransform->imaginary[i];
		magnitude = ds_real * ds_real + ds_imag * ds_imag;
		fft_outputBuffer[i] = magnitude;
		sum_fft += magnitude;
	}

	free(z_pad_inputBuffer);

	// -- end of PowSpec

	nfreqs = select + 1;
	nfilts = 40;

	aspectrum = audspec(fft_outputBuffer,nfreqs,nfilts);

	cepstra = spec2cep(aspectrum);

	lifterexp = 0.600;

	cepstra_lift = lifter(cepstra,lifterexp);

	destroyTransform(&inTransform);

	return(cepstra_lift);
}



/*float hz2mel(int f)
{
	int f_0;
	float f_sp;
	int brkfrq;
	float brkpt;
	float logstep;
	float z;

	f_0 = 0;
	f_sp = 200.0/3.0;
	brkfrq = 1000;
	brkpt = ((float)(brkfrq - f_0))/f_sp;
	logstep = exp(logf(6.4)/27);
	if(f<brkfrq)
	{
		z = (float)(f - f_0)/f_sp;
	}
	else
	{
		z = brkpt+(logf(f)/brkfrq)/(logf(logstep));
	}

	return z;
}

float* mel2hz(int nfilts,float minmel,float maxmel)
{
	int f_0;
	float f_sp;
	int brkfrq;
	float brkpt;
	float logstep;
	int i;
	int num_nfilts;
	float* f = (float*)calloc(nfilts+2,sizeof(float));

	num_nfilts = nfilts + 2;
	f_0 = 0;
	f_sp = 200.0/3.0;
	brkfrq = 1000;
	brkpt = ((float)(brkfrq - f_0))/f_sp;
	logstep = exp(logf(6.4)/27);

	for(i = 0;i < num_nfilts;i ++)
	{
		f[i] = minmel + i * (maxmel - minmel)/(nfilts + 1);
		if(f[i] < brkpt)
		{
			f[i] = f_0 + f_sp * f[i];
		}
		else
		{
			f[i] = brkfrq * exp(logf(logstep)*(f[i] - brkpt));
		}
	}

	return f;
}*/
