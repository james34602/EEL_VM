#include <stdlib.h>
#define DENORMAL_BUFFER (128)
typedef struct
{
	unsigned int N, m, L, N2;
	float *subbandData, *freqLabel;
	float *channelMatrix, *h;
	float *allpass_delay_chain, *virtualHilbertTransformDelay;
	float *APC_delay_1, *APC_delay_2, *Xk2;
	unsigned int *Sk, *decimationCounter;
	float alpha, postGain;
	// Denormal buster
	float noiseBuffer[DENORMAL_BUFFER];
	unsigned int noiseLoop;
} WarpedPFB;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
size_t getMemSizeWarpedPFB(unsigned int N, unsigned int m)
{
	unsigned int L = 2 * m * N;
	return sizeof(WarpedPFB) + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)) + (N * sizeof(float)) +
		(L * sizeof(float)) + (L * sizeof(float)) + (2 * m * N * sizeof(float));
}
float add_denormal_prevention_white_noise(unsigned int *rand_state)
{
	*rand_state = *rand_state * 1234567UL + 890123UL;
	int mantissa = *rand_state & 0x007F0000; // Keep only most significant bits
	int flt_rnd = mantissa | 0x1E000000; // Set exponent
	return *((float*)(&flt_rnd));
}
#include <math.h>
float* allpass_char(double alpha, unsigned int L, unsigned int *CFiltLen);
extern void cos_fib_paraunitary1(unsigned int N, unsigned int m, unsigned int L, double df, double *h_opt);
extern void subsamplingCal(unsigned int M, double alpha, double *f_def, unsigned int *Sk);
void initWarpedPFB(WarpedPFB *pfb, double fs, unsigned int N, unsigned int m)
{
	char *memBuf = (char*)(pfb)+sizeof(WarpedPFB);
	unsigned int i, j;
	pfb->noiseLoop = 0;
	unsigned int seed = 1337;
	for (i = 0; i < DENORMAL_BUFFER; i++)
		pfb->noiseBuffer[i] = add_denormal_prevention_white_noise(&seed);
	pfb->N = N;
	pfb->N2 = pfb->N * 2;
	pfb->m = m;
	unsigned int L = 2 * m * N;
	pfb->L = L;
	double df = 1.0 / (1.5 * max(N, m));
	double *h_opt = (double*)malloc(L * sizeof(double));
	cos_fib_paraunitary1(N, m, L, df, h_opt);
	pfb->postGain = (float)(1.0 / (double)(2 * 2 * N));
	double alpha = -(0.1957 - 1.048* sqrt((2.0 / M_PI) * atan(0.07212 * (fs / 1000.0))));
	if (alpha > 0.99)
		alpha = 0.99;
	pfb->alpha = (float)alpha;
	double *f_def = (double*)malloc((N + 1) * sizeof(double));
	pfb->Sk = (unsigned int*)(memBuf);
	subsamplingCal(N, -alpha, f_def, pfb->Sk);
	pfb->freqLabel = (float*)(memBuf + (N * sizeof(unsigned int))); // fcentre
	for (i = 1; i < N - 1; i++)
		pfb->freqLabel[i] = (float)((f_def[i] + f_def[i + 1]) * 0.5 * fs);
	pfb->freqLabel[0] = 0.0f;
	pfb->freqLabel[N - 1] = fs * 0.5;
	free(f_def);
	//for (i = 0; i < N; i++)
	//	printf("%1.7f %d\n", freqLabel[i], Sk[i]);
	unsigned char maximumDownsampling = 1;
	if ((N == 2) && (alpha == 0.0) && maximumDownsampling) // Half band mode
		pfb->Sk[1] = pfb->Sk[0] = 2;
	pfb->virtualHilbertTransformDelay = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float))); // fcentre
	for (i = 0; i < N; i++)
	{
		if (pfb->freqLabel[i] > 0.0f)
			pfb->virtualHilbertTransformDelay[i] = (float)(((fs / (double)pfb->Sk[i]) / (double)pfb->freqLabel[i]) / 4.0);
		else
			pfb->virtualHilbertTransformDelay[i] = 0.0f;
	}
	float ovpRatio = 0.0f;
	for (i = 0; i < N; i++)
		ovpRatio += 1.0f / pfb->Sk[i];
	//printf("Total oversampling ratio = %1.3f\n", ovpRatio);
	/*FILE *fp1 = fopen("c.dat", "wb");
	fwrite(C, sizeof(float), CFiltLen, fp1);
	fclose(fp1);*/
	// Analysis preparation
	pfb->allpass_delay_chain = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)));
	memset(pfb->allpass_delay_chain, 0, L * sizeof(float));
	pfb->h = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)));
	for (i = 0; i < m; i++)
	{
		if (((i + 1) % 2) == 0)
		{
			for (j = 0; j < 2 * N; j++)
				pfb->h[i * 2 * N + j] = (float)(-h_opt[i * 2 * N + j]);
		}
		else
		{
			for (j = 0; j < 2 * N; j++)
				pfb->h[i * 2 * N + j] = (float)(h_opt[i * 2 * N + j]);
		}
	}
	free(h_opt);
	//for (i = 0; i < L; i++)
	//	printf("%d %1.7f\n", i + 1, h[i]);
	pfb->channelMatrix = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)));
	for (unsigned int k = 0; k < N; k++)
	{
		if ((k % 2) == 0)
		{
			for (unsigned int l = 0; l < 2 * N; l++)
				pfb->channelMatrix[k * (2 * N) + l] = (float)(2.0 * cos((2.0 * k + 1.0) * (M_PI / (2.0 * N)) * ((double)l - (L - 1.0) / 2.0) + (M_PI / 4.0)));
		}
		else
		{
			for (unsigned int l = 0; l < 2 * N; l++)
				pfb->channelMatrix[k * (2 * N) + l] = (float)(2.0 * cos((2.0 * k + 1.0) * (M_PI / (2.0 * N)) * ((double)l - (L - 1.0) / 2.0) - (M_PI / 4.0)));
		}
	}
	//printFloatMatrix2File("res.txt", c, N, 2 * N);
	pfb->decimationCounter = (unsigned int*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)));
	for (i = 0; i < N; i++)
		pfb->decimationCounter[i] = pfb->Sk[i];
	pfb->subbandData = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)));
	pfb->APC_delay_1 = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)) + (N * sizeof(float)));
	pfb->APC_delay_2 = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)) + (N * sizeof(float)) +
		(L * sizeof(float)));
	memset(pfb->APC_delay_1, 0, L * sizeof(float));
	memset(pfb->APC_delay_2, 0, L * sizeof(float));
	pfb->Xk2 = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)) + (N * sizeof(float)) +
		(L * sizeof(float)) + (L * sizeof(float)));
}
void assignPtrWarpedPFB(WarpedPFB *pfb, unsigned int N, unsigned int m)
{
	char *memBuf = (char*)(pfb)+sizeof(WarpedPFB);
	unsigned int L = 2 * m * N;
	pfb->allpass_delay_chain = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)));
	memset(pfb->allpass_delay_chain, 0, L * sizeof(float));
	pfb->subbandData = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)));
	pfb->APC_delay_1 = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)) + (N * sizeof(float)));
	pfb->APC_delay_2 = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)) + (N * sizeof(float)) +
		(L * sizeof(float)));
	memset(pfb->APC_delay_1, 0, L * sizeof(float));
	memset(pfb->APC_delay_2, 0, L * sizeof(float));
	pfb->Xk2 = (float*)(memBuf + (N * sizeof(unsigned int)) + ((N + 1) * sizeof(float)) + (N * sizeof(float)) + (L * sizeof(float)) +
		(L * sizeof(float)) + (N * 2 * N * sizeof(float)) + (N * sizeof(unsigned int)) + (N * sizeof(float)) +
		(L * sizeof(float)) + (L * sizeof(float)));
}
void changeWarpingFactorWarpedPFB(WarpedPFB *pfb, float fs, float pfb_log_grid_den)
{
	double alpha = (-(0.1957 - 1.048* sqrt((2.0 / M_PI) * atan(0.07212 * (fs / 1000.0))))) * pfb_log_grid_den;
	if (alpha < -0.99f)
		alpha = -0.99f;
	if (alpha > 0.99f)
		alpha = 0.99f;
	if (pfb->alpha != alpha)
	{
		pfb->alpha = alpha;
		double *f_def = (double*)malloc((pfb->N + 1) * sizeof(double));
		subsamplingCal(pfb->N, -pfb->alpha, f_def, pfb->Sk);
		for (int i = 1; i < pfb->N - 1; i++)
			pfb->freqLabel[i] = (float)((f_def[i] + f_def[i + 1]) * 0.5 * fs);
		pfb->freqLabel[0] = 0.0f;
		pfb->freqLabel[pfb->N - 1] = fs * 0.5;
		free(f_def);
		unsigned char maximumDownsampling = 1;
		if ((pfb->N == 2) && (pfb->alpha == 0.0) && maximumDownsampling) // Half band mode
			pfb->Sk[1] = pfb->Sk[0] = 2;
		for (int i = 0; i < pfb->N; i++)
		{
			if (pfb->freqLabel[i] > 0.0f)
				pfb->virtualHilbertTransformDelay[i] = (float)(((fs / (double)pfb->Sk[i]) / (double)pfb->freqLabel[i]) / 4.0);
			else
				pfb->virtualHilbertTransformDelay[i] = 0.0f;
		}
		float ovpRatio = 0.0f;
		for (int i = 0; i < pfb->N; i++)
			ovpRatio += 1.0f / pfb->Sk[i];
	}
}
void getSubbandDatWarpedPFB(WarpedPFB *pfb, float *subbands, float *curSk)
{
	for (unsigned int i = 0; i < pfb->N; i++)
	{
		subbands[i] = pfb->subbandData[i];
		curSk[i] = (float)pfb->decimationCounter[i];
	}
}
void getSubbandDatWarpedPFBStereo(WarpedPFB *pfb1, WarpedPFB *pfb2, float *subbands1, float *subbands2, float *curSk)
{
	for (unsigned int i = 0; i < pfb1->N; i++)
	{
		subbands1[i] = pfb1->subbandData[i];
		subbands2[i] = pfb2->subbandData[i];
		curSk[i] = (float)pfb1->decimationCounter[i];
	}
}
void writeSubbandDatWarpedPFB(WarpedPFB *pfb, float *subbands)
{
	for (unsigned int i = 0; i < pfb->N; i++)
		pfb->subbandData[i] = subbands[i];
}
void writeSubbandDatWarpedPFBStereo(WarpedPFB *pfb1, WarpedPFB *pfb2, float *subbands1, float *subbands2)
{
	for (int i = 0; i < pfb1->N; i++)
	{
		pfb1->subbandData[i] = subbands1[i];
		pfb2->subbandData[i] = subbands2[i];
	}
}
float* getPhaseCorrFilterWarpedPFB(WarpedPFB *pfb, float phCorrAlpha, unsigned int *CFiltLen)
{
	if (phCorrAlpha < -0.99f)
		phCorrAlpha = -0.99f;
	if (phCorrAlpha > 0.99f)
		phCorrAlpha = 0.99f;
	return allpass_char(pfb->alpha * phCorrAlpha, pfb->L, CFiltLen);
}
void analysisWarpedPFB(WarpedPFB *pfb, float x)
{
	unsigned int i, j;
	// Warping network
	// Complexity: O(n)
	float R1_a = x + pfb->noiseBuffer[pfb->noiseLoop];
	pfb->noiseLoop = (pfb->noiseLoop + 1) & (DENORMAL_BUFFER - 1);
	float R2_a = pfb->allpass_delay_chain[0];
	pfb->allpass_delay_chain[0] = R1_a;
	for (i = 1; i < pfb->L; i++)
	{
		float R3_a = R2_a;
		R2_a = pfb->allpass_delay_chain[i];
		R1_a = (R2_a - R1_a) * pfb->alpha + R3_a;
		pfb->allpass_delay_chain[i] = R1_a;
	}
	// Analysis, multiply by the coefficients of the prototype filter
	// Complexity: O(2 * m * M + 2 * M * sum(1 ./ Sk))
	memset(pfb->subbandData, 0, pfb->N * sizeof(float));
	for (i = 0; i < pfb->N2; i++)
	{
		float Xk = 0.0f;
		for (j = 0; j < pfb->m; j++)
			Xk += pfb->allpass_delay_chain[i + j * pfb->N2] * pfb->h[i + j * pfb->N2];
		//printf("%1.8f\n", Xk[i]);
		// Cosine modulation
		for (j = 0; j < pfb->N; j++)
		{
			if (pfb->decimationCounter[j] == pfb->Sk[j])
				pfb->subbandData[j] += pfb->channelMatrix[j * pfb->N2 + i] * Xk;
		}
	}
	for (j = 0; j < pfb->N; j++)
		pfb->subbandData[j] *= pfb->Sk[j];
}
float synthesisWarpedPFB(WarpedPFB *pfb)
{
	unsigned int i, j;
	// Synthesis
	// Complexity same as analysis
	for (i = 0; i < pfb->N2; i++)
	{
		float Curr_X = 0.0f;
		// Cosine demodulation
		for (j = 0; j < pfb->N; j++)
		{
			if (pfb->decimationCounter[j] == pfb->Sk[j])
				Curr_X += pfb->channelMatrix[j * pfb->N2 + i] * pfb->subbandData[j];
		}
		//printf("%1.8f\n", Curr_X[i]);
		// Polyphase filtering
		for (j = 0; j < pfb->m; j++)
		{
			pfb->Xk2[j * pfb->N2 + i] = Curr_X * pfb->h[j * pfb->N2 + i];
		}
	}
	//for (i = 0; i < L; i++)
	//	printf("%1.8f\n", Xk2[i]);
	// Dewarping network
	// Complexity: O(n)
	float R1_b = pfb->Xk2[0];
	for (i = 0; i < pfb->L - 1; i++)
	{
		float R3_b = pfb->APC_delay_1[i];
		pfb->APC_delay_1[i] = R1_b;
		float R2_b = pfb->APC_delay_2[i];
		R1_b = (R2_b - R1_b) * pfb->alpha + R3_b;
		pfb->APC_delay_2[i] = R1_b;

		R1_b = R1_b + pfb->Xk2[i + 1];
	}
	// Subsampling counter
	for (i = 0; i < pfb->N; i++)
	{
		if (pfb->decimationCounter[i] >= pfb->Sk[i])
			pfb->decimationCounter[i] = 1;
		else
			pfb->decimationCounter[i]++;
	}
	return R1_b * pfb->postGain;
}
void analysisWarpedPFBStereo(WarpedPFB *pfb1, WarpedPFB *pfb2, float *x1, float *x2)
{
	unsigned int i, j;
	// Warping network
	// Complexity: O(n)
	float R1_a = *x1 + pfb1->noiseBuffer[pfb1->noiseLoop];
	float R1_a_x2 = *x2 + pfb1->noiseBuffer[pfb1->noiseLoop];
	pfb1->noiseLoop = (pfb1->noiseLoop + 1) & (DENORMAL_BUFFER - 1);
	float R2_a = pfb1->allpass_delay_chain[0];
	float R2_a_x2 = pfb2->allpass_delay_chain[0];
	pfb1->allpass_delay_chain[0] = R1_a;
	pfb2->allpass_delay_chain[0] = R1_a_x2;
	for (i = 1; i < pfb1->L; i++)
	{
		float R3_a = R2_a;
		float R3_a_x2 = R2_a_x2;
		R2_a = pfb1->allpass_delay_chain[i];
		R2_a_x2 = pfb2->allpass_delay_chain[i];
		R1_a = (R2_a - R1_a) * pfb1->alpha + R3_a;
		R1_a_x2 = (R2_a_x2 - R1_a_x2) * pfb1->alpha + R3_a_x2;
		pfb1->allpass_delay_chain[i] = R1_a;
		pfb2->allpass_delay_chain[i] = R1_a_x2;
	}
	// Analysis, multiply by the coefficients of the prototype filter
	// Complexity: O(2 * m * M + 2 * M * sum(1 ./ Sk))
	memset(pfb1->subbandData, 0, pfb1->N * sizeof(float));
	memset(pfb2->subbandData, 0, pfb1->N * sizeof(float));
	for (i = 0; i < pfb1->N2; i++)
	{
		float Xk = 0.0f;
		float Xk_x2 = 0.0f;
		for (j = 0; j < pfb1->m; j++)
		{
			Xk += pfb1->allpass_delay_chain[i + j * pfb1->N2] * pfb1->h[i + j * pfb1->N2];
			Xk_x2 += pfb2->allpass_delay_chain[i + j * pfb1->N2] * pfb1->h[i + j * pfb1->N2];
		}
		//printf("%1.8f\n", Xk[i]);
		// Cosine modulation
		for (j = 0; j < pfb1->N; j++)
		{
			if (pfb1->decimationCounter[j] == pfb1->Sk[j])
			{
				pfb1->subbandData[j] += pfb1->channelMatrix[j * pfb1->N2 + i] * Xk;
				pfb2->subbandData[j] += pfb1->channelMatrix[j * pfb1->N2 + i] * Xk_x2;
			}
		}
	}
	for (j = 0; j < pfb1->N; j++)
	{
		pfb1->subbandData[j] *= pfb1->Sk[j];
		pfb2->subbandData[j] *= pfb1->Sk[j];
	}
}
void synthesisWarpedPFBStereo(WarpedPFB *pfb1, WarpedPFB *pfb2, float *y1, float *y2)
{
	unsigned int i, j;
	// Synthesis
	// Complexity same as analysis
	for (i = 0; i < pfb1->N2; i++)
	{
		float Curr_X = 0.0f;
		float Curr_X_x2 = 0.0f;
		// Cosine demodulation
		for (j = 0; j < pfb1->N; j++)
		{
			if (pfb1->decimationCounter[j] == pfb1->Sk[j])
			{
				Curr_X += pfb1->channelMatrix[j * pfb1->N2 + i] * pfb1->subbandData[j];
				Curr_X_x2 += pfb1->channelMatrix[j * pfb1->N2 + i] * pfb2->subbandData[j];
			}
		}
		//printf("%1.8f\n", Curr_X[i]);
		// Polyphase filtering
		for (j = 0; j < pfb1->m; j++)
		{
			pfb1->Xk2[j * pfb1->N2 + i] = Curr_X * pfb1->h[j * pfb1->N2 + i];
			pfb2->Xk2[j * pfb1->N2 + i] = Curr_X_x2 * pfb1->h[j * pfb1->N2 + i];
		}
	}
	//for (i = 0; i < L; i++)
	//	printf("%1.8f\n", Xk2[i]);
	// Dewarping network
	// Complexity: O(n)
	float R1_b = pfb1->Xk2[0];
	float R1_b_x2 = pfb2->Xk2[0];
	for (i = 0; i < pfb1->L - 1; i++)
	{
		float R3_b = pfb1->APC_delay_1[i];
		float R3_b_x2 = pfb2->APC_delay_1[i];
		pfb1->APC_delay_1[i] = R1_b;
		pfb2->APC_delay_1[i] = R1_b_x2;
		float R2_b = pfb1->APC_delay_2[i];
		float R2_b_x2 = pfb2->APC_delay_2[i];
		R1_b = (R2_b - R1_b) * pfb1->alpha + R3_b;
		R1_b_x2 = (R2_b_x2 - R1_b_x2) * pfb1->alpha + R3_b_x2;
		pfb1->APC_delay_2[i] = R1_b;
		pfb2->APC_delay_2[i] = R1_b_x2;

		R1_b = R1_b + pfb1->Xk2[i + 1];
		R1_b_x2 = R1_b_x2 + pfb2->Xk2[i + 1];
	}
	// Subsampling counter
	for (i = 0; i < pfb1->N; i++)
	{
		if (pfb1->decimationCounter[i] >= pfb1->Sk[i])
			pfb1->decimationCounter[i] = 1;
		else
			pfb1->decimationCounter[i]++;
	}
	*y1 = R1_b * pfb1->postGain;
	*y2 = R1_b_x2 * pfb1->postGain;
}