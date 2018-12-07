#include "mex.h"

#include <cmath>
#include <cstring>

#define IN_LOGRHO prhs[0]
#define IN_R      prhs[1]
#define IN_LAMBDA prhs[2]
#define IN_DATAM prhs[3]
#define IN_DATAN prhs[4]
#define IN_FRAMEM prhs[5]
#define IN_FRAMEN prhs[6]
#define IN_K prhs[7]
#define IN_ITER prhs[8]

#define OUT_LOGRHO plhs[0]
#define OUT_R plhs[1]

double logsumexp(double a[], int k)
{
	double amax = a[0];
	int i;
	for (i = 1; i < k; i++)
	{
		amax = amax>a[i] ? amax : a[i];
	}
	double s = 0;
	for (i = 0; i < k; i++ )
	{
		s = s + exp(a[i] - amax);
	}
	s = log(s) + amax;
    return s;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	//mexPrintf("Hello World!\n");
	double *oldLogRho = mxGetPr(IN_LOGRHO);
	double *oldR = mxGetPr(IN_R);
	double lambda = mxGetScalar(IN_LAMBDA);
	int m = mxGetScalar(IN_DATAM), n = mxGetScalar(IN_DATAN);
	int row = mxGetScalar(IN_FRAMEM), col = mxGetScalar(IN_FRAMEN);
	int kMog = mxGetScalar(IN_K);
	int iterN = mxGetScalar(IN_ITER);

	OUT_LOGRHO = mxCreateNumericArray(mxGetNumberOfDimensions(IN_LOGRHO), mxGetDimensions(IN_LOGRHO), mxDOUBLE_CLASS, mxREAL);
	OUT_R = mxCreateNumericArray(mxGetNumberOfDimensions(IN_R), mxGetDimensions(IN_R), mxDOUBLE_CLASS, mxREAL);

	double *newLogRho = mxGetPr(OUT_LOGRHO);
	double *newR = mxGetPr(OUT_R);

	double zNeighborSum, v;
	double tempLogSumExpRho;
	double tempLogRho[100];
	int iChannelLen = m * n;
	int iTolalLen = iChannelLen * kMog;
	int iFramelen = row * col;

	memcpy(newLogRho, oldLogRho, sizeof(double) * iTolalLen);
	memcpy(newR, oldR, sizeof(double) * iTolalLen);
	
	int i, j, k, temp;
	double *pR1, *pR2, *pR3, *pR4;
	double *pNewRho1, *pNewRho2, *pNewRho3, *pNewRho4;
	double *pOldRho1, *pOldRho2, *pOldRho3, *pOldRho4;
	int nFramePerMog = iTolalLen / (iFramelen*kMog);
	for (int iter = 0; iter < iterN; iter++)
	{
		// compute neighbor accumulation and save into dAccu; 
		for (pR1 = newR, pNewRho1 = newLogRho, pOldRho1 = oldLogRho,
			n = 0; n < nFramePerMog; n++,
			pR1 += iFramelen, pNewRho1 += iFramelen, pOldRho1 += iFramelen)
		{
			// computer inner pixels
			for (pR2 = pR1, pNewRho2 = pNewRho1 , pOldRho2 = pOldRho1,
				j = 0; j < col ; j++,
				pR2 += row, pNewRho2 += row, pOldRho2 += row)
				for (pR3 = pR2, pNewRho3 = pNewRho2, pOldRho3 = pOldRho2,
					i = 0; i < row; i++,
					pR3++, pNewRho3++, pOldRho3++)
				{
					for (pR4 = pR3, pNewRho4 = pNewRho3, pOldRho4 = pOldRho3,
						k = 0; k < kMog; k++,
						pR4 += iChannelLen, pNewRho4 += iChannelLen, pOldRho4 += iChannelLen)
					{
						// 4-neigbor
						v = i == 0 ? 0 : *(pR4 - 1);
						v = i == row - 1 ? v : v + *(pR4 + 1);
						v = j == 0 ? v : v + *(pR4 - row);
						v = j == col - 1 ? v : v + *(pR4 + row);
					//	v = *(pR4 - 1) + *(pR4 + 1) + *(pR4 - row) + *(pR4 + row); 
						*pNewRho4 = *pOldRho4 + lambda * v;
						tempLogRho[k] = *pNewRho4;
					}
					tempLogSumExpRho = logsumexp(tempLogRho, kMog);
					for (pR4 = pR3, pNewRho4 = pNewRho3,
						k = 0; k < kMog; k++,
						pR4 += iChannelLen, pNewRho4 += iChannelLen)
					{
						*pR4 = exp(*pNewRho4 - tempLogSumExpRho);
					}
				}
		}
	}
    return;
}