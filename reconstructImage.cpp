#include "mex.h"

#include <cmath>
#include <cstring>
#include <algorithm>

#define IN_PATCHIMAGE prhs[0]
#define IN_LOCATIONS  prhs[1]
#define IN_PATCHSIZE prhs[2]
#define IN_DATAM prhs[3]
#define IN_DATAN prhs[4]
#define IN_FRAMEM prhs[5]
#define IN_FRAMEN prhs[6]
#define IN_RECORDLEN prhs[7]

#define OUT_IMAGE plhs[0]

int partition(double *values, int left, int right)
{
	int pos = right;
	double temp;
	right--;
	while (left <= right)
	{
		while (left < pos && values[left] <= values[pos])	left++;
		while (right >= 0 && values[right] > values[pos])	right--;
		if (left >= right)		break;
		temp = values[left]; values[left] = values[right];	values[right] = temp;
	}
	temp = values[left]; values[left] = values[pos]; values[pos] = temp;
	return left;
}

double getMidIndex(double *values, int size)
{
	int left = 0;
	int right = size - 1;
	int midPos = right >> 1;
	int index = -1;
	while (index != midPos)
	{
		index = partition(values, left, right);

		if (index < midPos)
			left = index + 1;
		else if (index > midPos)
			right = index - 1;
		else
			break;
	}
//	assert(index == midPos);
	return values[index];
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	double *pPatchImage = mxGetPr(IN_PATCHIMAGE);
	double *pLocations = mxGetPr(IN_LOCATIONS);
	int iPatchSize = mxGetScalar(IN_PATCHSIZE);
	int m = mxGetScalar(IN_DATAM), n = mxGetScalar(IN_DATAN);
	int row = mxGetScalar(IN_FRAMEM), col = mxGetScalar(IN_FRAMEN);
	int iRecordTolalLen = mxGetScalar(IN_RECORDLEN);
	int iPatchLen = iPatchSize*iPatchSize;

	int iFramePerPatch = m / iPatchLen;
	int iFrameLen = row*col;
	int iTotalLen = iFrameLen*iFramePerPatch;
	
	mwSize dims[3];
	dims[0] = row, dims[1] = col, dims[2] = iFramePerPatch;
	OUT_IMAGE = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    
	//plhs[1] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);

	double *pImage = mxGetPr(OUT_IMAGE);
	memset(pImage, 0, sizeof(double)*iTotalLen);

	// allocate memory
	int i, j, k, x, y, c, offset;
	double **pRecord = new double*[iTotalLen];
	for (i = 0; i < iTotalLen; i++)
	{
		pRecord[i] = new double[iRecordTolalLen];
	}
	int *pRecordLen = new int[iTotalLen];
	memset(pRecordLen, 0, sizeof(int)*iTotalLen);

	// assign each value into corresponding address
	double *pPatchImage1, *pPatchImage2, *pPatchImage3, *pPachImage4;
	double **pRecord1, **pRecord2, **pRecord3, **pRecord4;
	int  *pRecordLen1, *pRecordLen2, *pRecordLen3, *pRecordLen4;
	double *pImage1, *pImage2, *pImage3;
	for (pPatchImage1 = pPatchImage,
		c = 0; c < n; c++, // the column of patch image
		pPatchImage1 +=m)
	{
		x = pLocations[c+n];
		y = pLocations[c];
		offset = x*row + y;
		for (pPatchImage2 = pPatchImage1, pRecord1 = pRecord + offset, pRecordLen1 = pRecordLen+offset,
			k = 0; k < iFramePerPatch; k++,  // the temporal length of 3D patch image
			pPatchImage2+=iPatchLen, pRecord1 += iFrameLen, pRecordLen1 += iFrameLen)
		{
			for (pPatchImage3 = pPatchImage2, pRecord2 = pRecord1, pRecordLen2 = pRecordLen1,
				j = 0; j < iPatchSize; j++,  // the row of one 2D patch image
				pPatchImage3 += iPatchSize, pRecord2 += row, pRecordLen2 += row)
			{
				for (pPachImage4 = pPatchImage3, pRecord3 = pRecord2, pRecordLen3 = pRecordLen2,
					i = 0; i < iPatchSize; i++,  // one column of one 2D patch image
					pPachImage4++, pRecord3++, pRecordLen3++)
				{
					(*pRecord3)[*pRecordLen3] = *pPachImage4;
					(*pRecordLen3) = (*pRecordLen3)+1;
					
				}
	
			}

		}
	}

	// convert into image by medium filter
	for (pImage1=pImage, pRecord1 = pRecord, pRecordLen1 = pRecordLen,
		i = 0; i < iTotalLen; i++, pImage1++, pRecord1++, pRecordLen1++)
	{
		if (*pRecordLen1 > 0)
		{
			if ((*pRecordLen1) % 2 == 0)
			{
				k = (*pRecordLen1) / 2;
				std::sort(*pRecord1, (*pRecord1) + (*pRecordLen1));
				*pImage1 = ((*pRecord1)[k - 1] + (*pRecord1)[k])/2.0;

			}
			else
				*pImage1 = getMidIndex(*pRecord1, *pRecordLen1);
		} 			
	}

	// release memory
	for (i = 0; i < iTotalLen; i++)
	{
		delete[] pRecord[i];
	}
	delete[] pRecord;
	delete[] pRecordLen;
    return;
}
