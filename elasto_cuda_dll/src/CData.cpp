#include "stdafx.h"
#include "CData.h"
#include <fstream>
#include <iostream>
//#include "highgui.h"

CData::CData(int rows = 0, int cols = 0)
{
	dataMat = cvCreateMat(rows, cols, CV_32FC1);	//initiating the data
}

void CData::readData(float *input)
{
	for (int i = 0; i < dataMat->rows; i++)
	{
		for (int j = 0; j < dataMat->cols; j++)
		{
			*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*dataMat, i, j)))) = input[i * dataMat->cols + j];
		}
	}
}






CvMat* CData::getData()
{
	return dataMat;
}

CvMat *CData::getSubData(int x, int y, int w, int h)
{
	CvMat *pmat = cvCreateMat(h, w, dataMat->type);

	CvMat  sub_mat;
	cvGetSubRect(dataMat, &sub_mat, cvRect(x, y, w, h));

	cvCopy(&sub_mat, pmat);

	return pmat;
}