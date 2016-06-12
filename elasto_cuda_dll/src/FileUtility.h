//////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "opencv/cv.h"

int  ReadRFData(const char *file_path, float *rf, int rows, int cols);
int  ReadRFDataT(const char *file_path, short *rf, int rows, int cols);
int  ReadRFDataB(const char *file_path, short *rf, int rows, int cols);
int  ReadRFDataB(const char *file_path, float *rf, int rows, int cols);//add wangxiaomeng
int  ReadMatFile(const char *file_path, float *rf, int rows, int cols);

void  SaveDataFile(const char *filename, CvMat *pmat);

void  MakeBmpAndShow(const char *filename, const CvMat *pmat);
