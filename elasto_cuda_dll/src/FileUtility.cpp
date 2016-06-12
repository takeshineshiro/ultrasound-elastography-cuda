//////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "FileUtility.h"
#include "opencv/highgui.h"

#include <time.h>
#include <conio.h>
#include <fstream>
#include <io.h>

int  ReadRFData(const char *file_path, float *rf, int rows, int cols)
{
	FILE *file = fopen(file_path, "r");
	int  ok = -1;
	if (file)
	{
		ok = 0;
		int n;
		int ret;
		int i;
		int j;
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				ret = fscanf(file, "%d", &n);
				if (ret == EOF) break;
				rf[i * cols + j] = (float)n;
			}
		}

		fclose(file);
	}

	return 0;
}


int  ReadRFDataT(const char *file_path, short *rf, int rows, int cols)
{
	FILE *file = fopen(file_path, "r");
	int  ok = -1;
	if (file)
	{
		ok = 0;
		int n;
		int ret;
		int i;
		int j;
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				ret = fscanf(file, "%d", &n);
				if (ret == EOF) break;
				rf[i * cols + j] = (short)n;
			}
		}

		fclose(file);
	}

	return 0;
}

//creator wangxiaomeng
//date 20160131
int   ReadRFDataB(const char *file_path, float *rf, int rows, int cols)
{
	FILE *file = fopen(file_path, "rb");
	int  ok = -1;
	if (file)
	{
		ok = 0;
		short n;
		int ret;
		int i;
		int j;
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				ret = fread(&n, sizeof(n), 1, file);
				if (ret == 0) break;
				rf[i * cols + j] = (float)n; 
			}
		}
		fclose(file);
	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////
// 读取二进制RF数据， 数据元素，float
//
//////////////////////////////////////////////////////////////////////////
int   ReadRFDataB(const char *file_path, short *rf, int rows, int cols)
{
	FILE *file = fopen(file_path, "rb");
	int  ok = -1;
	if (file)
	{
		ok = 0;
		short n;
		int ret;
		int i;
		int j;
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				ret = fread(&n, sizeof(n), 1, file);
				if (ret == 0) break;
				rf[i * cols + j] = n;
			}
		}
		fclose(file);
	}

	return 0;
}


int  ReadMatFile(const char *file_path, float *rf, int rows, int cols)
{
	std::fstream infile;
	infile.open(file_path, std::ios_base::out | std::ios_base::in | std::ios_base::binary);

	char matheader[124];                 //mat file header, 124bytes
	infile.read(matheader, 124);

	//short version;
	char tmpver[4];
	infile.read(tmpver, 4);          //version, 2bytes; //endian indicator, 2bytes

	// #ifdef _DEBUG
	// 	version = *static_cast<short*>(static_cast<void*>(tmpver));
	// 	std::cout << version << std::endl;
	// 	version = *static_cast<short*>(static_cast<void*>(tmpver + 2));
	// 	std::cout << version << std::endl;
	// #endif //_DEBUG
	// 
	// 	int datanum;
	// 	char tmpnum[8];
	// 	infile.read(tmpnum, 8);            //datatype, 4bytes; datanumber, 4bytes
	// #ifdef _DEBUG
	// 	datanum = *static_cast<int*>(static_cast<void*>(tmpnum));
	// 	std::cout << datanum << std::endl;
	// 	datanum = *static_cast<int*>(static_cast<void*>(tmpnum + 4));
	// 	std::cout << datanum << std::endl;
	// #endif //_DEBUG

	char dataheader[56];
	infile.read(dataheader, 56);        //unknown header, 56bytes

	char    data_ch[8];
	double  f;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			infile.read(data_ch, 8);
			f = *(static_cast<double*>(static_cast<void*>(data_ch)));
			rf[i * cols + j] = static_cast<float>(f);
		}
	}
	infile.close();
	return 0;
}


void  SaveDataFile(const char *filename, CvMat *pmat)
{
	FILE *file = NULL;
	errno_t err = fopen_s(&file, filename, "wt");
	if (err == 0)
	{
		int i = 0;
		int j = 0;
		for (i = 0; i < pmat->rows; i++)
		{
			for (j = 0; j < pmat->cols; j++)
			{
				fprintf_s(file, "%f ", CV_MAT_ELEM(*pmat, float, i, j));
			}
			fprintf_s(file, "\n");
		}

		fclose(file);
	}
}


void  MakeBmpAndShow(const char * filename, const CvMat *pmat)
{
	IplImage *pimage = cvCreateImage(cvGetSize(pmat), IPL_DEPTH_32F, 3);

	cvCvtColor(pmat, pimage, CV_GRAY2BGR);
	cvNamedWindow(filename, CV_WINDOW_AUTOSIZE);
	cvShowImage(filename, pimage);  
	cvWaitKey(0);
	cvDestroyWindow(filename);

	cvSaveImage(filename, pimage);

	cvReleaseImage(&pimage);
}