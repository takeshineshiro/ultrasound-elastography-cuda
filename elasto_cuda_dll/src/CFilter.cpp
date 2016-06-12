#include "stdafx.h"
#include "CFilter.h"
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream>



 void SaveDispDataFileFilter(const char *filename, CvMat *pmat);

CFilter::CFilter(std::string filename)
{
	if (filename.size() == 0)
	{
		return;
	}

	std::fstream paramFile(filename.c_str());
	if (!paramFile)
	{
		return;
	}

	std::stringstream ss;
	float tmp;
	std::string str;

	//for (vector<double>::size_type ix = 0; ix != param.size(); ++ix)
	param.clear();
	while (!paramFile.eof())
	{
		ss.clear();
		paramFile >> str;
		//str.erase(str.length()-1);
		ss << str;
		ss >> tmp;
		param.push_back(tmp);
	}
	steps = param.size();
	paramFile.close();
	ss.clear();
}

void CFilter::Do()
{
	filterAlgorithm();
}



void CFilter::DoTimeFieldFilter(const CvMat *pSrc)
{
	double tmp = 0;

	const CvMat *pmat = (pSrc) ? pSrc : outDataMat;


	for (int i = 0; i < pmat->cols; ++i)
	{
		for (int j = 0; j < pmat->rows; ++j)
		{

			for (int k = 0; k <= j; ++k)
			{
				tmp += cvmGet(pmat, k, i) * param[(j - k) < steps ? (j - k) : 0];  /*卷积*/
			}


			cvmSet(outDataMat, j, i, tmp);



			tmp = 0;

		}
	}


}




void CFilter::filterAlgorithm()

{

	CvMat* tmpMat = cvCreateMat(outDataMat->rows, outDataMat->cols, outDataMat->type);

	cvCopy(outDataMat, tmpMat);

	float filttmp = 0.0;


	//正滤波
	for (int k = 0; k != outDataMat->rows; ++k)      //针对每根扫描线
	{
		for (int i = 0; i != steps; ++i)         //数据长度小于等于滤波器抽头长度
		{
			filttmp = 0.0;

			for (int j = 0; j <= i; ++j)
			{
				filttmp += param[j] * (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*tmpMat, k, i - j)))));
			}
			for (int j = i + 1; j < steps; j++)
			{
				filttmp += param[j] * (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*tmpMat, k, 0)))));
				//filttmp += param[j] * 0;
			}
			*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*outDataMat, k, i)))) = filttmp;
		}


		for (int i = steps; i != outDataMat->cols; ++i)        //数据长度大于滤波器抽头长度
		{
			filttmp = 0.0;
			for (std::vector<float>::size_type ix = 0; ix != steps; ++ix)
			{
				filttmp += param[ix] * (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*tmpMat, k, i - ix)))));
				//float tmp = (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*outDataMat, k, i - ix)))));
			}
			*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*outDataMat, k, i)))) = filttmp;
			//cvSetReal2D(outDataMat, k, i, static_cast<float>(filttmp));
		}
	}

	SaveDispDataFileFilter("positivefilterData.dat", outDataMat);

	//逆滤波，wxm注释掉了，为了测试20160225

	cvCopy(outDataMat, tmpMat);

	for (int k = 0; k != outDataMat->rows; ++k)     //针对每根扫描线
	{


		for (int i = 0; i != steps; ++i)         //数据长度小于等于滤波器抽头长度
		{

			filttmp = 0.0;

			for (int j = 0; j <= i; ++j)
			{
				filttmp += param[j] * (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*tmpMat, k, outDataMat->cols - i - 1 + j)))));

			}
			for (int j = i + 1; j < steps; j++)

			{
				filttmp += param[j] * (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*tmpMat, k, outDataMat->cols - 1)))));
				//filttmp += param[j] * 0;
			}
			*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*outDataMat, k, outDataMat->cols - i - 1)))) = filttmp;
		}



		for (int i = steps; i != outDataMat->cols; ++i)   //数据长度大于滤波器抽头长度
		{
			filttmp = 0.0;

			for (std::vector<double>::size_type ix = 0; ix != steps; ++ix)
				//for (std::vector<double>::size_type ix = 0; ix != param.size(); ++ix)
			{
				filttmp += param[ix] * (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*tmpMat, k, outDataMat->cols - i - 1 + ix)))));
				//filttmp += param[ix] * (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*tmpMat, k, outDataMat->cols - i + ix)))));
				//float tmp = (*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*outDataMat, k, i - ix)))));
			}
			*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*outDataMat, k, outDataMat->cols - i - 1)))) = filttmp;
			//*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*outDataMat, k, outDataMat->cols - i )))) = filttmp;
			//cvSetReal2D(outDataMat, k, outDataMat->cols - i - 1, static_cast<float>(filttmp));
		}



	}


	SaveDispDataFileFilter("negativefilterData.dat", outDataMat);


	cvReleaseMat(&tmpMat);


	// #ifdef _DEBUG
	// 	FILE *outfile;    //建立文件流存储位移
	// 	if(fopen_s(&outfile, "dis50overlap90f.dat", "wb")) std::cout << "can not open file" << std::endl;
	// 	double fordis=0;
	// 	for(int i = 0; i < outDataMat->rows; ++i)
	// 	{
	// 		for(int j = 0; j < outDataMat->cols; ++j)
	// 		{
	// 			fordis = *(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*outDataMat, i, j))));
	// 			fwrite(&fordis, sizeof(double), 1, outfile);
	// 		}
	// 	}
	// 	fclose(outfile);
	// #endif //_DEBUG
}




void  SaveDispDataFileFilter(const char *filename, CvMat *pmat)
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