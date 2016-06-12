// test_elasto_cuda_dll.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#include "CElasto.h"
#include "SysConfig.h"
#include "devMatchCuda.cuh"
#include "elasto_cuda.h"
#include <time.h>
#include <conio.h>
#include <fstream>
  

// #pragma   comment（lib，"elasto_cuda_dll.lib"）

int  readData(const char *file_path, float *rf, int rows, int cols)
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
				rf[i * cols + j] = n;
			}
		}
	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////
// 读取二进制RF数据， 数据元素，short
//
//////////////////////////////////////////////////////////////////////////
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
				rf[i * cols + j] = n;
			}
		}
		fclose(file);
	}

	return 0;
}


int  readMatFile(const char *file_path, float *rf, int rows, int cols)
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

int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "DEMO Elasto Graphy" << std::endl;

	ElastoInit();

	// read RF data

	EInput  in;

	EOutput out;

	in.filepath_d = "displace.bmp";
	in.filepath_s = "strain.bmp";

#if 1

	// 培田采集的数据
	/*
	const char *rf_filepath = "rf.dat";
	const int  line_num = 128;
	const int  sample_num_per_line = 2048;
	*/
	const char *rf_filepath = "rf.dat";
	in.rows = 300;
	in.cols = 4096 * 2;
	

	in.pDatas = new float[in.rows * in.cols];
	int ok = readData(rf_filepath, in.pDatas, in.rows, in.cols);

#else

	// 原来的数据
	const char *rf_filepath = "ph1.mat";
	in.rows = 300;
	in.cols = 4096 * 2;
	in.pDatas = new float[in.rows * in.cols];
	int ok = readMatFile(rf_filepath, in.pDatas, in.rows, in.cols);

#endif 

	clock_t start, finish;
	double total;
	while (1)
	{
		printf("Please Select:\n");
		printf("\t1: Do Demo\n");
		printf("\tq or Q: quit!\n");
		bool  exit = false;
        int  ch = _getche();



		switch (ch)
		{
		case '1':

			if (!initCUDA())  {                              //cpu  platform   
				// get strain image and modulus
				start = ::clock();




				ok = ElastoProcess(in, out);

				finish = ::clock();
				total = (double)(finish - start) / CLOCKS_PER_SEC;
				printf("\nTotalTime is %fs!\n", total);
				if (ok == 0)
				{
					//std::cout << "E = " << e << " kPa" << std::endl;
					printf("\tE=%fkPa,V=%fm/s\n", out.e, out.v);
				}
				break;

			}

			else  {                                      //gpu  platform 

				ElastoCuda * testCuda = new ElastoCuda();


				testCuda->process(in, out);
			
			
			
			
			}


		case 'q':
		case 'Q':
			exit = true;
			break;

		default:
			break;
		}
		if (exit) break;
	}

	ElastoRelease();

#if 0

	// display Strain Image 
	IplImage * image = cvLoadImage(in.filepath_d, CV_LOAD_IMAGE_UNCHANGED);
	assert(image);

	const char * window_name = "Strain Image";// 窗口的名字&标题

	cvNamedWindow(window_name);
	cvShowImage(window_name, image);
	cvWaitKey(0);
	cvDestroyWindow(window_name);
	cvReleaseImage(&image);

#else

	//printf("Press Any Key to Quit!\n");
	getchar();

#endif

	ElastoRelease();

	return 0;
}

