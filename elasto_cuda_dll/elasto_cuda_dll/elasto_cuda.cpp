#include "stdafx.h"
#define ELASTO__EXPORTS   

#ifndef   ELASTO_CUDA
#include "elasto_cuda.cuh"
#endif

#include "stdio.h"

#include  "lib_add_cuda.cuh"

#include  "CElasto.h"  

#include  "SysConfig.h"

#include "CData.h"

#include  "cuda_main.cuh"



ElastoCuda::ElastoCuda()  {

	initFile   = new  std::string();

	cudaMain   = new  CudaMain();

}





ElastoCuda::~ElastoCuda()  {

	delete    initFile;

	delete    cudaMain;

}



bool  ElastoCuda::isAvailable()  {

	return    cudaMain->isAvailable();

}



void  ElastoCuda::init(const std::string &ini_file) {

	*initFile = ini_file;

	ReadSysConfig( initFile,  *config);                                             //读取文件获取参数


	cudaMain->inputConfigParam(config);                                            // 配置参数


	cudaMain->getFilterParam( config->lpfilt_file);                                // 读取滤波器低通参数


	cudaMain->getFilterParam(config->bpfilt_file);                                // 读取滤波器带通参数


	cudaMain->getFilterParam(config->matchfilt_file);                             // 读取滤波器匹配参数

	cudaMain->mallocMem();                                                       //  分配内存


}


void   ElastoCuda::readRFData(const EInput & in)   {

	std::string filename;
	filename = in.filepath_s;

	CData*test = new CData(in.rows, in.cols);
	
	test->readData(in.pDatas);

}










bool  ElastoCuda::process(const EInput &input, EOutput &output)  {

	    bool       

	cudaMain->process(input, output);





}

















