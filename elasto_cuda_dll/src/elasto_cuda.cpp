
#include "stdafx.h"

#define ELASTO__EXPORTS   

#include "elasto_cuda.h"

#include "stdio.h"

#include  "lib_add_cuda.cuh"

#include  "CElasto.h"  

#include  "SysConfig.h"

#include "CData.h"

#include  "cuda_main.cuh"





ElastoCuda::ElastoCuda()  {

	initFile   = defaultElastoConfFile;

	  config   = new  ConfigParam;

	cudaMain   = new  CudaMain();

}





ElastoCuda::~ElastoCuda()  {

	delete    initFile;

	delete    cudaMain;

}



bool  ElastoCuda::isAvailable()  {

	return    cudaMain->isAvailable();

}



void  ElastoCuda::init(const EInput & in) {


	initFile = defaultElastoConfFile;


	ReadSysConfig( initFile,  *config);                                            // 读取文件获取参数


//	readRFData(in);                                                                //  读取rf数据

	

//	cudaMain->inputRfData( in);                                                    //  读取RF数据,给cpu_inputMat 

	cudaMain->inputConfigParam(config);                                            //  配置参数 给cpu_config


	cudaMain->getlowFilterParam( config->lpfilt_file);                            // 读取滤波器低通参数，给


	cudaMain->getbandFilterParam(config->bpfilt_file);                           // 读取滤波器带通参数，给


	cudaMain->getmatchFilterParam(config->matchfilt_file);                        // 读取滤波器匹配参数，给

	

	cudaMain->mallocMem();                                                         //  分配内存


	cudaMain->inputRfData(in);                                                    //  读取RF数据,给cpu_inputMat 


}


void   ElastoCuda::readRFData(const EInput & in)   {

	std::string filename;

	filename = in.filepath_s;

	CData*test = new CData(in.rows, in.cols);
	
	test->readData(in.pDatas);


}










bool  ElastoCuda::process(const EInput &input, EOutput &output)  {

	init(input);
	   
	cudaMain->process(input, output);



	return  true;


}

















