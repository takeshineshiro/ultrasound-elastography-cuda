#include  "cuda_main.cuh"


#ifdef  _CUDA_MAIN_CUH
#include  "SysConfig.h"
#include  "CElasto.h"
#include  "FileUtility.h"

#include <fstream>
#include <string>
#include <iostream>
#include <time.h>
#include <device_launch_parameters.h>
#include <device_functions.h>
#include <math_functions.h>
#include <string.h>
#include <cstdio>

#endif


//内核函数及device函数


//需要更改，该款GPU芯片只能支持1024个threads per  block !!!      changed  by  wong   2016/06/08

//零相移数字滤波器    正滤波           changed  by  wong    2016/5/13

__global__ void Bandpass_front_1(Complex* tInput, int iWidth, float* param, int iParaLen, Complex* tOutPut)
{
	//float muData[40];

	//k  = blockIdx.x   i = blockIdx.y;
	float data_sum;

	float data_1;

	data_sum = 0.0;

	//	__shared__     float     data_sum[8192];



	/*for (int i = 0; i < iParaLen; i++)
	{
	muData[i] = *(param + i);
	}*/

	if (threadIdx.x <= iParaLen - 1)                                         //数据数目小于等于抽头数目
	{

		for (int i = 0; i <= threadIdx.x; i++)
		{


			data_1 = *(tInput + blockIdx.x*iWidth + threadIdx.x - i);      //b(0)*x(n-0)+b(1)*x(n-1)+...+b(n)*x(0)   

			data_sum += (data_1*param[i]);

		}

		data_1 = *(tInput + blockIdx.x * iWidth);                          // x(0)


		for (int j = threadIdx.x + 1; j <= iParaLen - 1; j++)
		{
			data_sum += (data_1*param[j]);                                 //b(n+1)*x(0)+...+b(nb-2)*x(0)
		}

		*(tOutPut + blockIdx.x * iWidth + threadIdx.x) = data_sum;



	}
	else                                                                  //数据数目大于抽头数目            
	{
		//data_1 = (tInput + blockIdx.x*iWidth + blockIdx.y - threadIdx.x)->x;
		for (int i = 0; i <= iParaLen - 1; i++)
		{

			data_1 = *(tInput + blockIdx.x*iWidth + threadIdx.x - i);   //b(0)*x(n-0)+b(1)*x(n-1)+...+b(nb-2)*x(n-(nb-2))  

			data_sum += (data_1*param[i]);

		}

		*(tOutPut + blockIdx.x * iWidth + threadIdx.x)= data_sum;

	}

}



//零相移数字滤波器   逆滤波                          changed  by  wong     2016/5/13
__global__ void Bandpass_back_1(Complex* tInput, int iWidth, float* param, int iParaLen, Complex* tOutPut)
{
	/*
	if (threadIdx.x < iParaLen - 1)                                     //  n  >= nb-1
	{
	return;                                                            //  此处有问题     changed  by  wong
	}
	//   changed   by  wong     2016/5/11
	*/                                                                 // 此处没有考虑小于nb-1的情况，需要保留输出


	float data_1;

	float data_sum;

	data_sum = 0.0;



	if (threadIdx.x <= iParaLen - 1)   {                                //  数据长度小于等于滤波器抽头长度


		for (int i = 0; i <= threadIdx.x; i++)
		{


			data_1 = *(tInput + blockIdx.x*iWidth + iWidth - 1 - threadIdx.x + i);      //

			data_sum += (data_1*param[i]);

		}


		data_1 = *(tInput + blockIdx.x * iWidth + iWidth - 1);                        //  x(N-1) 


		for (int j = threadIdx.x + 1; j <= iParaLen - 1; j++)
		{

			data_sum += (data_1*param[j]);                                            // b(n+1)*x(N-1)+...+b(nb-1)*x(N-1) 

		}


		*(tOutPut + blockIdx.x * iWidth + threadIdx.x) = data_sum;                 // y(n)

	}

	else    {                                                                       //数据长度大于滤波器抽头长度         


		for (int i = 0; i <= iParaLen - 1; i++)
		{
			data_1 = *(tInput + blockIdx.x*iWidth + iWidth - 1 - threadIdx.x + i);

			data_sum += (data_1*param[i]);                                         //  y(N-1-n) = b(0)*x(N-1-n+0) +b(1)*x(N-1-n+1)+...+b(nb-1)*x(N-1-n+nb-1)

		}

		*(tOutPut + blockIdx.x * iWidth + iWidth - 1 - threadIdx.x)   = data_sum;




	}


}


//changed   by  wong      2016/5/13

__device__      void   xcorr_cuda(const  Complex* templateMat_startID, const Complex* objectMat_startID, Complex*resultMat_startID)     {


	for (int i = 0; i < 101; i++)   {


		Complex     sum_object = 0;

		Complex     frac_object = 0;


		Complex     pow_template = 0;


		Complex     pow_object = 0;


		Complex     result = 0;


		//sum_object 

		for (int j = 0; j < 100; j++)  {


			sum_object += *(objectMat_startID + i + j);


		}

		//  ave

		Complex   ave_object = sum_object / 200;


		//fraction

		for (int j = 0; j < 100; j++)  {

			Complex    tmp = *(templateMat_startID + j) *  (*(objectMat_startID + i + j) - ave_object);


			frac_object += tmp;

		}


		//pow   temp

		for (int j = 0; j < 100; j++)  {


			pow_template += *(templateMat_startID + j) * *(templateMat_startID + j);

		}

		//pow   objectMat 

		for (int j = 0; j < 100; j++)  {


			pow_object += *(objectMat_startID + i + j)* * (objectMat_startID + i + j);

		}

		//result

		result = sqrt(pow_template*pow_object);

		//output

		*(resultMat_startID + i) = frac_object / result;

	}


}


//changed   by  wong    2016/5/13

__device__      void   minMax_cuda(Complex*resultMat_startID, Complex* min_value, Complex*  max_value, int  max_location)   {

	int      max_loc_temp = 0;

	int      min_loc_temp = 0;

	Complex* max_temp     = 0;

	Complex* min_temp    = 0;

	//求最大值及位置

	for (int i = 0; i < 101; i++)  {

		if (*(resultMat_startID + i) >= *max_temp)  {

			*max_temp = *(resultMat_startID + i);

			max_loc_temp = i;

		}

	}

	//求最小值及位置

	for (int i = 0; i < 101; i++)  {

		if (*(resultMat_startID + i) <= *min_temp)  {

			*min_temp = *(resultMat_startID + i);

			min_loc_temp = i;
		}


	}

	//输出

	*min_value = *min_temp;

	*max_value = *max_temp;

	max_location = max_loc_temp;

}


 
//changed  by   wong      2016/5/17

__device__    void    interp_cuda(Complex*resultMat_startID, int  max_loc, Complex*max_value, int multiWin, int  winSize, Complex*  displace)     {

	Complex*pre  = resultMat_startID + max_loc - 1;

	Complex*next = resultMat_startID + max_loc + 1;


	*displace   = (multiWin - 1) * winSize / 2 - max_loc - (*pre - *next) / (2 * (*pre - 2 * *max_value + *next));


}





//输出为位移299*799矩阵

__global__   void  displacement_api_cuda(Complex*disInputCuda, int rows, int cols, int  multiWin, int winSize, int  stepSize, templateMat*templateMatShare, objectMat* objectMatShare, resultMat*resultMatShare, Complex*min, Complex*max, int*max_location, Complex* displacement )      {


	int   out_offset = blockIdx.x *blockDim.x + threadIdx.x;                     // 输出位移矩阵偏移值

	int    bid       = blockIdx.x ;                                              //  对应block ID 
	
	int    tid       = threadIdx.x;                                             //   对应thread ID 


   //考虑使用3D数组，因为共享内存不够用，只有49152个字节each block 

	//共享内存

	//改用全局内存！！！

//	__shared__     Complex*     templateMatShare[THREAD_NUM];        //100首地址

//	__shared__     Complex*     objectMatShare[THREAD_NUM];          //200首地址

//	__shared__     Complex*     resultMatShare[THREAD_NUM];          //101首地址 


//	__shared__     templateMat   templateMatShare[THREAD_NUM];

//	__shared__     objectMat     objectMatShare[THREAD_NUM];

//	__shared__     resultMat     resultMatShare[THREAD_NUM];


//	  Complex*templateMatShare;                          //   模板内存在GPU分配         考虑局部分配                  


//	  Complex*objectMatShare;                           //    目标内存在GPU分配         考虑局部分配


//	  Complex*resultMatShare;                           //    匹配结果在GPU分配         考虑局部分配




//	cudaMalloc(&templateMatShare, winSize* sizeof(Complex));             //模板矩阵


//	cudaMalloc(&objectMatShare, winSize*multiWin* sizeof(Complex));     //目标矩阵


//	cudaMalloc(&resultMatShare, (winSize + 1)* sizeof(Complex));        //结果矩阵




//	 templateMatShare[out_offset].elem   = (Complex*)(disInputCuda + blockIdx.x*cols + (multiWin - 1) * winSize / 2 + threadIdx.x * stepSize);
		







	    Complex*templateMatID;                               //ID


	  Complex*objectMatID;                                //ID



	//12784字节

//	__shared__     Complex*     min[THREAD_NUM];

//	__shared__     Complex*     max[THREAD_NUM];

//	__shared__      int        max_location[THREAD_NUM];

//	__shared__    Complex*     displacement[THREAD_NUM];




	//准备相关数据      线程块有问题 ，采用共享内存


	//	(templateMat*)(templateMat_startID + threadIdx.x)->elem = (Complex*)(disInputCuda + blockIdx.x*cols + (multiWin - 1) * winSize / 2 + threadIdx.x * stepSize);

	      templateMatID = (Complex*)(disInputCuda + blockIdx.x*cols + (multiWin - 1) * winSize / 2 + threadIdx.x * stepSize);


//		  cudaMemcpy(templateMatShare[out_offset].elem, templateMatID, winSize, cudaMemcpyDeviceToDevice);

		  for (int i = 0; i < 100; i++) {

			  templateMatShare[out_offset].elem[i] = *(templateMatID+i);

		  }

 // templateMatShare[threadIdx.x].elem = (Complex*)(disInputCuda + blockIdx.x*cols + (multiWin - 1) * winSize / 2 + threadIdx.x * stepSize);


		  objectMatID   = (Complex*)(disInputCuda + (blockIdx.x + 1)*cols + threadIdx.x * stepSize);


//		  cudaMemcpy(objectMatShare[out_offset].elem, objectMatID, winSize*multiWin, cudaMemcpyDeviceToDevice);

             
		  for (int i = 0; i < 200; i++)  {
		  
			  objectMatShare[out_offset].elem[i] = *(objectMatID+i);
		  
		  }
  


	//		objectMat_startID[threadIdx.x].elem   = disInputCuda + (blockIdx.x + 1)*cols + threadIdx.x * stepSize;



//	__syncthreads();

//	cudaThreadSynchronize();

//	cudaDeviceSynchronize();



	//相关运算

		  xcorr_cuda(templateMatShare[out_offset].elem, objectMatShare[out_offset].elem, resultMatShare[out_offset].elem);


//		__syncthreads();
//	cudaThreadSynchronize();


	//查找最大值

		  minMax_cuda(resultMatShare[out_offset].elem, &min[out_offset], &max[out_offset], max_location[out_offset]);


//	__syncthreads();

//	cudaThreadSynchronize();

	//插值

		  interp_cuda(resultMatShare[out_offset].elem, max_location[out_offset], &max[out_offset], multiWin, winSize, &displacement[out_offset]);


//	__syncthreads();

//	cudaThreadSynchronize();


	//去奇异


	//位移叠加


	//均值滤波


	//输出赋值

	//  *（disOutputCuda+bid）     =    displacement   ；

//	disOutputCuda[out_offset] = *displacement[threadIdx.x];


}



//去奇异        changed  by wong    2016/5/18

__global__  void   remove_singular_cuda(Complex*disOutputCuda, Complex*singularOutputCuda)   {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                         // 输出位移矩阵偏移值

	int   offrow = (blockIdx.x - 1)*blockDim.x + threadIdx.x;                   // 上一行输出位移矩阵偏移值

	int    bid   = blockIdx.x;                                                  //  block   id

	int    tid   = threadIdx.x;                                                //  thread  id             


	if (bid > 0 && (disOutputCuda[offset] > 12))  {

		singularOutputCuda[offset] = disOutputCuda[offrow];

	}

	else  {

		singularOutputCuda[offset] = disOutputCuda[offset];

	}


}


//位移叠加       changed   by  wong    2016/5/18

__global__   void   displace_add_cuda(Complex*singularOutputCuda, Complex*addOutputCuda)   {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                        // 输出位移矩阵偏移值

	int   offrow = (blockIdx.x - 1)*blockDim.x + threadIdx.x;                   // 上一行输出位移矩阵偏移值

	int    bid   = blockIdx.x;                                                  // block   id

	int    tid   = threadIdx.x;                                                 // thread  id             

	if (bid > 0)  {

		addOutputCuda[offset] = singularOutputCuda[offset] + singularOutputCuda[offrow];


	}

	else   {

		addOutputCuda[offset] = singularOutputCuda[offset];

	}



}



//数据扩展N-1列，增加冗余

__global__   void   extend_data_cuda(Complex*addOutputCuda, Complex*extendOutputCuda)   {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                        // 输出位移矩阵偏移值

	int    bid   = blockIdx.x;                                                  // block   id

	int    tid   = threadIdx.x;                                                 // thread  id   


	if (tid<N - 1)  {

		extendOutputCuda[offset] = 0;                                          //  extend  0

	}

	else
	{

		int   extoff = blockIdx.x *(blockDim.x - (N - 1)) + threadIdx.x - (N - 1);

		extendOutputCuda[offset] = addOutputCuda[extoff];

	}

}


//后累加平均   

__global__ void  smooth_filter_cuda(Complex*extendOutputCuda, Complex* smoothOutputCuda)   {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                        // 输出位移矩阵偏移值

	int   extbase = blockIdx.x*(blockDim.x + N - 1) + threadIdx.x;              // 基址

	int    bid = blockIdx.x;                                                    // block   id

	int    tid = threadIdx.x;                                                   // thread  id  


	Complex*  sum = 0;


	for (int i = extbase; i < extbase + N; i++)  {


		*sum += extendOutputCuda[i];


	}

	smoothOutputCuda[offset] = *sum / N;

}




__global__  void   timeField_filter_cuda(const Complex* smoothOutputCuda, const float* param,  const int  steps, Complex* timeFilterOutputCuda)    {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                        // 输出位移矩阵偏移值


	int    bid = blockIdx.x;                                                    // block   id

	int    tid = threadIdx.x;                                                   // thread  id     

	Complex  sum_temp = 0;

	float    coeff   = 0;

	for (int i = 0; i <= bid; i++)   {

		if ((bid - i) < steps)

			coeff = param[bid - i];

		else

			coeff = param[0];


		sum_temp += smoothOutputCuda[i*blockDim.x + threadIdx.x] * coeff;


	}

	timeFilterOutputCuda[offset] = sum_temp;

}






bool    CudaMain::isAvailable()  {

	int   count = 0;

	printf("Start to detecte devices.........\n");                   //  显示检测到的设备数

	cudaGetDeviceCount(&count);                                     //   检测计算能力大于等于1.0的设备数

	if (count == 0){

		fprintf(stderr, "There is no device.\n");

		return false;

	}


	printf("%d device/s detected.\n", count);                      //   显示检测到的设备数


	int i;

	for (i = 0; i < count; i++){                                  //  依次验证检测到的设备是否支持CUDA

		cudaDeviceProp prop;

		if (cudaGetDeviceProperties(&prop, i) == cudaSuccess) {  //  获得设备属性并验证是否正确

			if (prop.major >= 1)                                 //  验证主计算能力，即计算能力的第一位数是否大于1

			{
				printf("Device %d: %s supports CUDA %d.%d.\n", i + 1, prop.name, prop.major, prop.minor);//显示检测到的设备支持的CUDA版本
				break;


			}
		}
	}

	if (i == count) {                                         //   没有支持CUDA1.x的设备
		fprintf(stderr, "There is no device supporting CUDA 1.x.\n");
		return false;
	}

	cudaSetDevice(i);                                       //    设置设备为主叫线程的当前设备

	return true;

}





CudaMain::CudaMain()  {

	


	cpu_inputMat      = NULL;
	
	cpu_SplineOutMat  = NULL ;

	cpu_RadonMat      = NULL;

	cpu_WaveRate      = 0  ;

	mallocFlag        = false;

	cpu_config        = new   ConfigParam ;

	cpu_disMat        = NULL;



//	memset(cpu_config, 0, sizeof(cpu_config));




	inputMat         = NULL;

	zeroFilterMat    = NULL;

	frontFilterMat   = NULL;

	disOutput        = NULL;

	bandfilterParam  = NULL;

	lowfilterParam   = NULL;

	matchfilterParam = NULL;

	lowFrontMat      = NULL;

	lowBackMat       = NULL;

	singularOutputCuda = NULL;

	addOutputCuda      = NULL;

	extendOutputCuda   = NULL;


	radonIn          = NULL;

	radonOut         = NULL;


}




CudaMain :: ~CudaMain()  {

	freeMem(); 



}



void   CudaMain::inputConfigParam( ConfigParam*config) {



	cpu_config = config;


}



void  CudaMain::inputRfData(const EInput& in) {     //读入数据到cpu_inputMat中

	float* input = in.pDatas;


	for (int i = 0; i < cpu_inputMat->rows; i++)
	{
		for (int j = 0; j < cpu_inputMat->cols; j++)
		{
			*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*cpu_inputMat, i, j)))) = input[i * cpu_inputMat->cols + j];
		}
	}


}




void  CudaMain::getbandFilterParam(std::string paramFileName) {

	if (paramFileName.size() == 0)
	{
		exit(1);
	}

	std::fstream paramFile(paramFileName.c_str());

	if (!paramFile)
	{
		exit(1);
	}

	float tmp;

	std::string str;


	cpu_bandfilterParam.clear();

	while (!paramFile.eof())
	{
		paramFile >> tmp;
		cpu_bandfilterParam.push_back(tmp);
	}
	paramFile.close();



}


void   CudaMain::getlowFilterParam(std::string paramFileName) {

	if (paramFileName.size() == 0)
	{
		exit(1);
	}

	std::fstream paramFile(paramFileName.c_str());

	if (!paramFile)
	{
		exit(1);
	}

	float tmp;

	std::string str;


	cpu_lowfilterParam.clear();

	while (!paramFile.eof())
	{
		paramFile >> tmp;
		cpu_lowfilterParam.push_back(tmp);
	}
	paramFile.close();


}



void  CudaMain::getmatchFilterParam(std::string paramFileName) {

	if (paramFileName.size() == 0)
	{
		exit(1);
	}

	std::fstream paramFile(paramFileName.c_str());

	if (!paramFile)
	{
		exit(1);
	}

	float tmp;

	std::string str;


	cpu_matchfilterParam.clear();

	while (!paramFile.eof())
	{
		paramFile >> tmp;
		cpu_matchfilterParam.push_back(tmp);
	}
	paramFile.close();

}











void  CudaMain::mallocMem(void)  {

	mallocMats();

	mallocGPUMem();

	

}



void CudaMain::freeMem(void)  {

	freeMats();
   
	deleteGPUMem();
}





void   CudaMain::mallocGPUMem() {



	int  MatRows = cpu_config->shearFrameLineNum;

	int  MatCols = cpu_config->sampleNumPerLine ;

	int windowHW = cpu_config->windowHW;

	int maxLag   = cpu_config->maxLag;

	int step     = cpu_config->step;


	int interpnum  = cpu_config->fitline_pts;

	int iBPParaLen = 40;                                                      // bandpassfilter的长度；

	iBPParaLen     = (iBPParaLen > cpu_bandfilterParam.size()) ? iBPParaLen : cpu_bandfilterParam.size();


	int iLPParaLen = 40;                                                      // lowpassfilter的长度；

	iLPParaLen     = (iBPParaLen > cpu_lowfilterParam.size()) ? iBPParaLen : cpu_lowfilterParam.size();


	int iMHParaLen = 40;                                                      // matchfilter的长度；

	iMHParaLen    = (iBPParaLen > cpu_matchfilterParam.size()) ? iBPParaLen : cpu_matchfilterParam.size();






	if (MatRows == 0 || MatCols == 0)
	{

		printf("  row  and col  is zero! call InputConfigParas first!\n");
		return;

	}

	cudaError cudaStatus = cudaSetDevice(0);                                 // 0是titan显卡

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return;
	}



	int  multiWin = 2;                                                      //  大窗口对小窗口的倍数

	int cxorrLines = MatRows - 1;                                           //  位移矩阵相关扫描线数目        299

	int iOutRows = (MatCols - multiWin*windowHW) / step;                    //  位移矩阵计算需要匹配的段数     799 

	int extRows = iOutRows + N - 1;                                         //  扩展矩阵  799+100-1

	cudaMalloc(&disOutput, cxorrLines *iOutRows* sizeof(Complex));          //  位移矩阵GPU内存分配


//	cudaMalloc(&templateMatShare, cxorrLines *iOutRows* sizeof(templateMat));       //模板矩阵在GPU全局内存分配


//	cudaMalloc(&objectMatShare, cxorrLines *iOutRows* sizeof(objectMat));         //目标矩阵在GPU全局内存分配


//	cudaMalloc(&resultMatShare, cxorrLines *iOutRows* sizeof(resultMat));        //结果矩阵在GPU全局内存分配






	cudaMalloc(&singularOutputCuda, cxorrLines *iOutRows* sizeof(Complex)); // 去奇异矩阵在GPU内存分配


	cudaMalloc(&addOutputCuda, cxorrLines *iOutRows* sizeof(Complex));      // 位移叠加在GPU内存分配


	cudaMalloc(&extendOutputCuda, cxorrLines *extRows* sizeof(Complex));     // 扩展矩阵在GPU内存分配



	cudaMalloc(&inputMat, MatRows * MatCols * sizeof(Complex));             //   输入矩阵在GPU上对应的内存；


	cudaMalloc(&zeroFilterMat, MatRows * MatCols * sizeof(Complex));       //   带通零相位滤波输出在GPU内存分配；


	cudaMalloc(&frontFilterMat, MatRows * MatCols * sizeof(Complex));     //  带通零相位正滤波在GPU内存分配



	cudaMalloc(&lowBackMat, cxorrLines * iOutRows * sizeof(Complex));       //   低通零相位滤波输出在GPU内存分配；


	cudaMalloc(&lowFrontMat, cxorrLines * iOutRows * sizeof(Complex));     //  低通零相位正滤波在GPU内存分配

	


	cudaMalloc(&bandfilterParam, iBPParaLen * sizeof(float));                // iBPParaLen滤波器长度40


	cudaMalloc(&lowfilterParam, iLPParaLen * sizeof(float));                // iLPParaLen滤波器长度40


	cudaMalloc(&matchfilterParam, iMHParaLen * sizeof(float));              // iMHParaLen滤波器长度40



	int RadonInputCols      = 1961;                                     // 1961

	int RadonInputRows      = 4;                                       // 4

	cudaMalloc(&radonIn, sizeof(float) * RadonInputCols * RadonInputRows);                //拉东变换GPU输入

	cudaMalloc(&radonOut, sizeof(float) * RadonInputCols * (RadonInputCols - 1));        //拉东变换GPU输出  


	mallocFlag             = true;




}





void  CudaMain::deleteGPUMem()  {


	if (inputMat != NULL)
	{
		cudaFree(inputMat);

		inputMat = NULL;
	}

	
	if (zeroFilterMat != NULL)
	{
		cudaFree(zeroFilterMat);
		zeroFilterMat = NULL;
	}


	if (frontFilterMat != NULL)
	{
		cudaFree(frontFilterMat);
		frontFilterMat = NULL;
	}


	if (lowBackMat != NULL)
	{
		cudaFree(lowBackMat);
		lowBackMat = NULL;
	}


	if (lowFrontMat != NULL)
	{
		cudaFree(lowFrontMat);
		lowFrontMat = NULL;
	}




	if (disOutput != NULL)
	{
		cudaFree(disOutput);
		disOutput = NULL;
	}




	if (singularOutputCuda != NULL)
	{
		cudaFree(singularOutputCuda);

		singularOutputCuda = NULL;
	}


	if (addOutputCuda != NULL)
	{
		cudaFree(addOutputCuda);

		addOutputCuda = NULL;
	}


	if (extendOutputCuda != NULL)
	{
		cudaFree(extendOutputCuda);

		extendOutputCuda = NULL;
	}







	if (bandfilterParam != NULL)
	{
		cudaFree(bandfilterParam);
		bandfilterParam = NULL;
	}

	
	if (lowfilterParam != NULL)
	{
		cudaFree(lowfilterParam);
		lowfilterParam = NULL;
	}



	if (matchfilterParam != NULL)
	{
		cudaFree(matchfilterParam);
		matchfilterParam = NULL;
	}






	if (radonIn != NULL)
	{
		cudaFree(radonIn);
	}

	if (radonIn != NULL)
	{
		cudaFree(radonIn);
	}


	cudaDeviceReset();

	mallocFlag = false;


}




void  CudaMain::mallocMats() {


	cpu_inputMat    =   cvCreateMat(cpu_config->shearFrameLineNum, cpu_config->sampleNumPerLine, CV_32FC1);         //输入矩阵

	int  MatRows    = cpu_config->shearFrameLineNum;

	int  MatCols    = cpu_config->sampleNumPerLine;

	int windowHW    = cpu_config->windowHW;

	int maxLag      = cpu_config->maxLag;

	int step        = cpu_config->step;


	int  multiWin   = 2;                                                    //  大窗口对小窗口的倍数

	int cxorrLines  = MatRows - 1;                                         //   位移矩阵相关扫描线数目        299

	int iOutRows    = (MatCols - multiWin*windowHW) / step;               //    位移矩阵计算需要匹配的段数     799 

	cpu_disMat      = cvCreateMat(cxorrLines, iOutRows, CV_32FC1);       //     位移矩阵     


	cpu_SplineOutMat = cvCreateMat(1962, 4, CV_32FC1);                  //    SplineOutMat输出，便于画图，比较结果  

		
	cpu_RadonMat    = cvCreateMat(1962, 4, CV_32FC1);                  //     radon输出，比较计算结果  



	mallocFlag     = false; 
	

//	cpu_config     = (ConfigParam*)malloc(1 * sizeof(ConfigParam));     

	
//	memset(cpu_config, 0, sizeof(cpu_config));

}



void   CudaMain::freeMats() {

	if (cpu_inputMat != NULL)
	{
		cvReleaseMat(&cpu_inputMat);
		cpu_inputMat = NULL;
	}
	

	if (cpu_disMat != NULL)
	{
		cvReleaseMat(&cpu_disMat);
		cpu_disMat = NULL;
	}
	

	if (cpu_SplineOutMat != NULL)
	{
		cvReleaseMat(&cpu_SplineOutMat);
		cpu_SplineOutMat = NULL;
	}

	if (cpu_RadonMat != NULL)
	{
		cvReleaseMat(&cpu_RadonMat);
		cpu_RadonMat = NULL;
	}
	

	mallocFlag = NULL;


	free(cpu_config);


	cpu_config = NULL;

	  
}










CvMat*  CudaMain::bandpassFilt_cuda(CvMat* rawMat)  {


	Complex* h_MatData = (Complex*)rawMat->data.fl;

	cudaMemsetAsync(frontFilterMat, 0, sizeof(Complex)*rawMat->cols*rawMat->rows);

	cudaMemcpyAsync(zeroFilterMat, h_MatData, sizeof(Complex)*rawMat->cols*rawMat->rows, cudaMemcpyHostToDevice);    //拷贝CPU中RF数据到GPU

	int steps = cpu_bandfilterParam.size();

	cudaMemcpyAsync(bandfilterParam, &cpu_bandfilterParam[0], sizeof(float)*steps, cudaMemcpyHostToDevice);                  //拷贝CPU中抽头数据到GPU 





	dim3 blockID, threadID;

	blockID.x  = rawMat->rows;

	threadID.x = rawMat->cols;

	cudaThreadSynchronize();

	Bandpass_front_1 <<<blockID, threadID >> >(zeroFilterMat, rawMat->cols, bandfilterParam, steps, frontFilterMat);

	cudaThreadSynchronize();


	cudaMemcpy(zeroFilterMat, frontFilterMat, sizeof(Complex)*rawMat->cols*rawMat->rows, cudaMemcpyDeviceToDevice);


	Bandpass_back_1 << <blockID, threadID >> >(zeroFilterMat, rawMat->cols, bandfilterParam, steps, frontFilterMat);


	cudaThreadSynchronize();

	   
	cudaFree(bandfilterParam);

	cudaMemcpy(h_MatData, frontFilterMat, sizeof(Complex)*rawMat->cols*rawMat->rows, cudaMemcpyDeviceToHost);   //拷贝GPU处理后数据到	CPU

	cudaFree(zeroFilterMat);

	cudaFree(frontFilterMat);


	SaveDataFile("bpfilt.dat", rawMat);


	return rawMat;


}







void  CudaMain::zeroFilter_cuda(CvMat* rawMat, Complex*filterOutput) {







}





CvMat*  CudaMain::computeDisplacement_cuda(CvMat* filtOutMat, int  multiWin, int winSize, int stepSize){

//	CvMat*outputMat = 0;

	int     WinNum    = (filtOutMat->cols - multiWin*winSize) / stepSize;       //  一维位移矩阵

	Complex* hInput   = (Complex*)filtOutMat->data.fl;                         //   数据位置；

	//Complex*hOutput = (Complex*)outputMat->data.fl;                        //   数据位置；


	cudaMemcpy(inputMat, hInput, filtOutMat->cols*filtOutMat->rows*sizeof(Complex), cudaMemcpyHostToDevice);   //  CPU-GPU

	dim3 dBlock;

	dim3 dThread;

	dBlock.x = filtOutMat->rows - 1;                                 // 输出矩阵行数 ,块数        299

	dThread.x = WinNum;                                             // 输出矩阵列数 , 线程数      799


//	__device__   Complex*templateMatShare;                          //   模板内存在GPU分配         考虑局部分配                  


//	__device__   Complex*objectMatShare;                           //    目标内存在GPU分配         考虑局部分配


//	__device__   Complex*resultMatShare;                           //    匹配结果在GPU分配         考虑局部分配




	templateMat*templateMatShare;                                 //   模板内存在GPU分配 


	objectMat* objectMatShare;                                   //    目标内存在GPU分配 



	resultMat*resultMatShare;                                   //    匹配结果在GPU分配 



	Complex*      min;


	Complex*      max;

	int*          max_location;


	Complex*      displacement;








	cudaMalloc(&templateMatShare, dBlock.x*dThread.x* sizeof(templateMat));             //模板矩阵在GPU全局内存分配


	cudaMalloc(&objectMatShare,  dBlock.x*dThread.x* sizeof(objectMat));               //目标矩阵在GPU全局内存分配


	cudaMalloc(&resultMatShare,  dBlock.x*dThread.x* sizeof(resultMat));             //结果矩阵在GPU全局内存分配



	cudaMalloc(&min, dBlock.x*dThread.x* sizeof(Complex));                           // min在GPU全局内存分配


	cudaMalloc(&max, dBlock.x*dThread.x* sizeof(Complex));                          // max在GPU全局内存分配


	cudaMalloc(&max_location, dBlock.x*dThread.x* sizeof(int));                     // max_location在GPU全局内存分配


	cudaMalloc(&displacement, dBlock.x*dThread.x* sizeof(Complex));                // max_location在GPU全局内存分配



	
	//求位移矩阵  

	displacement_api_cuda << < dBlock, dThread >> >   (inputMat, filtOutMat->rows, filtOutMat->cols, multiWin, winSize, stepSize, templateMatShare, objectMatShare, resultMatShare, min, max, max_location, displacement);

	cudaThreadSynchronize();


	cudaFree(templateMatShare);

	cudaFree(objectMatShare);

	cudaFree(resultMatShare);

	cudaFree(min);

	cudaFree(max);

	cudaFree(max_location);








	//去奇异                                   

	remove_singular_cuda << <dBlock, dThread >> >   (displacement, singularOutputCuda);

	cudaThreadSynchronize();

	//位移叠加                 

	displace_add_cuda << <dBlock, dThread >> >  (singularOutputCuda, addOutputCuda);

	cudaThreadSynchronize();

	//前N-1列补0    

	int  ext_threads = dThread.x + N - 1;

	extend_data_cuda << < dBlock, ext_threads >> > (addOutputCuda, extendOutputCuda);

	cudaThreadSynchronize();

	cudaFree(addOutputCuda);

	//平滑滤波  

	smooth_filter_cuda << <dBlock, dThread >> >   (extendOutputCuda, disOutput);

	cudaThreadSynchronize();

	cudaFree(extendOutputCuda);

	//时域滤波，匹配滤波，50Hz增强     这里 param, iParaLen, steps 使用常量内存 

	int steps = cpu_matchfilterParam.size();

	cudaMemcpyAsync(matchfilterParam, &cpu_matchfilterParam[0], sizeof(float)*steps, cudaMemcpyHostToDevice);             //拷贝CPU中抽头数据到GPU 


	timeField_filter_cuda << <dBlock, dThread >> > (disOutput, matchfilterParam,  steps, singularOutputCuda);

	cudaThreadSynchronize();

	cudaFree(disOutput);

	//从GPU拷贝到CPU内存


	cudaMemcpy(hInput, singularOutputCuda, dBlock.x  * dThread.x*sizeof(Complex), cudaMemcpyDeviceToHost);   //  GPU-CPU


	cudaFree(singularOutputCuda);

	return    filtOutMat;



}






void   CudaMain::zeroDisplacement_cuda(CvMat* inputMat, int  multiWin, int winSize, int stepSize, Complex*disOutput){





}





CvMat*  CudaMain::lowpassFilt_cuda(CvMat* disMat)  {


	Complex* h_MatData = (Complex*)disMat->data.fl;

	cudaMemsetAsync(lowBackMat, 0, sizeof(Complex)*disMat->cols*disMat->rows);

	cudaMemcpyAsync(lowFrontMat, h_MatData, sizeof(Complex)*disMat->cols*disMat->rows, cudaMemcpyHostToDevice);           //拷贝CPU中RF数据到GPU

	int steps = cpu_lowfilterParam.size();

	cudaMemcpyAsync(lowfilterParam, &cpu_lowfilterParam[0], sizeof(float)*steps, cudaMemcpyHostToDevice);                  //拷贝CPU中抽头数据到GPU 





	dim3 blockID, threadID;

	blockID.x = disMat->rows;

	threadID.x = disMat->cols;

	cudaThreadSynchronize();

	Bandpass_front_1 << <blockID, threadID >> >(lowFrontMat, disMat->cols, lowfilterParam, steps, lowBackMat);

	cudaThreadSynchronize();


	cudaMemcpy(lowFrontMat, lowBackMat, sizeof(Complex)*disMat->cols*disMat->rows, cudaMemcpyDeviceToDevice);


	Bandpass_back_1 << <blockID, threadID >> >(lowFrontMat, disMat->cols, lowfilterParam, steps, lowBackMat);


	cudaThreadSynchronize();


	cudaFree(lowfilterParam);

	cudaMemcpy(h_MatData, lowBackMat, sizeof(Complex)*disMat->cols*disMat->rows, cudaMemcpyDeviceToHost);   //拷贝GPU处理后数据到	CPU

	cudaFree(lowFrontMat);

	cudaFree(lowBackMat);


	return disMat;


}









void  CudaMain::process(const EInput &input, EOutput& output) {

//	mallocMem();                                                                                 // 分配内存

//	inputRfData(input);                                                                          // 读取RF数据,给cpu_inputMat 

//	inputConfigParam(config);                                                                    // 配置参数 给cpu_config

// getFilterParam(config->bpfilt_file);                                                          // 读取滤波器带通参数，给cpu_filterParam




	   bandpassFilt_cuda(cpu_inputMat);                                                          //    带通滤波       


	   int  multiWin    = 2;

	   int winSize      = cpu_config->windowHW;

	   int  stepSize    = cpu_config->step;


	   
	   cpu_disMat       = computeDisplacement_cuda(cpu_inputMat,  multiWin,  winSize, stepSize);   //  位移计算       
	

	   lowpassFilt_cuda(cpu_disMat);                                                              //    低通滤波     




	   










}



