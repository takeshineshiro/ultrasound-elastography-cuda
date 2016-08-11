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
#include  <math.h>
#include <string.h>
#include <cstdio>
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include "ImageFunc.h"


#endif



__global__ void Bandpass_front_1(Complex* tInput, int iWidth, float* param, int iParaLen, Complex* tOutPut)
{
	
	float data_sum;

	float data_1;

	data_sum = 0.0;

	

	if (threadIdx.x <= iParaLen - 1)                                     
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
	else                                                                        
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






__global__ void Bandpass_front_1024(Complex* tInput, int iWidth, float* param, int iParaLen, Complex* tOutPut)     {

	 float data_sum;

	 float data_1;

	 data_sum = 0.0;

	 int   line_serial;

	 int    bid = blockIdx.x;


	   line_serial = bid /16;



	 int  line_mod = bid % 16;






	 if ((0 == line_mod))    {                                                                                    

		 if ((threadIdx.x <= iParaLen - 1))                                                       
		 {

			 for (int i = 0; i <= threadIdx.x; i++)
			 {


				 data_1 = *(tInput + line_serial * iWidth + threadIdx.x - i);                  //b(0)*x(n-0)+b(1)*x(n-1)+...+b(n)*x(0)   

				 data_sum += (data_1*param[i]);

			 }


			 data_1 = *(tInput + line_serial * iWidth);                                      // x(0)


			 for (int j = threadIdx.x + 1; j <= iParaLen - 1; j++)
			 {
				 data_sum += (data_1*param[j]);                                                    //b(n+1)*x(0)+...+b(nb-2)*x(0)
			 }

			 *(tOutPut + line_serial * iWidth + threadIdx.x) = data_sum;



		 }

		 else  if ((threadIdx.x > iParaLen - 1))   {                                               

			 
			 for (int i = 0; i <= iParaLen - 1; i++)
			 {

				 data_1 = *(tInput + line_serial *iWidth + threadIdx.x - i);                 //b(0)*x(n-0)+b(1)*x(n-1)+...+b(nb-2)*x(n-(nb-2))  

				 data_sum += (data_1*param[i]);

			 }

			 *(tOutPut + line_serial * iWidth + threadIdx.x) = data_sum;


		 }


	 }

	else                                                                                     
	{

		for (int i = 0; i <= iParaLen - 1; i++)
		{

			data_1    = *(tInput + blockIdx.x*blockDim.x+ threadIdx.x - i);   //b(0)*x(n-0)+b(1)*x(n-1)+...+b(nb-2)*x(n-(nb-2))  

			data_sum += (data_1*param[i]);

		}

		*(tOutPut + blockIdx.x*blockDim.x + threadIdx.x) = data_sum;

	}


 


}




__global__ void Bandpass_back_1(Complex* tInput, int iWidth, float* param, int iParaLen, Complex* tOutPut)
{
	

	float data_1;

	float data_sum;

	data_sum = 0.0;



	if (threadIdx.x <= iParaLen - 1)   {                               


		for (int i = 0; i <= threadIdx.x; i++)
		{


			data_1 = *(tInput + blockIdx.x*iWidth + iWidth - 1 - threadIdx.x + i);     

			data_sum += (data_1*param[i]);

		}


		data_1 = *(tInput + blockIdx.x * iWidth + iWidth - 1);                        //  x(N-1) 


		for (int j = threadIdx.x + 1; j <= iParaLen - 1; j++)
		{

			data_sum += (data_1*param[j]);                                            // b(n+1)*x(N-1)+...+b(nb-1)*x(N-1) 

		}


		*(tOutPut + blockIdx.x * iWidth + iWidth - 1 - threadIdx.x) = data_sum;                 // y(n)

	}

	else    {                                                                             


		for (int i = 0; i <= iParaLen - 1; i++)
		{
			data_1 = *(tInput + blockIdx.x*iWidth + iWidth - 1 - threadIdx.x + i);

			data_sum += (data_1*param[i]);                                         //  y(N-1-n) = b(0)*x(N-1-n+0) +b(1)*x(N-1-n+1)+...+b(nb-1)*x(N-1-n+nb-1)

		}

		*(tOutPut + blockIdx.x * iWidth + iWidth - 1 - threadIdx.x)   = data_sum;




	}



	__syncthreads();  

}






__global__ void Bandpass_back_1024(Complex* tInput, int iWidth, float* param, int iParaLen, Complex* tOutPut)   {


	float  data_1;

	float  data_sum;

	data_sum = 0.0;


	int   line_serial;

	int    bid = blockIdx.x;


	line_serial = bid / 16;



	int  line_mod = bid % 16;










	if ((0 == line_mod))    {                                                                   


		if (threadIdx.x <= iParaLen - 1)   {                                                      


			for (int i = 0; i <= threadIdx.x; i++)
			{

				 
				data_1 = *(tInput + line_serial*iWidth + iWidth - 1 - threadIdx.x + i);     

				data_sum += (data_1*param[i]);

			}


			data_1 = *(tInput + line_serial * iWidth + iWidth - 1);                          // x(N-1) 


			for (int j = threadIdx.x + 1; j <= iParaLen - 1; j++)
			{

				data_sum += (data_1*param[j]);                                                   // b(n+1)*x(N-1)+...+b(nb-1)*x(N-1) 

			}


			*(tOutPut + line_serial * iWidth + iWidth-1 - threadIdx.x) = data_sum;              // y(n)      

		}

		else  if (threadIdx.x  >iParaLen - 1)   {                                                 

			data_sum = 0;

			for (int i = 0; i <= iParaLen - 1; i++)
			{
				data_1 = *(tInput + line_serial*iWidth + iWidth - 1 - threadIdx.x + i);

				data_sum += (data_1*param[i]);                                               //  y(N-1-n) = b(0)*x(N-1-n+0) +b(1)*x(N-1-n+1)+...+b(nb-1)*x(N-1-n+nb-1)

			}

			*(tOutPut + line_serial * iWidth + iWidth - 1 - threadIdx.x) = data_sum;




		}



	}  

	else  {                                                                                     

		    data_sum = 0;
		  

		for (int i = 0; i <= iParaLen - 1; i++)
		{
			data_1 = *(tInput + line_serial*iWidth + iWidth - 1 - (threadIdx.x + line_mod*blockDim.x) + i);

			data_sum += (data_1*param[i]);                                               //  y(N-1-n) = b(0)*x(N-1-n+0) +b(1)*x(N-1-n+1)+...+b(nb-1)*x(N-1-n+nb-1)

		}

		*(tOutPut + line_serial * iWidth + iWidth - 1 - (threadIdx.x + line_mod*blockDim.x)) = data_sum;




	}
	
	
	
	
	}
















__device__      void   xcorr_cuda(const  Complex* templateMat_startID, const Complex* objectMat_startID, Complex*resultMat_startID)     {


         #pragma unroll

	for (int i = 0; i < 101; i++)   {


		Complex     sum_object = 0;

		Complex     frac_object = 0;


		Complex     pow_template = 0;


		Complex      pow_object = 0;


		Complex     result = 0;


	

		for (int j = 0; j < 100; j++)  {


			sum_object += *(objectMat_startID + i + j);


		}

		

		Complex   ave_object =   sum_object / 100;


		

		for (int j = 0; j < 100; j++)  {

			Complex    tmp = *(templateMat_startID + j) *  (*(objectMat_startID + i + j) - ave_object);


			frac_object += tmp;

		}


	

		for (int j = 0; j < 100; j++)  {


			pow_template += *(templateMat_startID + j) * *(templateMat_startID + j);

		}

	

		for (int j = 0; j < 100; j++)  {


			pow_object += *(objectMat_startID + i + j)* * (objectMat_startID + i + j);

		}

		

		result = sqrt(pow_template*pow_object);

		

		*(resultMat_startID + i) = frac_object / result;

	}


}



__device__      void   minMax_cuda(Complex*resultMat_startID, Complex* min_value, Complex*  max_value, int * max_location)   {



	for (int i = 0; i < 101; i++)  {

		if (*(resultMat_startID + i) >= *max_value)  {

			*max_location   = i;

			*max_value = *(resultMat_startID + i);



		}

	}




}


 


__device__    void    interp_cuda(Complex*resultMat_startID, int *  max_loc, Complex*max_value, int * multiWin, int * winSize, Complex*  displace)     {

	Complex*pre = (Complex*)resultMat_startID + *max_loc - 1;

	Complex*next = (Complex*)resultMat_startID + *max_loc + 1;


	*displace   = (*multiWin - 1) * *winSize / 2 - *max_loc - (*pre - *next) / (2 * (*pre - 2 * *max_value + *next));


}






__global__   void  displacement_api_cuda(Complex*disInputCuda, int rows, int cols, int  multiWin, int winSize, int  stepSize, templateMat*templateMatShare, objectMat* objectMatShare, resultMat*resultMatShare, Complex*min, Complex*max, int*max_location, Complex* displacement )      {


	int   out_offset = blockIdx.x *blockDim.x + threadIdx.x;                     

	int    bid       = blockIdx.x ;                                              
	
	int    tid       = threadIdx.x;                                        




	    Complex*templateMatID;                               //ID


	    Complex*objectMatID;                                //ID




	      templateMatID = (Complex*)(disInputCuda + blockIdx.x*cols + (multiWin - 1) * winSize / 2 + threadIdx.x * stepSize);




		  for (int i = 0; i < 100;i++)  {

			  if (i < 64)    {
			  
				  templateMatShare[out_offset].tempData.elem[i]= *(templateMatID + i);
			  
			  }
		  
			  else
				      

			    templateMatShare[out_offset].tempData.atom[i-64] = *(templateMatID + i);
		  
		  
		  }



		  objectMatID   = (Complex*)(disInputCuda + (blockIdx.x + 1)*cols + threadIdx.x * stepSize);


   
		  for (int i = 0; i < 200; i++)  {

			  if (i<64)
				  objectMatShare[out_offset].objData.elem_0[i]     = *(objectMatID + i);
			  else if (i<128)
				  objectMatShare[out_offset].objData.elem_1[i - 64] = *(objectMatID + i);

			  else if (i<192)
				  objectMatShare[out_offset].objData.elem_2[i - 128] = *(objectMatID + i);
			  else
				  objectMatShare[out_offset].objData.atom[i - 192]   = *(objectMatID + i);

			

		  }

       



		  for (int i = 0; i < 101; i++)  {
		   
			  if (i<64)
				  resultMatShare[out_offset].resData.elem[i]   = 0;
		   
			  else
				  resultMatShare[out_offset].resData.atom[i-64] = 0;
		  
		  }
		 





		  xcorr_cuda(templateMatShare[out_offset].tempData.elem, objectMatShare[out_offset].objData.elem_0, resultMatShare[out_offset].resData.elem);



		minMax_cuda(resultMatShare[out_offset].resData.elem, &min[out_offset], &max[out_offset], &max_location[out_offset]);


		interp_cuda(resultMatShare[out_offset].resData.elem, &max_location[out_offset], &max[out_offset], &multiWin, &winSize, &displacement[out_offset]);


	

}





__global__  void   remove_singular_cuda(Complex*disOutputCuda, Complex*singularOutputCuda)   {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                                           

	int    bid   = blockIdx.x;                                                                     //  block   id

	int    tid   = threadIdx.x;                                                                    //  thread  id    

	int    offrow = 0;

	if (bid  > 0 && bid < gridDim.x - 1 && tid < blockDim.x-1 )   {
	
		    offrow = (blockIdx.x - 1)*blockDim.x + threadIdx.x;
	
	} 




	if (bid > 0 && bid < gridDim.x - 1 && tid < blockDim.x - 1 && (abs(disOutputCuda[offset]) > 12))  {

		singularOutputCuda[offset] = disOutputCuda[offrow];

	}

	else  {

		singularOutputCuda[offset] = disOutputCuda[offset];

	}


}



__global__   void   displace_add_cuda(Complex*singularOutputCuda, Complex*addOutputCuda)   {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                               

	int    bid   = blockIdx.x;                                                          // block   id

	int    tid   = threadIdx.x;                                                         // thread  id   

	int   offrow = (bid >0 ) ? ( (blockIdx.x - 1)*blockDim.x + threadIdx.x)  :0 ;       

	int   nextoff =  (blockIdx.x + 1)*blockDim.x + threadIdx.x;


	Complex  sum = 0.0;




	if (bid > 0)  {

		     

		for (int i = 0; i < bid; i++)   {

			int  off = i*blockDim.x + threadIdx.x;


			sum = sum + singularOutputCuda[off];

		}


		addOutputCuda[offset] = singularOutputCuda[offset] + sum;


	}

	else   {

		addOutputCuda[offset] = singularOutputCuda[offset];

	}













}





__global__   void   extend_data_cuda(Complex*addOutputCuda, Complex*extendOutputCuda)   {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                        

	int    bid   = blockIdx.x;                                                  // block   id

	int    tid   = threadIdx.x;                                                 // thread  id   


	if (tid<N - 1)  {

		int   add_base = blockIdx.x *(blockDim.x - (N - 1));

		extendOutputCuda[offset] = addOutputCuda[add_base];                    //  extend  primites

	}

	else
	{

		int   extoff = blockIdx.x *(blockDim.x - (N - 1)) + threadIdx.x - (N - 1);

		extendOutputCuda[offset] = addOutputCuda[extoff];

	}

}


  

__global__ void  smooth_filter_cuda(Complex*extendOutputCuda, Complex* smoothOutputCuda)   {

	int   offset  = blockIdx.x *blockDim.x + threadIdx.x;

	int   extbase = blockIdx.x*(blockDim.x + N - 1) + threadIdx.x;              

	int    bid = blockIdx.x;                                                    // block   id

	int    tid = threadIdx.x;                                                   // thread  id  


	Complex   sum = 0;


	for (int i = extbase; i < extbase + N; i++)  {


		Complex  temp = *(extendOutputCuda + i);

		sum = sum + temp;


	}

	smoothOutputCuda[offset] = sum / N;

}




__global__  void   timeField_filter_cuda(const Complex* smoothOutputCuda, const float* param,  const int  steps, Complex* timeFilterOutputCuda)    {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                       


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

	printf("Start to detecte devices.........\n");                   

	cudaGetDeviceCount(&count);                                    

	if (count == 0){

		fprintf(stderr, "There is no device.\n");

		return false;

	}


	printf("%d device/s detected.\n", count);                      


	int i;

	for (i = 0; i < count; i++){                                 

		cudaDeviceProp prop;

		if (cudaGetDeviceProperties(&prop, i) == cudaSuccess) {  

			if (prop.major >= 1)                                

			{
				printf("Device %d: %s supports CUDA %d.%d.\n", i + 1, prop.name, prop.major, prop.minor);
				break;


			}
		}
	}

	if (i == count) {                                         
		fprintf(stderr, "There is no device supporting CUDA 1.x.\n");
		return false;
	}

	cudaSetDevice(i);                                       

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



void  CudaMain::inputRfData(const EInput& in) {    

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

	int iBPParaLen = 40;                                                     

	iBPParaLen     = (iBPParaLen > cpu_bandfilterParam.size()) ? iBPParaLen : cpu_bandfilterParam.size();


	int iLPParaLen = 40;                                                    

	iLPParaLen     = (iBPParaLen > cpu_lowfilterParam.size()) ? iBPParaLen : cpu_lowfilterParam.size();


	int iMHParaLen = 40;                                                     

	iMHParaLen    = (iBPParaLen > cpu_matchfilterParam.size()) ? iBPParaLen : cpu_matchfilterParam.size();






	if (MatRows == 0 || MatCols == 0)
	{

		printf("  row  and col  is zero! call InputConfigParas first!\n");
		return;

	}

	cudaError cudaStatus = cudaSetDevice(0);                                

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return;
	}



	int  multiWin = 2;                                                      //  大窗口对小窗口的倍数

	int cxorrLines = MatRows - 1;                                           //  位移矩阵相关扫描线数目        299

	int iOutRows = (MatCols - multiWin*windowHW) / step;                    //  位移矩阵计算需要匹配的段数     799 

	int extRows = iOutRows + N - 1;                                         //  扩展矩阵  799+100-1

	cudaMalloc(&disOutput, cxorrLines *iOutRows* sizeof(Complex));          //  位移矩阵GPU内存分配







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



	cudaMalloc(&fit_IN, cxorrLines *iOutRows* sizeof(Complex));          //  位移矩阵GPU内存分配


	int   points = 5;


	int   strain_col = iOutRows - points + 1;

	cudaMalloc(&fit_Out, cxorrLines *strain_col* sizeof(Complex));          //  位移矩阵GPU内存分配






	int RadonInputCols      = 1961;                                    

	int RadonInputRows      = 4;                                       

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



	int  fit_point  = 5;

	
	int  fit_cols = iOutRows - fit_point + 1;

	cpu_fitMat = cvCreateMat(cxorrLines, fit_cols, CV_32FC1);                 





	cpu_SplineOutMat = cvCreateMat(1962, 4, CV_32FC1);                  //    SplineOutMat输出，便于画图，比较结果  

		
	cpu_RadonMat    = cvCreateMat(1962, 4, CV_32FC1);                  //     radon输出，比较计算结果  



	mallocFlag     = false; 
	


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

	cudaMemcpyAsync(zeroFilterMat, h_MatData, sizeof(Complex)*rawMat->cols*rawMat->rows, cudaMemcpyHostToDevice);    

	int steps = cpu_bandfilterParam.size();

	cudaMemcpyAsync(bandfilterParam, &cpu_bandfilterParam[0], sizeof(float)*steps, cudaMemcpyHostToDevice);                 





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

	cudaMemcpy(h_MatData, frontFilterMat, sizeof(Complex)*rawMat->cols*rawMat->rows, cudaMemcpyDeviceToHost);   

	cudaFree(zeroFilterMat);

	cudaFree(frontFilterMat);


	return rawMat;


}





CvMat*  CudaMain::bandpassFilt_1024_cuda(CvMat* rawMat)  {


	Complex* h_MatData = (Complex*)rawMat->data.fl;

	cudaMemsetAsync(frontFilterMat, 0, sizeof(Complex)*rawMat->cols*rawMat->rows);

	cudaMemcpyAsync(zeroFilterMat, h_MatData, sizeof(Complex)*rawMat->cols*rawMat->rows, cudaMemcpyHostToDevice);    

	int steps = cpu_bandfilterParam.size();

	cudaMemcpyAsync(bandfilterParam, &cpu_bandfilterParam[0], sizeof(float)*steps, cudaMemcpyHostToDevice);                 


	dim3 blockID, threadID;

	blockID.x = rawMat->rows*16;                           





	threadID.x = rawMat->cols/16;                      



	cudaThreadSynchronize();

	Bandpass_front_1024 << <blockID, threadID >> >(zeroFilterMat, rawMat->cols, bandfilterParam, steps, frontFilterMat);

	cudaThreadSynchronize();



	cudaMemcpy(zeroFilterMat, frontFilterMat, sizeof(Complex)*rawMat->cols*rawMat->rows, cudaMemcpyDeviceToDevice);


	Bandpass_back_1024 << <blockID, threadID >> >(zeroFilterMat, rawMat->cols, bandfilterParam, steps, frontFilterMat);



	cudaThreadSynchronize();


	cudaFree(bandfilterParam);

	cudaMemcpy(h_MatData, frontFilterMat, sizeof(Complex)*rawMat->cols*rawMat->rows, cudaMemcpyDeviceToHost);   

	cudaFree(zeroFilterMat);

	cudaFree(frontFilterMat);



	return rawMat;
 






}










void  CudaMain::zeroFilter_cuda(CvMat* rawMat, Complex*filterOutput) {







}





void   CudaMain::computeDisplacement_cuda(CvMat* filtOutMat, int  multiWin, int winSize, int stepSize, CvMat*outputMat){



	int     WinNum    = (filtOutMat->cols - multiWin*winSize) / stepSize;     

	Complex* hInput   = (Complex*)filtOutMat->data.fl;                        

	Complex*hOutput  = (Complex*)outputMat->data.fl;                        


	cudaMemcpy(inputMat, hInput, filtOutMat->cols*filtOutMat->rows*sizeof(Complex), cudaMemcpyHostToDevice);   

	dim3 dBlock;

	dim3 dThread;

	dBlock.x = filtOutMat->rows - 1;                             


	dThread.x = WinNum;                                           


	templateMat*templateMatShare;                               


	objectMat* objectMatShare;                                  



	resultMat*resultMatShare;                                  



	Complex*      min;


	Complex*      max;

	int*          max_location;


	Complex*      displacement;



	cudaMalloc(&templateMatShare, dBlock.x*dThread.x* sizeof(templateMat));             


	cudaMalloc(&objectMatShare,  dBlock.x*dThread.x* sizeof(objectMat));               


	cudaMalloc(&resultMatShare,  dBlock.x*dThread.x* sizeof(resultMat));            



	cudaMalloc(&min, dBlock.x*dThread.x* sizeof(Complex));                          


	cudaMalloc(&max, dBlock.x*dThread.x* sizeof(Complex));                         


	cudaMalloc(&max_location, dBlock.x*dThread.x* sizeof(int));                     


	cudaMalloc(&displacement, dBlock.x*dThread.x* sizeof(Complex));               



	displacement_api_cuda << < dBlock, dThread >> >   (inputMat, filtOutMat->rows, filtOutMat->cols, multiWin, winSize, stepSize, templateMatShare, objectMatShare, resultMatShare, min, max, max_location, displacement);

	cudaThreadSynchronize();


	cudaFree(templateMatShare);

	cudaFree(objectMatShare);

	cudaFree(resultMatShare);

	cudaFree(min);

	cudaFree(max);

	cudaFree(max_location);

                               

	remove_singular_cuda << <dBlock, dThread >> >   (displacement, singularOutputCuda);

	cudaThreadSynchronize();

               

	displace_add_cuda << <dBlock, dThread >> >  (singularOutputCuda, addOutputCuda);

	cudaThreadSynchronize();
   

	int  ext_threads = dThread.x + N - 1;

	extend_data_cuda << < dBlock, ext_threads >> > (addOutputCuda, extendOutputCuda);

	cudaThreadSynchronize();

	cudaFree(addOutputCuda);
 

	smooth_filter_cuda << <dBlock, dThread >> >   (extendOutputCuda, disOutput);

	cudaThreadSynchronize();

	cudaFree(extendOutputCuda);


	int steps = cpu_matchfilterParam.size();

	cudaMemcpyAsync(matchfilterParam, &cpu_matchfilterParam[0], sizeof(float)*steps, cudaMemcpyHostToDevice);             


	timeField_filter_cuda << <dBlock, dThread >> > (disOutput, matchfilterParam,  steps, singularOutputCuda);

	cudaThreadSynchronize();

	cudaFree(disOutput);


	cudaMemcpy(hOutput, singularOutputCuda, dBlock.x  * dThread.x*sizeof(Complex), cudaMemcpyDeviceToHost);   


	cudaFree(singularOutputCuda);




}






void   CudaMain::zeroDisplacement_cuda(CvMat* inputMat, int  multiWin, int winSize, int stepSize, Complex*disOutput){





}




__global__ void     lowpass_front_799(Complex* tInput, int iWidth, float* param, int iParaLen, Complex* tOutPut)    {


	float data_sum;

	float data_1;

	data_sum = 0.0;


	if (threadIdx.x <= iParaLen - 1)                                   
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
	else                                                                         
	{
		
		for (int i = 0; i <= iParaLen - 1; i++)
		{

			data_1 = *(tInput + blockIdx.x*iWidth + threadIdx.x - i);   //b(0)*x(n-0)+b(1)*x(n-1)+...+b(nb-2)*x(n-(nb-2))  

			data_sum += (data_1*param[i]);

		}

		*(tOutPut + blockIdx.x * iWidth + threadIdx.x) = data_sum;

	}




}




__global__  void   lowpass_back_799(Complex* tInput, int iWidth, float* param, int iParaLen, Complex* tOutPut)   {


	

	float data_1;

	float data_sum;

	data_sum = 0.0;



	if (threadIdx.x <= iParaLen - 1)   {                              


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


		*(tOutPut + blockIdx.x * iWidth + iWidth - 1 -threadIdx.x) = data_sum;                 // y(n)

	}

	else    {                                                                             


		for (int i = 0; i <= iParaLen - 1; i++)
		{
			data_1 = *(tInput + blockIdx.x*iWidth + iWidth - 1 - threadIdx.x + i);

			data_sum += (data_1*param[i]);                                         //  y(N-1-n) = b(0)*x(N-1-n+0) +b(1)*x(N-1-n+1)+...+b(nb-1)*x(N-1-n+nb-1)

		}

		*(tOutPut + blockIdx.x * iWidth + iWidth - 1 - threadIdx.x) = data_sum;




	}



}









CvMat*  CudaMain::lowpassFilt_799_cuda(CvMat* disMat)  { 


	Complex* h_MatData = (Complex*)disMat->data.fl;

	cudaMemsetAsync(lowBackMat, 0, sizeof(Complex)*disMat->rows*disMat->cols);

	cudaMemcpyAsync(lowFrontMat, h_MatData, sizeof(Complex)*disMat->rows*disMat->cols, cudaMemcpyHostToDevice);          

	int steps = cpu_lowfilterParam.size();

	cudaMemcpyAsync(lowfilterParam, &cpu_lowfilterParam[0], sizeof(float)*steps, cudaMemcpyHostToDevice);                





	dim3 blockID, threadID;

	blockID.x = disMat->rows;

	threadID.x = disMat->cols;


	lowpass_front_799 << <blockID, threadID >> >(lowFrontMat, disMat->cols, lowfilterParam, steps, lowBackMat);

	cudaThreadSynchronize();

	
	cudaMemcpy(lowFrontMat, lowBackMat, sizeof(Complex)*disMat->rows*disMat->cols, cudaMemcpyDeviceToDevice);


	lowpass_back_799 << <blockID, threadID >> >(lowFrontMat, disMat->cols, lowfilterParam, steps, lowBackMat);


	cudaThreadSynchronize();


	cudaFree(lowfilterParam);

	cudaMemcpy(h_MatData, lowBackMat, sizeof(Complex)*disMat->cols*disMat->rows, cudaMemcpyDeviceToHost);   

	cudaFree(lowFrontMat);

	cudaFree(lowBackMat);


	return disMat;




}


__device__  void    fitLine_cv_func(Complex*xx_tmp, Complex*yy_tmp, Complex* result)   {


	Complex xmean = 0.0f;

	Complex ymean = 0.0f;

	for (int i = 0; i < 5; i++)
	{
		xmean += xx_tmp[i];

		ymean += yy_tmp[i];

	}


	xmean /= 5;

	ymean /= 5;


	Complex sumx2 = 0.0f;

	Complex sumxy = 0.0f;

	for (int i = 0; i < 5; i++)
	{

		sumx2 += (xx_tmp[i] - xmean) * (xx_tmp[i] - xmean);

		sumxy += (yy_tmp[i] - ymean) * (xx_tmp[i] - xmean);

	}


	*result = (Complex)(sumxy / sumx2);

}




__global__  void  fitLine_L2_cuda(Complex*strain_IN, Complex*xx_IN,  Complex*strainOut)   {

	int   offset = blockIdx.x *blockDim.x + threadIdx.x;                       


	int    bid = blockIdx.x;                                                    // block   id

	int    tid = threadIdx.x;                                                   // thread  id     


	int    act_off = blockIdx.x *(blockDim.x + 5 - 1) + threadIdx.x;  // input



	Complex   xx_tmp[5];

	Complex  yy_tmp[5];


	for (int i = 0; i < 5; i++)  {

		xx_tmp[i] = xx_IN[act_off + i];

		yy_tmp[i] = strain_IN[act_off + i];


	}


	fitLine_cv_func(xx_tmp, yy_tmp, &strainOut[offset]);


}







void  CudaMain::strainCalculate_cuda(CvMat*disMat,    CvMat* fitMat)  {


	Complex* h_MatData = (Complex*)disMat->data.fl;


	Complex* out_MatData = (Complex*)fitMat->data.fl;


	cudaMemcpyAsync(fit_IN, h_MatData, sizeof(Complex)*disMat->rows*disMat->cols, cudaMemcpyHostToDevice);          


	int   fit_points = 5;


	dim3 blockID, threadID;

	blockID.x = disMat->rows;

	threadID.x = disMat->cols - fit_points +1;




	int   xx_rows = disMat->rows;

	int  xx_cols  = disMat->cols;


	CvMat*    xx_mat = cvCreateMat(xx_rows, xx_cols, CV_32FC1);


	   
	for (int i = 0; i < xx_rows; i++)   {

		for (int j = 0; j < xx_cols; j++)  {
		
		
			*(static_cast<float*>(static_cast<void*>(CV_MAT_ELEM_PTR(*xx_mat, i, j)))) = j*fit_points;
				
		}
	
	}


	Complex* xx_IN  ;

    cudaMalloc(&xx_IN, sizeof(Complex)*xx_rows*xx_cols);           


	Complex* xx_Data = (Complex*)xx_mat->data.fl;


	cudaMemcpyAsync(xx_IN, xx_Data, sizeof(Complex)*xx_mat->rows*xx_mat->cols, cudaMemcpyHostToDevice);



	fitLine_L2_cuda << <blockID, threadID >> >  (fit_IN, xx_IN,  fit_Out);


	cudaMemcpy(out_MatData, fit_Out, sizeof(Complex)*fitMat->rows*fitMat->cols, cudaMemcpyDeviceToHost);   




	cudaFree(fit_IN);

	cudaFree(xx_IN);

	cudaFree(fit_Out);







}




void  CudaMain::ImagePostProc(IplImage *strImage, const char *filename, const CvPoint &start, const CvPoint &end)
{

	const char * gray_file = "strain_gpu_gray.bmp";


	{
		IplImage *pimgStrain = cvCreateImage(cvGetSize(strImage), strImage->depth, 3);

		cvCvtColor(strImage, pimgStrain, CV_GRAY2BGR);

		cvSaveImage(gray_file, pimgStrain);

		cvReleaseImage(&pimgStrain);

	}

	{

		IplImage *pImage = cvLoadImage(gray_file, 0);

		IplImage *pimgStrain = cvCreateImage(cvGetSize(pImage), pImage->depth, 3);

		pimgStrain = cvCreateImage(cvGetSize(pImage), pImage->depth, 3);


	
		ImageAdjust(pImage, pImage, 0, 0.5, 0, 0.5, 0.6);

		

		cvNot(pImage, pImage);

		
		cvCvtColor(pImage, pimgStrain, CV_GRAY2BGR);


		ChangeImgColor(pimgStrain);


		cvLine(pimgStrain, start, end, CV_RGB(255, 0, 0), 2, CV_AA, 0);  


		cvSaveImage(filename, pimgStrain);


		
		cvReleaseImage(&pImage);

		cvReleaseImage(&pimgStrain);

	}

}






//////////////////////////////////////////////////////////////////////////
// 拉东变换
// pmatDisplacement,   rows: disp;  cols: time-extent( lines)
//     列,表示一条线, 也就是时间 轴
//     行,表示应变的值
//////////////////////////////////////////////////////////////////////////

void   CudaMain::RadonSum(const CvMat *pmatDisplacement, CvMat **ppmatRodan) {


	int xstart          = 0;

	int xend            = pmatDisplacement->rows;                     

	int t               = pmatDisplacement->cols;                    

	CvMat *pmatRodan    = cvCreateMat(t - 1, t, pmatDisplacement->type);

	cvZero(pmatRodan);

	int tstart          = 0;

	int tend            = 0;

	int dx              = 0;

	float dt            = 0.0f;

	float c             = 0.0f;


	for (tstart = 0; tstart < t - 1; tstart++)
	{

		for (tend = tstart + 1; tend < t; tend++)
		{

			c = (float)(xend - xstart) / (tend - tstart);                     //k

			for (dx = xstart; dx < xend; dx++)
			{

				dt = tstart + (dx - xstart) / c;                             //

				CV_MAT_ELEM(*pmatRodan, float, tstart, tend) = CV_MAT_ELEM(*pmatRodan, float, tstart, tend)
					+ CV_MAT_ELEM(*pmatDisplacement, float, dx, (int)dt);

			}
		}
	}


	*ppmatRodan = pmatRodan;






}






void  CudaMain::RadonProcess2(CvPoint &s, CvPoint &e, ConfigParam*config, const CvRect &sub_rc, const CvMat &matStrain)
{

	int  radon_num = config->radon_num;                    
	 

	int  radon_step = config->radon_step;                  



	int  intpl_multiple = 1;                                



	std::vector<RadonParam> array_params;



	for (int i = 0; i < radon_num; i++)                     
	{


		RadonParam param;

		param.rc.x = sub_rc.x;

		param.rc.y = sub_rc.y + i*radon_step;

		param.rc.width = sub_rc.width;

		param.rc.height = sub_rc.height;


		CvMat *pmatSub = cvCreateMatHeader(param.rc.height-1, param.rc.width-1, matStrain.type);


		cvGetSubRect(&matStrain, pmatSub, cvRect(param.rc.x, param.rc.y, param.rc.width-1, param.rc.height-1));


		CvMat *pmatRadon = 0;


		CvMat *pmatMultiple = cvCreateMat(pmatSub->rows, pmatSub->cols * intpl_multiple, pmatSub->type);


		cvResize(pmatSub, pmatMultiple);


		RadonSum(pmatMultiple, &pmatRadon);


		double  min_val;


		double  max_val;


		CvPoint min_loc;


		CvPoint max_loc;


		cvMinMaxLoc(pmatRadon, &min_val, &max_val, &min_loc, &max_loc);


		param.pt = max_loc;


		param.xWidth = param.pt.y - param.pt.x;//add by wxm


		array_params.push_back(param);


		cvReleaseMat(&pmatRadon);


		cvReleaseMat(&pmatMultiple);


		cvReleaseMatHeader(&pmatSub);


	}


	std::sort(array_params.begin(), array_params.end(), MyLessThan2());



	if (config->calc_type.compare("middle") == 0)
	{

		int size = array_params.size();


		s.x = array_params[size / 2].pt.y / intpl_multiple;


		s.y = array_params[size / 2].rc.y;


		e.x = array_params[size / 2].pt.x / intpl_multiple;


		e.y = array_params[size / 2].rc.y + array_params[size / 2].rc.height-1;

	}

	else if (config->calc_type.compare("max") == 0)
	{

		int size = array_params.size();

		s.x = array_params[0].pt.y / intpl_multiple;

		s.y = array_params[0].rc.y;


		e.x = array_params[0].pt.x / intpl_multiple;

		e.y = array_params[0].rc.y + array_params[0].rc.height-1;

	}

	else if (config->calc_type.compare("min") == 0)

	{

		int size = array_params.size();

		s.x = array_params[size - 1].pt.y / intpl_multiple;

		s.y = array_params[size - 1].rc.y;


		e.x = array_params[size - 1].pt.x / intpl_multiple;

		e.y = array_params[size - 1].rc.y + array_params[size - 1].rc.height-1;

	}

	else
	{
		//

	}










}








void    CudaMain::random_proess_cuda(CvMat*fitMat, ConfigParam*config, EOutput &output)  {

	
	

		int    win_size          = config->windowHW;                                             


		double overlap           = (config->windowHW - config->step) / (float)config->windowHW;  

		double sound_velocity    = config->acousVel;                                            


		double sample_frq        = config->sampleFreqs;                                                             

		double prf               = 1 / 300e-6;                                                  


		int    dep_start         = (config->sb_x < 0) ? 0 : config->sb_x;

		int    dep_size          = (config->sb_w < 0) ? fitMat->width : config->sb_w;

		int    dep_end           = dep_start + dep_size - 1;

		int    t_start           = (config->sb_y < 0) ? 0 : config->sb_y;

		int    t_size            = (config->sb_h < 0) ? fitMat->rows : config->sb_h;

		int    t_end             = t_start + t_size - 1;


		CvMat *pmatStrainTran    = cvCreateMat(fitMat->cols, fitMat->rows, fitMat->type);     


		cvTranspose(fitMat, pmatStrainTran);


		CvPoint                   start;

		CvPoint                    end;
		
		CvRect                     rect;


		rect.x                    = t_start;

		rect.y                    = dep_start;

		rect.width                = t_size;

		rect.height               = dep_size;



		
#if 1
		RadonProcess2(start, end, config ,rect, *pmatStrainTran);
#endif






		double v                  = ((end.y - start.y) * win_size * (1 - overlap) * sound_velocity / sample_frq / 2)
			/ ((end.x - start.x) / prf);



		double e                  = v * v * 3;


		output.v                  = (float)v;


		output.e                  = (float)e;

		cvReleaseMat(&pmatStrainTran);


		
    

}












void  CudaMain::process(const EInput &input, EOutput& output) {




	   bandpassFilt_1024_cuda(cpu_inputMat);                                                      


	   int  multiWin    = 2;

	   int winSize      = cpu_config->windowHW;

	   int  stepSize    = cpu_config->step;

   
	    computeDisplacement_cuda(cpu_inputMat, multiWin, winSize, stepSize, cpu_disMat);            
	

		lowpassFilt_799_cuda (cpu_disMat);                                                    


		strainCalculate_cuda(cpu_disMat,   cpu_fitMat);                                     


		random_proess_cuda(cpu_fitMat, cpu_config, output);                                




		int   ss = 0;




}



