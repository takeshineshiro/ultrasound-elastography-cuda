#ifndef  _CUDA_MAIN_CUH

#define  _CUDA_MAIN_CUH

#include <opencv\cxcore.h>

#include <cuda_runtime.h>

#include <device_launch_parameters.h>

//#include "SysConfig.h"

#include <vector>


struct   EInput;

struct   EOutput;

struct   ConfigParam;


typedef float Complex;



#define    THREAD_NUM      799

#define    N               100



typedef  struct  templateMat{

	Complex  elem[100];


}  templateMat;                                         //changed   by   wong    2016/5/19




typedef  struct  objectMat{

	Complex  elem[200];


}  objectMat;                                            //changed   by   wong    2016/5/19



typedef  struct resultMat{

	Complex  elem[101];


}  resultMat;                                         //changed   by   wong    2016/5/19




class   CudaMain  {


private:

//cpu内存分配

	CvMat*    cpu_inputMat;                      // 输入数据存放矩阵；

	CvMat*    cpu_disMat;                       //  位移矩阵    

	CvMat*    cpu_SplineOutMat;                //  SplineOutMat输出，便于画图，比较结果         

	CvMat*    cpu_RadonMat;                    //  radon输出，比较计算结果       
    
	float     cpu_WaveRate;                    //  最终结果，波速；



	std::vector<float> cpu_lowfilterParam;     // 零相移滤波器抽头


	std::vector<float> cpu_bandfilterParam;    // 零相移滤波器抽头

	
	std::vector<float> cpu_matchfilterParam;   // 零相移滤波器抽头





	bool               mallocFlag;            // 内存分配完毕

	ConfigParam*        cpu_config;           // 配置参数





//gpu内存分配

	Complex* inputMat;                       // 输入矩阵在GPU分配                         

	Complex* zeroFilterMat;                 // 零相位滤波输出在GPU内存分配


	Complex* frontFilterMat;               //  零相位正滤波在GPU内存分配\


	Complex* lowFrontMat;                 //   零相位滤波输出在GPU内存分配


	Complex* lowBackMat;                 //   零相位正滤波在GPU内存分配\




	Complex* disOutput;                    //  位移矩阵在GPU分配


//	templateMat*templateMatShare;             //   模板内存在GPU分配         


//	objectMat* objectMatShare;             //    目标内存在GPU分配         


//	resultMat*resultMatShare;              //    匹配结果在GPU分配        



	Complex*singularOutputCuda;            // 去奇异


	Complex*addOutputCuda;                 //位移叠加        


	Complex*extendOutputCuda;              // 前N-1列补0 


	float*   lowfilterParam;              // 滤波抽头在GPU内存分配  


	float*   bandfilterParam;             // 滤波抽头在GPU内存分配 


	float*   matchfilterParam;            // 滤波抽头在GPU内存分配 


	float*    radonIn;                     //拉东变换在GPU内存分配

	float*    radonOut;                   //拉东变换在GPU内存分配




private:                                                                       // 私有函数


	void mallocGPUMem(void);                                                  // 分配GPU内存


	void deleteGPUMem(void);                                                 // 释放GPU内存


	void mallocMats(void);                                                  // 分配cpu内存


	void freeMats(void);                                                   // 释放CPU内存





 virtual	CvMat* bandpassFilt_cuda(CvMat* rawMat);                      //计算零相移滤波（带通），输出矩阵为了便于画图，保存，比较结果




 virtual    CvMat* bandpassFilt_1024_cuda(CvMat* rawMat);               //计算零相移滤波（带通）1024线程，输出矩阵为了便于画图，保存，比较结果       



 virtual    void   zeroFilter_cuda(CvMat* rawMat,Complex*filterOutput);   //计算零相移滤波（带通或低通）,输出保留在GPU中


	                                                                     // 计算一维位移矩阵，输出矩阵为了便于画图，保存，比较结果            
  virtual  CvMat*computeDisplacement_cuda(CvMat* inputMat, int  multiWin, int winSize, int stepSize);


                                                                        // 计算一维位移矩阵，输出保留在GPU中
  virtual  void  zeroDisplacement_cuda(CvMat* inputMat, int  multiWin, int winSize, int stepSize, Complex*disOutput);



  virtual	CvMat* lowpassFilt_cuda(CvMat* disMat);                      //计算零相移滤波（低通），输出矩阵为了便于画图，保存，比较结果




public:                                  //共有函数



	 CudaMain();


	~CudaMain();


	void  inputConfigParam( ConfigParam*config);              // 读取配置


	void  inputRfData(  const EInput& in);                   //  读取RF数据  


	void  getlowFilterParam(std::string paramFileName);      //   读取滤波器抽头


	void  getbandFilterParam(std::string paramFileName);   //   读取滤波器抽头


	void  getmatchFilterParam(std::string paramFileName);   //   读取滤波器抽头

	  
	void mallocMem(void);                                   //     分配CPU和GPU内存
	 
	void freeMem(void);                                     //     释放内存


	float getRate(void) const {                           //    获取速度
		 
		return cpu_WaveRate; 

	}  


	bool isAvailable();                                    //  是否存在GPU模块


	void process(const EInput &input, EOutput &output);    //  主要处理函数


















};





























#endif