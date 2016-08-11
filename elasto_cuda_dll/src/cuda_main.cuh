#ifndef  _CUDA_MAIN_CUH

#define  _CUDA_MAIN_CUH

#include <opencv\cxcore.h>

#include <cuda_runtime.h>

#include <device_launch_parameters.h>


#include "ImageFunc.h"

#include <vector>


struct   EInput;

struct   EOutput;

struct   ConfigParam;


typedef float Complex;



#define    THREAD_NUM      799

#define    N               100



typedef   struct    templateData  {

	Complex  elem[64];

	Complex  atom[36];


}  templateData ;



typedef   struct    objectData  {

	Complex  elem_0[64];

	Complex  elem_1[64];

	Complex  elem_2[64];

	Complex  atom[8];


} objectData ;




typedef   struct     resultData  {

	Complex  elem[64];

	Complex  atom[37];



} resultData;





typedef  struct  templateMat{



	templateData  tempData;



}  templateMat;                                         




typedef  struct  objectMat{

 

 objectData    objData;


}  objectMat;                                            



typedef  struct resultMat{

	resultData   resData;


}  resultMat;                                         





typedef struct
{
	CvRect     rc;

	CvPoint    pt;

	float      xWidth;//横坐标间距，越小代表斜率越大

} RadonParam;





struct MyLessThan2
{
	bool operator()(const RadonParam &x, const RadonParam &y)
	{
		return x.xWidth > y.xWidth;
	}
};












class   CudaMain  {


private:


	CvMat*    cpu_inputMat;                      // 输入数据存放矩阵；

	CvMat*    cpu_disMat;                       //  位移矩阵    

	CvMat*    cpu_fitMat;                       


	CvMat*    cpu_SplineOutMat;                //  SplineOutMat输出，便于画图，比较结果         

	CvMat*    cpu_RadonMat;                    //  radon输出，比较计算结果       
    
	float     cpu_WaveRate;                    //  最终结果，波速；



	std::vector<float> cpu_lowfilterParam;     // 零相移滤波器抽头


	std::vector<float> cpu_bandfilterParam;    // 零相移滤波器抽头

	
	std::vector<float> cpu_matchfilterParam;   // 零相移滤波器抽头





	bool                mallocFlag;            // 内存分配完毕

	ConfigParam*        cpu_config;           // 配置参数




	Complex* inputMat;                       // 输入矩阵在GPU分配                         

	Complex* zeroFilterMat;                 // 零相位滤波输出在GPU内存分配


	Complex* frontFilterMat;               //  零相位正滤波在GPU内存分配\


	Complex* lowFrontMat;                 //   零相位滤波输出在GPU内存分配


	Complex* lowBackMat;                 //   零相位正滤波在GPU内存分配\




	Complex* disOutput;                    //  位移矩阵在GPU分配


	Complex*singularOutputCuda;            // 去奇异


	Complex*addOutputCuda;                 //位移叠加        


	Complex*extendOutputCuda;              // 前N-1列补0 


	float*   lowfilterParam;              // 滤波抽头在GPU内存分配  


	float*   bandfilterParam;             // 滤波抽头在GPU内存分配 


	float*   matchfilterParam;            // 滤波抽头在GPU内存分配 



	Complex*   fit_IN;                         // 应变输入
                           

	Complex*   fit_Out;                       //  应变输出



	float*    radonIn;                     //拉东变换在GPU内存分配

	float*    radonOut;                   //拉东变换在GPU内存分配




private:                                                                       // 私有函数


	void mallocGPUMem(void);                                                  // 分配GPU内存


	void deleteGPUMem(void);                                                 // 释放GPU内存


	void mallocMats(void);                                                  // 分配cpu内存


	void freeMats(void);                                                   // 释放CPU内存





 virtual	CvMat* bandpassFilt_cuda(CvMat* rawMat);                     




 virtual    CvMat* bandpassFilt_1024_cuda(CvMat* rawMat);                 



 virtual    void   zeroFilter_cuda(CvMat* rawMat,Complex*filterOutput);   


	                                                                                
 virtual  void computeDisplacement_cuda(CvMat* inputMat, int  multiWin, int winSize, int stepSize, CvMat*outputMat);


                                                                       
  virtual  void  zeroDisplacement_cuda(CvMat* inputMat, int  multiWin, int winSize, int stepSize, Complex*disOutput);





  virtual	CvMat* lowpassFilt_799_cuda(CvMat* disMat);                 



  virtual  void   strainCalculate_cuda(CvMat*dis,  CvMat* fitMat);     //计算应变值



  virtual  void   random_proess_cuda(CvMat*fitMat, ConfigParam*config, EOutput &output);


  virtual  void    RadonProcess2(CvPoint &s, CvPoint &e, ConfigParam*config,  const CvRect &rect, const CvMat &matStrain);


  virtual  void    RadonSum(const CvMat *pmatDisplacement, CvMat **ppmatRadan);  


  virtual  void  ImagePostProc(IplImage *pImg, const char *filename, const CvPoint &s, const CvPoint &e);




public:                                  



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