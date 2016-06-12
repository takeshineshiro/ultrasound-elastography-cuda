#ifndef     ELASTO_CUDA

#define     ELASTO_CUDA


#include <stdio.h>                //   引入C函数库-实际上本程序就是应该以C的方式编译，尽管其后缀为cpp类型
#include <stdlib.h>
#include <cuda_runtime.h>         //   引入CUDA运行时库头文件


#ifdef __cplusplus                //   指明函数的编译方式，以得到没有任何修饰的函数名


extern "C"
{
#endif

#ifdef CUDADLLTEST_EXPORTS

#define CUDADLLTEST_API __declspec(dllexport) //导出符号宏定义

#else

#define CUDADLLTEST_API __declspec(dllimport)

#endif


	CUDADLLTEST_API bool initCUDA();    //要导出的CUDA初始化函数



#ifdef __cplusplus
}
#endif


#endif