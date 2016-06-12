// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the ELASTO_CUDA_DLL_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// ELASTO_CUDA_DLL_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef ELASTO_CUDA_DLL_EXPORTS
#define ELASTO_CUDA_DLL_API __declspec(dllexport)
#else
#define ELASTO_CUDA_DLL_API __declspec(dllimport)
#endif

// This class is exported from the elasto_cuda_dll.dll
class ELASTO_CUDA_DLL_API Celasto_cuda_dll {
public:
	Celasto_cuda_dll(void);
	// TODO: add your methods here.
};

extern ELASTO_CUDA_DLL_API int nelasto_cuda_dll;

ELASTO_CUDA_DLL_API int fnelasto_cuda_dll(void);
