// elasto_cuda_dll.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "elasto_cuda_dll.h"


// This is an example of an exported variable
ELASTO_CUDA_DLL_API int nelasto_cuda_dll=0;

// This is an example of an exported function.
ELASTO_CUDA_DLL_API int fnelasto_cuda_dll(void)
{
	return 42;
}

// This is the constructor of a class that has been exported.
// see elasto_cuda_dll.h for the class definition
Celasto_cuda_dll::Celasto_cuda_dll()
{
	return;
}
