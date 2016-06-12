#include "stdafx.h"
#include "CDataProcess.h"

CDataProcess::CDataProcess()
{

}

CvMat* CDataProcess::doProcess(CDataProcess *dptr)
{
	dptr->Do();
	return outDataMat;
}

void CDataProcess::Do()
{
	
}

CvMat* CDataProcess::inDataMat = 0;
CvMat* CDataProcess::outDataMat = 0;

