#ifndef CDISPLACEMENT_H_H_H
#define CDISPLACEMENT_H_H_H
#pragma   once 

//#include "cv.h"
#include "CDataProcess.h"

class CDisplacement : public CDataProcess{

	friend class CStrain;
public:
	CDisplacement(CvMat*, int, int);
	void Do();
private:
	void displacementAlgorithm();
	void displacementAlgorithm2();

	void ComputeDispalcement();

	CvMat * DoCalcDisp(CvMat* inputMat, int windowHW, int maxLog, int step, float fs, float c);
	int winSize;
	int stepSize;
};
#endif //define CDISPLACEMENT_H_H_H