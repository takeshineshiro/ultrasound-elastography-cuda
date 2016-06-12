#ifndef CFILTER_H_H_H
#define CFILTER_H_H_H
#pragma   once 

#include "opencv/cv.h"
#include "CDataProcess.h"
#include <string>
#include <vector>

//#define  PI 3.1415926

class CFilter : public CDataProcess{

public:
	CFilter(std::string);
	void Do();

	void DoTimeFieldFilter(const CvMat *pSrc);

private:
	void filterAlgorithm();
	std::vector<float> param;
	int		steps;
};
#endif  //define CFILTER_H_H_H