#ifndef CDATA_H_H_H
#define CDATA_H_H_H
#pragma   once 

#include "opencv\cv.h"
#include <string>
class CData{

private:
	CvMat *dataMat;
public:
	CData(int rows, int cols);				//set the size of the data
	//void readData(std::string filename, float* line, int length);	//read from file
	//void readData();						//read from DAQ Card
	void  readData(float *input);
	CvMat *getData();
	CvMat *getSubData(int x, int y, int w, int h);
};

#endif	//define CDATA_H_H_H