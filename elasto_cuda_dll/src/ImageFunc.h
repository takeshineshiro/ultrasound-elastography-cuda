//////////////////////////////////////////////////////////////////////////
//

#pragma  once

#include "opencv\cv.h"

void MakeImage(const CvMat *psrc, const char *filename);

int ImageStretchByHistogram(IplImage *src, IplImage *dst);

int ImageAdjust(IplImage* src, IplImage* dst,
	double low, double high,   // X方向：low and high are the intensities of src
	double bottom, double top, // Y方向：mapped to bottom and top of dst
	double gamma);

int ImageStretchByHistogram2(IplImage *src, IplImage *dst);

int ChangeImgColor(IplImage *scr);
