//////////////////////////////////////////////////////////////////////////
//

#include "stdafx.h"
#include "ImageFunc.h"
#include "opencv/highgui.h"

//////////////////////////////////////////////////////////////////////////
// 内部对psrc进行了转置操作
//////////////////////////////////////////////////////////////////////////
void MakeImage(const CvMat *psrc, const char *filename)
{
	CvMat *pmat = cvCreateMat(psrc->cols, psrc->rows, psrc->type);

	cvTranspose(psrc, pmat);
	IplImage *pimage = cvCreateImage(cvGetSize(pmat), IPL_DEPTH_32F, 3);
	cvCvtColor(pmat, pimage, CV_GRAY2BGR);

	cvSaveImage(filename, pimage);

	cvReleaseImage(&pimage);
	cvReleaseMat(&pmat);
}


/*************************************************
Function:      通过直方图变换进行图像增强，将图像灰度的域值拉伸到0-255
src1:               单通道灰度图像
dst1:              同样大小的单通道灰度图像
*************************************************/
int ImageStretchByHistogram(IplImage *src1, IplImage *dst1)
{
	//p[]存放图像各个灰度级的出现概率；
	//p1[]存放各个灰度级之前的概率和，用于直方图变换；
	//num[]存放图象各个灰度级出现的次数;
	assert(src1->width == dst1->width);
	double p[256], p1[256], num[256];

	//清空三个数组
	memset(p, 0, sizeof(p));
	memset(p1, 0, sizeof(p1));
	memset(num, 0, sizeof(num));
	int height = src1->height;
	int width = src1->width;
	long wMulh = height * width;

	//求存放图象各个灰度级出现的次数
	//statistics   
	for (int x = 0; x < src1->width; x++)
	{
		for (int y = 0; y < src1->height; y++){
			uchar v = ((uchar*)(src1->imageData + src1->widthStep*y))[x];
			num[v]++;
		}
	}

	//求存放图像各个灰度级的出现概率
	//calculate probability   
	for (int i = 0; i < 256; i++)
	{
		p[i] = num[i] / wMulh;
	}

	//求存放各个灰度级之前的概率和
	//p1[i]=sum(p[j]);  j<=i;   
	for (int i = 0; i < 256; i++)
	{
		for (int k = 0; k <= i; k++)
			p1[i] += p[k];
	}

	//直方图变换
	// histogram transformation   
	for (int x = 0; x < src1->width; x++)
	{
		for (int y = 0; y < src1->height; y++){
			uchar v = ((uchar*)(src1->imageData + src1->widthStep*y))[x];
			((uchar*)(dst1->imageData + dst1->widthStep*y))[x] = p1[v] * 255 + 0.5;
		}
	}
	return 0;
}

//图像亮度变换
int ImageAdjust(IplImage* src, IplImage* dst,
	double low, double high,   // X方向：low and high are the intensities of src
	double bottom, double top, // Y方向：mapped to bottom and top of dst
	double gamma)
{
	if (low < 0 && low>1 && high < 0 && high>1 &&
		bottom < 0 && bottom>1 && top < 0 && top>1 && low > high)
		return -1;
	double low2 = low * 255;
	double high2 = high * 255;
	double bottom2 = bottom * 255;
	double top2 = top * 255;
	double err_in = high2 - low2;
	double err_out = top2 - bottom2;

	int x, y;
	double val;

	// intensity transform
	for (y = 0; y < src->height; y++)
	{
		for (x = 0; x < src->width; x++)
		{
			val = ((uchar*)(src->imageData + src->widthStep*y))[x];
			val = pow((val - low2) / err_in, gamma) * err_out + bottom2;
			if (val > 255) val = 255; if (val < 0) val = 0; // Make sure src is in the range [low,high]
			((uchar*)(dst->imageData + dst->widthStep*y))[x] = (uchar)val;
		}
	}
	return 0;
}



/*************************************************
Function:
Description:     因为摄像头图像质量差，需要根据直方图进行图像增强，
将图像灰度的域值拉伸到0-255
Calls:
Called By:
Input:           单通道灰度图像
Output:          同样大小的单通道灰度图像
Return:
Others:           http://www.xiaozhou.net/ReadNews.asp?NewsID=771
DATE:               2007-1-5
*************************************************/
int ImageStretchByHistogram2(IplImage *src, IplImage *dst)
{
	//p[]存放图像各个灰度级的出现概率；
	//p1[]存放各个灰度级之前的概率和，用于直方图变换；
	//num[]存放图象各个灰度级出现的次数; 

	assert(src->width == dst->width);
	float p[256], p1[256], num[256];
	//清空三个数组
	memset(p, 0, sizeof(p));
	memset(p1, 0, sizeof(p1));
	memset(num, 0, sizeof(num));

	int height = src->height;
	int width = src->width;
	long wMulh = height * width;

	//求存放图象各个灰度级出现的次数
	// to do use openmp
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar v = ((uchar*)(src->imageData + src->widthStep*y))[x];
			num[v]++;
		}
	}

	//求存放图像各个灰度级的出现概率
	for (int i = 0; i<256; i++)
	{
		p[i] = num[i] / wMulh;
	}

	//求存放各个灰度级之前的概率和
	for (int i = 0; i<256; i++)
	{
		for (int k = 0; k <= i; k++)
			p1[i] += p[k];
	}

	//直方图变换
	// to do use openmp
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar v = ((uchar*)(src->imageData + src->widthStep*y))[x];
			((uchar*)(dst->imageData + dst->widthStep*y))[x] = p1[v] * 255 + 0.5;
		}
	}

	return 0;
}

int ChangeImgColor(IplImage *scr)
{
	CvScalar avgChannels = cvAvg(scr);
	double avgB = avgChannels.val[0];//获取第一通道平均值
	double avgG = avgChannels.val[1];//获取第二通道平均值
	double avgR = avgChannels.val[2];//获取第三通道平均值

	CvScalar idensity;
	int i = 0, j = 0;
	for (; i < scr->height; i++)
	{
		for (j = 0; j < scr->width; j++)
		{
			idensity = cvGet2D(scr, i, j);
			idensity.val[0] = idensity.val[0] - avgB + 19;//修改色素值
			idensity.val[1] = idensity.val[1] - avgG + 79;
			idensity.val[2] = idensity.val[2] - avgR + 168;
			cvSet2D(scr, i, j, idensity);
		}
	}

	return 0;
}
