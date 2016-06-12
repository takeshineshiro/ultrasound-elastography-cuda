#ifndef CSTRAIN_H_H_H
#define CSTRAIN_H_H_H
//#pragma   once 

#include "opencv\cv.h"
#include "CElasto.h"
#include "CDataProcess.h"
#include <iostream>

//////////////////////////////////////////////////////////////////////////
//存放 Radon 拉东变换 次数的开始，结束点坐标
//////////////////////////////////////////////////////////////////////////
typedef struct RadonResultPoint
{
	CvPoint startPoint;
	CvPoint endPoint;

} RadonResultPoint;



class CStrain : public CDataProcess{


public:

	CStrain();

	virtual  void Do();

//private:
	float strainAlgorithm(const EInput &input, EOutput &output);

	void  CalcStrain(const EInput &input, EOutput &output);

	// 一次拉东变换
	void  CalcStrain2(const EInput &input, EOutput &output);
	
	// 多次拉东变换
	void  CalcStrain3(const EInput &input, EOutput &output);

	void  RadonProcess(CvPoint &s, CvPoint &e, const CRect &rect, const CvMat &matStrain);

	void  RadonProcess2(CvPoint &s, CvPoint &e, const CRect &rect, const CvMat &matStrain);

	void  RadonProcess3(CvPoint &s,
	                    CvPoint &e,
	                    const CRect &sub_rc,
	                    const CvMat &matStrain,
						std::vector<RadonResultPoint>& vecStartEnd);


	//////////////////////////////////////////////////////////////////////////
	// 拉东变换
	// 输入：
	//   pmatDisplacement,   矩阵-应变；
	// 输出：
	//   ppmatRadon,   指针的地址， 函数创建一个矩阵保存拉东变换的结果，并把这个矩阵的地址
	//                 保存在ppmatRadon中。用户在使用必须释放它。
	//////////////////////////////////////////////////////////////////////////

	static void  RadonSum(const CvMat *pmatDisplacement, CvMat **ppmatRadan);


private:

	// 图像后处理

	void  ImagePostProc(IplImage *pImg, const char *filename, const CvPoint &s, const CvPoint &e);


	//////////////////////////////////////////////////////////////////////////
	// 计算应变值和应变图的灰度值
	// 输入：
	//    count， 拟合的点数
	//    pmat，  应变的矩阵
	//    pimg，  应变图
	//////////////////////////////////////////////////////////////////////////
	void  ComputeStrainValueAndImage(int count, CvMat *pmat, IplImage *pimg);

};



typedef struct 
{
	CRect   rc;

	CvPoint pt;

	float xWidth;//横坐标间距，越小代表斜率越大

} RadonParam;



struct MyLessThan
{
	bool operator() (RadonParam &lhs, RadonParam &rhs)
	{
		return (lhs.pt.x - lhs.pt.y) < (rhs.pt.x - lhs.pt.y);
	}
};



//creator wangxiaomeng 
//date    20160214
// /*按照降序排列*/

struct MyLessThan2
{
	bool operator()(const RadonParam &x, const RadonParam &y)
	{
		return x.xWidth > y.xWidth;
	}
};




#endif //define CSTRAIN_H_H_H