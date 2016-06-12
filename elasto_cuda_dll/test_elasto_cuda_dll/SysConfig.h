//////////////////////////////////////////////////////////////////////////
//
//

#pragma once

#include <string>

const  char DefaultElastoConfFile[] = ".\\config.ini";

typedef struct ConfigParam
{
	//////////////////////////////////////////////////////////////////////////
	// 算法用的配置&&Start
	int     sampleFreqs;    //采样率;
	float   acousVel;       //声速, m/s;
	float   prf;

	float   threshold;      // 2012.9.5 奇异值滤波的阈值 王丛知
	int     windowHW;       // 2012.9.5 互相关窗口长度 王丛知
	int     maxLag;         // 2012.9.5 互相关计算的最大偏移（越大计算量越大） 王丛知
	int     step;           // 2012.9.5 互相关计算的步长（每次增加的新数据点数，因此窗口中大部分点都overlap） 王丛知

	int     fitline_pts;    // 最小二乘法计算直线拟合的点的数量

	std::string  lpfilt_file;
	std::string  bpfilt_file;
	std::string  matchfilt_file;

	int     box_x; // input data mat
	int     box_y;
	int     box_w;
	int     box_h;
	int     sb_x;  // strain data sub mat for ladong transform
	int     sb_y;
	int     sb_w;
	int     sb_h;

	int     radon_step;
	int     radon_num;
	std::string  calc_type;

	//  算法用的配置&&End
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	// 通讯协议相关&&Start
	int   lps;              // lines per second,每秒钟上传的线的数量
	int   sampleNumPerLine; // 每条线包含的样本点的数量
	int   elmSize;          // 每个样本点的字节长度
	int   shearFrameLineNum;// 剪切波数据帧包含的扫描线的数量

	// 通讯协议相关&&End
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	// UI相关&&Start
	int  mm_fr;

	int  m_rect_x;
	int  m_rect_y;
	int  m_rect_w;
	int  m_rect_h;

	int  e_rect_x;
	int  e_rect_y;
	int  e_rect_w;
	int  e_rect_h;

	int  s_rect_x;
	int  s_rect_y;
	int  s_rect_w;
	int  s_rect_h;

	// UI相关&&End
	//////////////////////////////////////////////////////////////////////////

	int nDyn;

} ConfigParam, SysConfig, *PSysConfig;

bool ReadSysConfig(const char *ini_file, ConfigParam &param);


