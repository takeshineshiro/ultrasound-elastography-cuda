//////////////////////////////////////////////////////////////////////////
//
//
#include "stdafx.h"

#include "SysConfig.h"

bool ReadSysConfig(const char *ini_file, ConfigParam &param)
{
	//FILE *file = fopen(ini_file, "r");

	if (ini_file)
	{
		char text[100];

		//////////////////////////////////////////////////////////////////////////
		// 下面是和算法相关的配置
		param.acousVel = (float) GetPrivateProfileInt("Parameters", "acousVel", 1000, ini_file);
		param.sampleFreqs = GetPrivateProfileInt("Parameters", "sampleFreqs", 3000000, ini_file);

		param.windowHW = GetPrivateProfileInt("Parameters", "windowHW", 10, ini_file);
		param.step     = GetPrivateProfileInt("Parameters", "step", 3, ini_file);
		param.maxLag         = GetPrivateProfileInt("Parameters", "maxLag", 30, ini_file);

		GetPrivateProfileString("Parameters", "threshold", "", text, 99, ini_file);
		param.threshold = (float) strtod(text, NULL);

		GetPrivateProfileString("Parameters", "bpfilt_file", "", text, 99, ini_file);
		param.bpfilt_file = text;

		GetPrivateProfileString("Parameters", "lpfilt_file", "", text, 99, ini_file);
		param.lpfilt_file = text;

		GetPrivateProfileString("Parameters", "matchfilt_file", "", text, 99, ini_file);
		param.matchfilt_file = text;

		param.fitline_pts = GetPrivateProfileInt("Parameters", "fitline_pts", 2, ini_file);

		param.box_x = GetPrivateProfileInt("Parameters", "box_x", 0, ini_file);
		param.box_y = GetPrivateProfileInt("Parameters", "box_y", 0, ini_file);
		param.box_w = GetPrivateProfileInt("Parameters", "box_w", 8192, ini_file);
		param.box_h = GetPrivateProfileInt("Parameters", "box_h", 300, ini_file);
		param.sb_x  = GetPrivateProfileInt("Parameters", "sb_x",  -1, ini_file);
		param.sb_y  = GetPrivateProfileInt("Parameters", "sb_y",  -1, ini_file);
		param.sb_w  = GetPrivateProfileInt("Parameters", "sb_w",  -1, ini_file);
		param.sb_h  = GetPrivateProfileInt("Parameters", "sb_h",  -1, ini_file);
		param.radon_num = GetPrivateProfileInt("Parameters", "radon_num",  1, ini_file);
		param.radon_step = GetPrivateProfileInt("Parameters", "radon_step",  1, ini_file);
		GetPrivateProfileString("Parameters", "calc_type", "middle", text, 99, ini_file);
		param.calc_type = text;

		// 算法相关的配置结束
		//////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		// 通讯协议相关的配置
		param.lps   = GetPrivateProfileInt("Protocol", "lps", 100, ini_file);
		param.sampleNumPerLine = GetPrivateProfileInt("Protocol", "sampleNumPerLine", 4096, ini_file);
		param.elmSize = GetPrivateProfileInt("Protocol", "elmSize", 1, ini_file);
		param.shearFrameLineNum = GetPrivateProfileInt("Protocol", "shearFrameLineNum", 100, ini_file);
		//通讯协议相关配置&&结束
		//////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		// UI相关的配置
		param.mm_fr = GetPrivateProfileInt("UI", "mm_fr", 5, ini_file);

		param.m_rect_x =  GetPrivateProfileInt("UI", "m_rect_x", 10, ini_file);
		param.m_rect_y =  GetPrivateProfileInt("UI", "m_rect_y", 20, ini_file);
		param.m_rect_w =  GetPrivateProfileInt("UI", "m_rect_w", 100, ini_file);
		param.m_rect_h =  GetPrivateProfileInt("UI", "m_rect_h", 200, ini_file);
		param.e_rect_x =  GetPrivateProfileInt("UI", "e_rect_x", 10, ini_file);
		param.e_rect_y =  GetPrivateProfileInt("UI", "e_rect_y", 20, ini_file);
		param.e_rect_w =  GetPrivateProfileInt("UI", "e_rect_w", 100, ini_file);
		param.e_rect_h =  GetPrivateProfileInt("UI", "e_rect_h", 200, ini_file);
		param.s_rect_x =  GetPrivateProfileInt("UI", "s_rect_x", 10, ini_file);
		param.s_rect_y =  GetPrivateProfileInt("UI", "s_rect_y", 20, ini_file);
		param.s_rect_w =  GetPrivateProfileInt("UI", "s_rect_w", 100, ini_file);
		param.s_rect_h =  GetPrivateProfileInt("UI", "s_rect_h", 200, ini_file);

		// UI相关配置&&结束
		//////////////////////////////////////////////////////////////////////////

		param.nDyn = GetPrivateProfileInt("MMode", "Dyn", 40, ini_file);

		return true;
	}
	else
	{
		return false;
	}
}
