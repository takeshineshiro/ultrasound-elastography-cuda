// TestTime.h: interface for the CTestTime class.
// 本文设计的类用于测试 程序的运行时间.
// 文本格式用ansi编译.不支持unicode编码.
//            
// history:
//     创建,  杨戈,  2012/3/26,
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TESTTIME_H__BED4EC0C_D99A_4741_AAF9_3A7256E418C2__INCLUDED_)
#define AFX_TESTTIME_H__BED4EC0C_D99A_4741_AAF9_3A7256E418C2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef  __AFX_H__
	#include <windows.h>

#endif

#include <time.h>
#include<string>


//获取的时间写入到mfc的调试窗口中，也就是调用了TRACE
class CTestTime
{
public:
	CTestTime();
	virtual ~CTestTime();

	virtual void run();
	virtual void stop(CString &text);//返回时间，字符表示
	virtual long stop();// 返回时间，ms

protected:
    virtual long getTimeTicks();

    long  getTimeTicksOnTimeCount();

    long  getTimeTicksOnThreadTime();

    long  getTimeTicksOnClock();

//private:
	int     m_startTickCount;//ms
	int     m_stopTickCount; //ms
};

#endif // !defined(AFX_TESTTIME_H__BED4EC0C_D99A_4741_AAF9_3A7256E418C2__INCLUDED_)
