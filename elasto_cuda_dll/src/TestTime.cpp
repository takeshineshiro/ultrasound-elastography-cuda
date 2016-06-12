// TestTime.cpp: implementation of the CTestTime class.
//
//////////////////////////////////////////////////////////////////////

//#include <afx.h>
#include "StdAfx.h"

#include "TestTime.h"

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <assert.h>
#include <mmsystem.h>

using namespace std;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


CTestTime::CTestTime(): m_startTickCount(0), m_stopTickCount(0)
{
}

CTestTime::~CTestTime()
{
}

long CTestTime::getTimeTicks()
{
    return getTimeTicksOnClock();
}

long CTestTime::getTimeTicksOnTimeCount()
{
    return ::GetTickCount();
}

long CTestTime::getTimeTicksOnClock()
{
    return ::clock();
}

//////////////////////////////////////////////////////////////////////////
//
//  GetTickCount,  据说精度达到18~20ms winbase.h
//  clock,         精度达到1ms, time.h
//  timeGetTime,   同clock, mmsystem.h,  winmm.lib 
//
//////////////////////////////////////////////////////////////////////////
long CTestTime::getTimeTicksOnThreadTime()
{
#define Li2Double(x) ((double)((x).HighPart) * 4.294967296E9 + (double)((x).LowPart))

	//time_t t;
	//clock_t t;
	//return (long)::GetTickCount();
	//return (long)::clock();
	//return (long)::timeGetTime();
	FILETIME createTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;
	LARGE_INTEGER   kernelTimeInt64;
	LARGE_INTEGER   userTimeInt64;

	BOOL ok = ::GetThreadTimes(::GetCurrentThread(), &createTime, &exitTime, &kernelTime, & userTime);

	if (ok)
	{
		kernelTimeInt64.LowPart = kernelTime.dwLowDateTime;
		kernelTimeInt64.HighPart = kernelTime.dwHighDateTime;

		userTimeInt64.LowPart  = userTime.dwLowDateTime;
		userTimeInt64.HighPart = userTime.dwHighDateTime;

		long  time_ms = (long) ((Li2Double(kernelTimeInt64) + Li2Double(userTimeInt64)) / 10000.0f);
		return time_ms;
	}
	else
	{
		return -1;
	}
}

void CTestTime::run()
{
	m_startTickCount = getTimeTicks();
	m_stopTickCount  = 0;
}

void CTestTime::stop(CString &text)
{
	m_stopTickCount = getTimeTicks();
	text.Format(_T("%dms"), m_stopTickCount - m_startTickCount);
}

long CTestTime::stop()
{
	m_stopTickCount = getTimeTicks();
	return (m_stopTickCount - m_startTickCount);
}
