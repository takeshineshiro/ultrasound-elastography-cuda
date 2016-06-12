#include "StdAfx.h"
#include "Log.h"


CLog *CLog::s_ptrInstance = NULL;

CLog::CLog(const char *file):m_strFilePath(file), m_pfLog(0)
{
	//BOOL ok = m_fileLog.Open(m_strFilePath, CFile::modeCreate | CFile::modeNoTruncate | CFile::modeWrite | CFile::typeText);
	//TRACE("%s Open is %s\n", ok ? "Ok" : "Fail");
	errno_t err = fopen_s(&m_pfLog, m_strFilePath, "a");
	TRACE("%s Open is %s\n", file, err == 0 ? "Ok" : "Fail");
}


CLog::~CLog(void)
{
	//m_fileLog.Close();
	if (m_pfLog)   fclose(m_pfLog);
}


CLog *CLog::Instance()
{
	if (s_ptrInstance == NULL)
	{
		s_ptrInstance = new CLog();
	}

	return s_ptrInstance;
}

BOOL CLog::Write(const char *str, int len)
{
	BOOL ok = FALSE;
	m_semAccess.Lock();
	if (m_pfLog)
	{
		fputs(str, m_pfLog);
		ok = TRUE;
	}
	m_semAccess.Unlock();
	return ok;
}
