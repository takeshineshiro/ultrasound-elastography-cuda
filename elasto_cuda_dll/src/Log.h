#pragma once

#include <afxmt.h>

const char DefaultLogFilePath[] = "elasto.log";

class CLog
{
public:
	~CLog(void);

	static CLog *Instance();

	BOOL   Write(const char *str, int len);

private:
	CLog(const char *file = DefaultLogFilePath);

	static CLog *s_ptrInstance;

	CString      m_strFilePath;

	CSemaphore   m_semAccess;

	//CFile        m_fileLog;
    FILE        *m_pfLog;
};

