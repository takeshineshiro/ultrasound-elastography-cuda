//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "ErrorText.h"

// usage
//     CHAR msgText[256];
//     getLastErrorText(msgText,sizeof(msgText));
//     return error message
//     message buffer
//     buffer size
char *  GetLastErrorText(char *pBuf, DWORD bufSize)                     
{
	DWORD retSize;
	LPTSTR pTemp = NULL;

	if (bufSize < 16) 
	{
		if (bufSize > 0) 
		{
			pBuf[0]='\0';
		}
		return(pBuf);
	}
	retSize = FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER
		| FORMAT_MESSAGE_FROM_SYSTEM
		| FORMAT_MESSAGE_ARGUMENT_ARRAY,
		NULL,
		GetLastError(),
		LANG_NEUTRAL,
		(LPTSTR)&pTemp,
		0,
		NULL );
	if (!retSize || pTemp == NULL) 
	{
		pBuf[0]='\0';
	}
	else 
	{
		pTemp[strlen(pTemp) - 2] = '\0'; //remove cr and newline character
		sprintf(pBuf, "%0.*s (0x%x)", bufSize-16, pTemp, GetLastError());
		LocalFree((HLOCAL)pTemp);
	}
	return(pBuf);
}