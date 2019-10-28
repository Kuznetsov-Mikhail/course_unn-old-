
// Signals.h: главный файл заголовка для приложения PROJECT_NAME
//

#pragma once

#ifndef __AFXWIN_H__
	#error "включить pch.h до включения этого файла в PCH"
#endif

#include "resource.h"		// основные символы


// CSignalsApp:
// Сведения о реализации этого класса: Signals.cpp
//

class CSignalsApp : public CWinApp
{
public:
	CSignalsApp();

// Переопределение
public:
	virtual BOOL InitInstance();

// Реализация

	DECLARE_MESSAGE_MAP()
};

extern CSignalsApp theApp;
