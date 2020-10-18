// SignalsDlg.h: файл заголовка
//

#pragma once
#include <vector>
#include <math.h>
#include "Signals_Processing.h"
#include "Signals_helper.h"
#include "ChartViewer.h"
#include <fstream>
#include <iostream>
#include "rwfile.h"
#include <cstring>
using namespace std;
#define LOGS_PATH "Logs/"

// Диалоговое окно CSignalsDlg
class CSignalsDlg : public CDialogEx
{
	// Создание
public:
	CSignalsDlg(CWnd* pParent = nullptr);	// стандартный конструктор

// Данные диалогового окна
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_SIGNALS_DIALOG };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// поддержка DDX/DDV

// Реализация
protected:
	HICON m_hIcon;

	// Созданные функции схемы сообщений
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	
	/*Объект отрисовки*/
	CChartViewer viewer1;
	/*Объект отрисовки*/
	CChartViewer viewer2;
	/*Объект отрисовки*/
	CChartViewer viewer3;
	/*Объект класса обработки сигналов*/
	Signals_Processing sp;
	/*Преобразование вектора double к вектору DoubleArray*/
	DoubleArray vectorToArray(vector<double>& v);
	/*Функции отрисовки*/
	void ViewerDraw(vector<complex<double>>& data, int Xmax, CChartViewer& viewer_num);
	void ViewerDrawSpectrum(vector<double>& data, int Xmax, CChartViewer& viewer_num);
	void ViewerDraw(vector<double>& data, int Xmax, CChartViewer& viewer_num);
	void ViewerDraw(vector<vector<double>>& data, double Xmin, double Xmax, CChartViewer& viewer_num, string PathPic);
	/*Обработчик кнопки Cancel*/
	afx_msg void OnBnClickedCancel();
	/*Частота передачи данных*/
	double bitrate;
	/*Частота дискретизации*/
	double sampling;
	/*Кол-во бит в принятом сигнале*/
	int bits_size;
	/*Кол-во бит на задержку*/
	int delay_size;
	/*True задержка в битах*/
	int delay_lama;
	/*Функция обновления данных Signals_Processing*/
	void updateSP();
	/*Логическая переменная для отрисовки*/
	BOOL Signals_or_Spectrs;
	/*Обработчик кнопки генерации сигналов*/
	afx_msg void OnBnClickedButton1();
	/*Шум в дБ*/
	double noize_lvl;
	int _k;

	////////////////////////
	vector <complex<double>> ImSignal1;
	vector <complex<double>> ImSignal2;
	vector <complex<double>> ImSpectr1;
	vector <complex<double>> ImSpectr2;
	vector<vector<complex<double>>> SignalS; //вектор комплексных сигналов, заполняемых из файла
	vector <complex<double>> ImSignal2Rch;
	vector <double> Rrr;
	vector <double> ResearchRrr;
	vector <vector<float>> ResearchRrr2D;
	vector <pair<double,double>> veroiatnosti_fhss;
	vector <pair<double, double>> veroiatnosti_un;
	double f_dop;
	BOOL Dopler_On;
	BOOL scramble;
	////////////////////////// fast convolution
	afx_msg void OnBnClickedCheck1();
	afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedButton3();
	afx_msg void OnBnClickedButton4();
	afx_msg void OnBnClickedButton5();
	afx_msg void OnBnClickedButton6();
	afx_msg void OnBnClickedButton7();
	afx_msg void OnBnClickedButton8();
	afx_msg void OnBnClickedButton9();
	////////////////////////////////////////////
	template <class T>
	string convertToStrPng(T* var) {
		ostringstream ss;
		ss << *var;
		ss << ".png";
		return ss.str();
	}
	int delay_lama_r;
	bool vectorDoubleToFile(const vector <pair <double, double> >& points, const string& file)
	{
		if (points.empty())
			return false;

		vector <pair <double, double> >::const_iterator i;

		if (file != "")
		{
			ofstream out(file.c_str());
			if (out.fail())
			{
				out.close();
				return false;
			}
			for (i = points.begin(); i != points.end(); ++i)
			{
				out << i->first << "dB: " << (i->second)*100 <<"%"<< "\n";
			}

			out.close();
		}
		return true;
	}
	BOOL Signals_generator_type;
	int test_time_cr;
	void Signals_Gen(int bits_size,int delay_size, double noize_v)
	{
		updateSP();
		ResearchRrr.clear();
		ImSpectr1.clear();
		ImSpectr2.clear();
		if (Signals_generator_type) sp.Link16_Signals_Generator(ImSignal1, ImSignal2, bits_size, delay_size, scramble);
		else sp.Simple_Signals_Generator(ImSignal1, ImSignal2, bits_size, delay_size);
		ImSignal2Rch.clear();
		ImSignal2Rch = ImSignal2;
		if (Dopler_On)	sp.Dopler(ImSignal2, f_dop, sp.average_frequency);
		sp.addNoize(ImSignal2, noize_v);
		ImSpectr1.clear();
		ImSpectr2.clear();
		sp.FAST_FUR(ImSignal1, ImSpectr1, -1);
		sp.FAST_FUR(ImSignal2, ImSpectr2, -1);
		sp.spVertex(ImSpectr1);
		sp.spVertex(ImSpectr2);
	}
	/*Объект класса File_helper*/
	File_helper fh;
	/// <summary>
	/// VI COURSE
	/// </summary>
	Signals_helper MySignals;
	void TrueViewerDraw(vector<vector<double>>& data, double Xmin, double Xmax, CChartViewer& viewer_num, string PathPic, bool podpisi);
	afx_msg void OnBnClickedButton10();
	afx_msg void OnBnClickedButton11();
	afx_msg void OnBnClickedButton12();
};
