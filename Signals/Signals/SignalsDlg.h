
// SignalsDlg.h: файл заголовка
//

#pragma once
#include <vector>
#include <math.h>
#include "Signals_Processing.h"
#include "ChartViewer.h"
#include <fstream>
using namespace std;

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
	CChartViewer viewer1;
	CChartViewer viewer2;
	CChartViewer viewer3;
	Signals_Processing sp;
	DoubleArray vectorToArray(vector<double>& v);
	void ViewerDraw(vector<complex<double>>& data, int Xmax, CChartViewer& viewer_num);
	void ViewerDrawSpectrum(vector<double>& data, int Xmax, CChartViewer& viewer_num);
	void ViewerDraw(vector<double>& data, int Xmax, CChartViewer& viewer_num);
	void ViewerDraw(vector<double>& data, double Xmin, double Xmax, CChartViewer& viewer_num, string PathPic);
	afx_msg void OnBnClickedCancel();
	double bitrate;
	double sampling;
	int bits_size;
	int delay_size;
	int delay_lama;
	void updateSP();
	BOOL Signals_or_Spectrs;
	afx_msg void OnBnClickedButton1();
	double noize_lvl;

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
	vector <pair<double,double>> veroiatnosti;
	double f_dop;
	double alfa;
	BOOL Dopler_On;
	afx_msg void OnBnClickedCheck1();
	afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedButton3();
	BOOL scramble;
	////////////////////////// fast convolution
	BOOL parallel;
	CProgressCtrl Prog_bar;
	CProgressCtrl ProgBarRes;
	afx_msg void OnBnClickedButton4();
	afx_msg void OnBnClickedButton5();
	afx_msg void OnBnClickedButton6();
	////////////////////////////////////////////
	/*template <class T>
	bool convertFromStr(string& str, T* var) {
		istringstream ss(str);
		return (ss >> *var);
	}*/

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
};
