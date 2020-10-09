
// SignalsDlg.cpp: файл реализации
//
#pragma once
#include "pch.h"
#include "framework.h"
#include "Signals.h"
#include "SignalsDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// Диалоговое окно CSignalsDlg



CSignalsDlg::CSignalsDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_SIGNALS_DIALOG, pParent)
	, bitrate(5e6)
	, sampling(35e7)
	, bits_size(1000)
	, delay_size(100)
	, delay_lama(0)
	, Signals_or_Spectrs(FALSE)
	, noize_lvl(100)
	, f_dop(3e4)
	, Dopler_On(TRUE)
	, scramble(TRUE)
	, Signals_generator_type(TRUE)
	, test_time_cr(0)
	, _k(32)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CSignalsDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_SIGNAL1, viewer1);
	DDX_Control(pDX, IDC_SIGNAL2, viewer2);
	DDX_Control(pDX, IDC_SIGNAL3, viewer3);
	DDX_Text(pDX, IDC_EDIT1, bitrate);
	DDX_Text(pDX, IDC_EDIT2, sampling);
	DDX_Text(pDX, IDC_EDIT3, bits_size);
	DDX_Text(pDX, IDC_EDIT4, delay_size);
	DDX_Text(pDX, IDC_EDIT5, delay_lama);
	DDX_Check(pDX, IDC_CHECK1, Signals_or_Spectrs);
	DDX_Text(pDX, IDC_EDIT6, noize_lvl);
	DDX_Text(pDX, IDC_EDIT7, f_dop);
	DDX_Check(pDX, IDC_CHECK2, Dopler_On);
	DDX_Check(pDX, IDC_CHECK3, scramble);
	DDX_Check(pDX, IDC_CHECK5, Signals_generator_type);
	DDX_Text(pDX, IDC_EDIT9, test_time_cr);
	DDX_Text(pDX, IDC_EDIT8, _k);
}

BEGIN_MESSAGE_MAP(CSignalsDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDCANCEL, &CSignalsDlg::OnBnClickedCancel)
	ON_BN_CLICKED(IDC_BUTTON1, &CSignalsDlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_CHECK1, &CSignalsDlg::OnBnClickedCheck1)
	ON_BN_CLICKED(IDC_BUTTON2, &CSignalsDlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_BUTTON3, &CSignalsDlg::OnBnClickedButton3)
	ON_BN_CLICKED(IDC_BUTTON4, &CSignalsDlg::OnBnClickedButton4)
	ON_BN_CLICKED(IDC_BUTTON5, &CSignalsDlg::OnBnClickedButton5)
	ON_BN_CLICKED(IDC_BUTTON6, &CSignalsDlg::OnBnClickedButton6)
	ON_BN_CLICKED(IDC_BUTTON7, &CSignalsDlg::OnBnClickedButton7)
	ON_BN_CLICKED(IDC_BUTTON8, &CSignalsDlg::OnBnClickedButton8)
	ON_BN_CLICKED(IDC_BUTTON9, &CSignalsDlg::OnBnClickedButton9)
	ON_BN_CLICKED(IDC_BUTTON10, &CSignalsDlg::OnBnClickedButton10)
	ON_BN_CLICKED(IDC_BUTTON11, &CSignalsDlg::OnBnClickedButton11)
END_MESSAGE_MAP()


// Обработчики сообщений CSignalsDlg

BOOL CSignalsDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Задает значок для этого диалогового окна.  Среда делает это автоматически,
	//  если главное окно приложения не является диалоговым
	SetIcon(m_hIcon, TRUE);			// Крупный значок
	SetIcon(m_hIcon, FALSE);		// Мелкий значок

	ViewerDraw(ImSignal1, ImSignal1.size(), viewer1);
	ViewerDraw(ImSignal2, ImSignal2.size(), viewer2);
	ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
	// TODO: добавьте дополнительную инициализацию
	return TRUE;  // возврат значения TRUE, если фокус не передан элементу управления
}

// При добавлении кнопки свертывания в диалоговое окно нужно воспользоваться приведенным ниже кодом,
//  чтобы нарисовать значок.  Для приложений MFC, использующих модель документов или представлений,
//  это автоматически выполняется рабочей областью.

void CSignalsDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // контекст устройства для рисования

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Выравнивание значка по центру клиентского прямоугольника
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Нарисуйте значок
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// Система вызывает эту функцию для получения отображения курсора при перемещении
//  свернутого окна.
HCURSOR CSignalsDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

DoubleArray CSignalsDlg::vectorToArray(std::vector<double>& v)
{
	return (v.size() == 0) ? DoubleArray() : DoubleArray(v.data(), (int)v.size());
}
void CSignalsDlg::ViewerDraw(vector<complex<double>>& data, int Xmax, CChartViewer& viewer_num)
{
	// In this example, we simply use random data for the 3 data series.
	vector<double> dataReal, dataImag;
	for (int i = 0; i < data.size(); i++)
	{
		dataReal.push_back(data[i].real());
		dataImag.push_back(data[i].imag());
	}
	DoubleArray Arr_dataReal = vectorToArray(dataReal);
	DoubleArray Arr_dataImag = vectorToArray(dataImag);
	vector<double>Datatime;
	for (int i = 0; i < Xmax; i++)Datatime.push_back(i);
	DoubleArray timeStamps = vectorToArray(Datatime);

	// Create a XYChart object of size 600 x 400 pixels
	XYChart* c = new XYChart(850, 240);

	// Add a title box using grey (0x555555) 20pt Arial font
	//c->addTitle("", "arial.ttf", 20, 0x555555);

	// Set the plotarea at (70, 70) and of size 500 x 300 pixels, with transparent background and
	// border and light grey (0xcccccc) horizontal grid lines
	c->setPlotArea(70, 50, 700, 120, Chart::Transparent, -1, Chart::Transparent, 0xcccccc);

	// Add a legend box with horizontal layout above the plot area at (70, 35). Use 12pt Arial font,
	// transparent background and border, and line style legend icon.
	LegendBox* b = c->addLegend(20, 5, false, "arial.ttf", 12);
	b->setBackground(Chart::Transparent, Chart::Transparent);
	b->setLineStyleKey();

	// Set axis label font to 12pt Arial
	c->xAxis()->setLabelStyle("arial.ttf", 12);
	c->yAxis()->setLabelStyle("arial.ttf", 12);

	// Set the x and y axis stems to transparent, and the x-axis tick color to grey (0xaaaaaa)
	c->xAxis()->setColors(Chart::TextColor, Chart::TextColor, Chart::TextColor, 0xaaaaaa);
	c->yAxis()->setColors(Chart::TextColor);

	// Set the major/minor tick lengths for the x-axis to 10 and 0.
	c->xAxis()->setTickLength(10, 0);

	// For the automatic axis labels, set the minimum spacing to 80/40 pixels for the x/y axis.
	c->xAxis()->setTickDensity(80);
	c->yAxis()->setTickDensity(40);

	// Add a title to the y axis using dark grey (0x555555) 14pt Arial font
	c->yAxis()->setTitle("", "arial.ttf", 14, 0x555555);

	// Add a line layer to the chart with 3-pixel line width
	LineLayer* layer = c->addLineLayer();
	layer->setLineWidth(3);

	// Add 3 data series to the line layer
	layer->addDataSet(Arr_dataReal, 0x5588cc, "Real");
	layer->addDataSet(Arr_dataImag, 0xee9944, "Imag");
	// The x-coordinates for the line layer
	layer->setXData(timeStamps);
	viewer_num.setChart(c);
	delete c;
}
void CSignalsDlg::ViewerDrawSpectrum(vector<double>& data, int Xmax, CChartViewer& viewer_num)
{
	// In this example, we simply use random data for the 3 data series.
	DoubleArray Arr_dataReal = vectorToArray(data);
	vector<double>Datatime;
	double step = sampling / data.size();
	for (int i = 0; i < Xmax; i++)Datatime.push_back(i * step);
	DoubleArray timeStamps = vectorToArray(Datatime);

	// Create a XYChart object of size 600 x 400 pixels
	XYChart* c = new XYChart(850, 240);

	// Add a title box using grey (0x555555) 20pt Arial font
	//c->addTitle("", "arial.ttf", 20, 0x555555);

	// Set the plotarea at (70, 70) and of size 500 x 300 pixels, with transparent background and
	// border and light grey (0xcccccc) horizontal grid lines
	c->setPlotArea(70, 50, 700, 120, Chart::Transparent, -1, Chart::Transparent, 0xcccccc);

	// Add a legend box with horizontal layout above the plot area at (70, 35). Use 12pt Arial font,
	// transparent background and border, and line style legend icon.
	LegendBox* b = c->addLegend(20, 5, false, "arial.ttf", 12);
	b->setBackground(Chart::Transparent, Chart::Transparent);
	b->setLineStyleKey();

	// Set axis label font to 12pt Arial
	c->xAxis()->setLabelStyle("arial.ttf", 12);
	c->yAxis()->setLabelStyle("arial.ttf", 12);

	// Set the x and y axis stems to transparent, and the x-axis tick color to grey (0xaaaaaa)
	c->xAxis()->setColors(Chart::TextColor, Chart::TextColor, Chart::TextColor, 0xaaaaaa);
	c->yAxis()->setColors(Chart::TextColor);

	// Set the major/minor tick lengths for the x-axis to 10 and 0.
	c->xAxis()->setTickLength(10, 0);

	// For the automatic axis labels, set the minimum spacing to 80/40 pixels for the x/y axis.
	c->xAxis()->setTickDensity(80);
	c->yAxis()->setTickDensity(40);

	// Add a title to the y axis using dark grey (0x555555) 14pt Arial font
	c->yAxis()->setTitle("", "arial.ttf", 14, 0x555555);

	// Add a line layer to the chart with 3-pixel line width
	LineLayer* layer = c->addLineLayer();
	layer->setLineWidth(3);

	// Add 3 data series to the line layer
	layer->addDataSet(Arr_dataReal, 0x5588cc, "Real");
	// The x-coordinates for the line layer
	layer->setXData(timeStamps);
	viewer_num.setChart(c);
	delete c;
}
void CSignalsDlg::ViewerDraw(vector<double>& data, int Xmax, CChartViewer& viewer_num)
{
	// In this example, we simply use random data for the 3 data series.
	DoubleArray Arr_dataReal = vectorToArray(data);
	vector<double>Datatime;
	int bits = 1;
	double step = 1;
	if (data.size() != 0)
	{
		bits = data.size() / sp.bit_time;
		step = (double)bits / data.size();
	}
	for (int i = 0; i < Xmax; i++)Datatime.push_back(i * step);
	DoubleArray timeStamps = vectorToArray(Datatime);

	// Create a XYChart object of size 600 x 400 pixels
	XYChart* c = new XYChart(850, 240);

	// Add a title box using grey (0x555555) 20pt Arial font
	//c->addTitle("", "arial.ttf", 20, 0x555555);

	// Set the plotarea at (70, 70) and of size 500 x 300 pixels, with transparent background and
	// border and light grey (0xcccccc) horizontal grid lines
	c->setPlotArea(70, 50, 700, 120, Chart::Transparent, -1, Chart::Transparent, 0xcccccc);

	// Add a legend box with horizontal layout above the plot area at (70, 35). Use 12pt Arial font,
	// transparent background and border, and line style legend icon.
	LegendBox* b = c->addLegend(20, 5, false, "arial.ttf", 12);
	b->setBackground(Chart::Transparent, Chart::Transparent);
	b->setLineStyleKey();

	// Set axis label font to 12pt Arial
	c->xAxis()->setLabelStyle("arial.ttf", 12);
	c->yAxis()->setLabelStyle("arial.ttf", 12);

	// Set the x and y axis stems to transparent, and the x-axis tick color to grey (0xaaaaaa)
	c->xAxis()->setColors(Chart::TextColor, Chart::TextColor, Chart::TextColor, 0xaaaaaa);
	c->yAxis()->setColors(Chart::TextColor);

	// Set the major/minor tick lengths for the x-axis to 10 and 0.
	c->xAxis()->setTickLength(10, 0);

	// For the automatic axis labels, set the minimum spacing to 80/40 pixels for the x/y axis.
	c->xAxis()->setTickDensity(80);
	c->yAxis()->setTickDensity(40);

	// Add a title to the y axis using dark grey (0x555555) 14pt Arial font
	c->yAxis()->setTitle("", "arial.ttf", 14, 0x555555);

	// Add a line layer to the chart with 3-pixel line width
	LineLayer* layer = c->addLineLayer();
	layer->setLineWidth(3);

	// Add 3 data series to the line layer
	layer->addDataSet(Arr_dataReal, 0x5588cc, "Real");
	// The x-coordinates for the line layer
	layer->setXData(timeStamps);
	viewer_num.setChart(c);
	delete c;
}
void CSignalsDlg::ViewerDraw(vector<vector<double>>& data, double Xmin, double Xmax, CChartViewer& viewer_num, string PathPic)
{
	if (data.empty())return;
	// In this example, we simply use random data for the 3 data series.
	vector<DoubleArray> Arr_dataReal; Arr_dataReal.resize(data.size());
	for (int i = 0; i < data.size(); i++)
	{
		Arr_dataReal[i] = vectorToArray(data[i]);
	}

	vector<double>Datatime;

	double OXstep = (Xmax - Xmin) / (data[0].size() - 1);
	for (double i = Xmin; i <= Xmax; i += OXstep)Datatime.push_back(i);
	DoubleArray timeStamps = vectorToArray(Datatime);

	// Create a XYChart object of size 600 x 400 pixels
	XYChart* c = new XYChart(850, 240);

	// Add a title box using grey (0x555555) 20pt Arial font
	//c->addTitle("", "arial.ttf", 20, 0x555555);

	// Set the plotarea at (70, 70) and of size 500 x 300 pixels, with transparent background and
	// border and light grey (0xcccccc) horizontal grid lines
	c->setPlotArea(70, 50, 700, 120, Chart::Transparent, -1, Chart::Transparent, 0xcccccc);

	// Add a legend box with horizontal layout above the plot area at (70, 35). Use 12pt Arial font,
	// transparent background and border, and line style legend icon.
	LegendBox* b = c->addLegend(20, 5, false, "arial.ttf", 12);
	b->setBackground(Chart::Transparent, Chart::Transparent);
	b->setLineStyleKey();

	// Set axis label font to 12pt Arial
	c->xAxis()->setLabelStyle("arial.ttf", 12);
	c->yAxis()->setLabelStyle("arial.ttf", 12);

	// Set the x and y axis stems to transparent, and the x-axis tick color to grey (0xaaaaaa)
	c->xAxis()->setColors(Chart::TextColor, Chart::TextColor, Chart::TextColor, 0xaaaaaa);
	c->yAxis()->setColors(Chart::TextColor);

	// Set the major/minor tick lengths for the x-axis to 10 and 0.
	c->xAxis()->setTickLength(10, 0);

	// For the automatic axis labels, set the minimum spacing to 80/40 pixels for the x/y axis.
	c->xAxis()->setTickDensity(80);
	c->yAxis()->setTickDensity(40);

	// Add a title to the y axis using dark grey (0x555555) 14pt Arial font
	c->yAxis()->setTitle("", "arial.ttf", 14, 0x555555);

	// Add a line layer to the chart with 3-pixel line width
	LineLayer* layer = c->addLineLayer();
	layer->setLineWidth(3);
	//
	layer->setDataLabelFormat("{value|1} ");
	// Add 3 data series to the line layer
	for (int i = 0; i < Arr_dataReal.size(); i++)
	{
		stringstream ss;
		ss << "Data " << i;
		string name = ss.str();
		int color = 0. + 16777216 * rand() / RAND_MAX;
		layer->addDataSet(Arr_dataReal[i], color, name.c_str());
	}
	// The x-coordinates for the line layer
	layer->setXData(timeStamps);
	viewer_num.setChart(c);
	const char* chPathPic = PathPic.c_str();
	c->makeChart(chPathPic);
	delete c;
}

//Обновление переменных класса sp - класса обработки генерируемых сигналов
void CSignalsDlg::updateSP() //Обновление переменных класса sp - класса обработки генерируемых сигналов
{
	sp.BrV = bitrate;
	sp.sampling = sampling;
}

//Exit
void CSignalsDlg::OnBnClickedCancel()
{
	CDialogEx::OnCancel();
}

//Генерация сигналов
void CSignalsDlg::OnBnClickedButton1()
{
	UpdateData(1);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	Signals_Gen(bits_size, delay_size, noize_lvl);
	OnBnClickedCheck1();
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(0);
}

//Сигналы/Спекрты
void CSignalsDlg::OnBnClickedCheck1() //Сигналы/Спекрты
{
	UpdateData(1);
	if (Signals_or_Spectrs)
	{
		vector <double> ImSpectr1Mod;
		vector <double> ImSpectr2Mod;
		for (int i = 0; i < ImSpectr1.size(); i++)
		{
			ImSpectr1Mod.push_back(sqrt(pow(ImSpectr1[i].real(), 2) + pow(ImSpectr1[i].real(), 2)));
		}
		for (int i = 0; i < ImSpectr2.size(); i++)
		{
			ImSpectr2Mod.push_back(sqrt(pow(ImSpectr2[i].real(), 2) + pow(ImSpectr2[i].real(), 2)));
		}
		ViewerDrawSpectrum(ImSpectr1Mod, ImSpectr1Mod.size(), viewer1);
		ViewerDrawSpectrum(ImSpectr2Mod, ImSpectr2Mod.size(), viewer2);
	}
	else
	{
		ViewerDraw(ImSignal1, ImSignal1.size(), viewer1);
		ViewerDraw(ImSignal2, ImSignal2.size(), viewer2);
	}
	UpdateData(0);
}

//Функция Uncertainty
void CSignalsDlg::OnBnClickedButton3()//Функция Uncertainty
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	ResearchRrr.clear();
	auto start = steady_clock::now();
	sp.Uncertainty(ResearchRrr, ImSignal1, ImSignal2, _k);
	auto end = steady_clock::now();
	auto elapsed = duration_cast<milliseconds>(end - start);
	test_time_cr = elapsed.count();
	if (ResearchRrr.size() != NULL)
	{
		double buff_ResearchRrr = 0;
		for (int i = 0; i < ResearchRrr.size(); i++)
		{
			if (ResearchRrr[i] > buff_ResearchRrr)
			{
				buff_ResearchRrr = ResearchRrr[i]; delay_lama = i;
			}
		}
	}
	delay_lama = int((double)delay_lama / sp.bit_time);
	ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
}

//Функция Uncertainty OMP
void CSignalsDlg::OnBnClickedButton5()
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	ResearchRrr.clear();

	auto start = steady_clock::now();
	sp.Uncertainty_omp(ResearchRrr, ImSignal1, ImSignal2, _k);
	auto end = steady_clock::now();
	auto elapsed = duration_cast<milliseconds>(end - start);
	test_time_cr = elapsed.count();
	if (ResearchRrr.size() != NULL)
	{
		double buff_ResearchRrr = 0;
		for (int i = 0; i < ResearchRrr.size(); i++)
		{
			if (ResearchRrr[i] > buff_ResearchRrr)
			{
				buff_ResearchRrr = ResearchRrr[i]; delay_lama = i;
			}
		}
	}
	delay_lama = int((double)delay_lama / sp.bit_time);
	ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
}

//Функция Uncertainty IPP
void CSignalsDlg::OnBnClickedButton9()
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	ResearchRrr.clear();

	auto start = steady_clock::now();
	sp.Uncertainty_ipp(ResearchRrr, ImSignal1, ImSignal2, _k);
	auto end = steady_clock::now();
	auto elapsed = duration_cast<milliseconds>(end - start);
	test_time_cr = elapsed.count();
	if (ResearchRrr.size() != NULL)
	{
		double buff_ResearchRrr = 0;
		for (int i = 0; i < ResearchRrr.size(); i++)
		{
			if (ResearchRrr[i] > buff_ResearchRrr)
			{
				buff_ResearchRrr = ResearchRrr[i]; delay_lama = i;
			}
		}
	}
	delay_lama = int((double)delay_lama / sp.bit_time);
	ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
}

//Генератор фильтров
void CSignalsDlg::OnBnClickedButton4()
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	//////////////////// Создание фильтров
	sp.fir_s.clear();
	sp.fir_s.resize(sp.operating_frequencies.size());
	for (int i = 0; i < sp.operating_frequencies.size(); i++)
	{
		double freq1 = ((sp.operating_frequencies[i] - 912.5) * pow(10, 6) - (sp.BrV / 4));
		double freq2 = ((sp.operating_frequencies[i] - 912.5) * pow(10, 6) + (sp.BrV / 4));
		//double trans_zone = (sp.BrV / 8.) * (3./2.); /// ПОЧЕМУ!!!! НЕ РАБОТАЕТ ПРИ sp.BrV / 8.
		double trans_zone = (sp.BrV / 4.);
		CreateFirFilter(freq1, freq2, trans_zone, sampling, sp.fir_s[i]);
	}
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
}

//ФН для ППРЧ сигналов
void CSignalsDlg::OnBnClickedButton2()
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	int found_delay;
	double pi = sp.Uncertainty_ipp_jtids(delay_size, ImSignal1, ImSignal2, _k, ResearchRrr, found_delay, delay_lama);
	ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
}

//Исследование
void CSignalsDlg::OnBnClickedButton6()
{
	UpdateData(1);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	veroiatnosti_fhss.clear();
	veroiatnosti_un.clear();
	double noize_min_r = -20;
	double noize_max_r = -10;
	int noize_dots_r = 10;
	double noize_step_r = (noize_max_r - noize_min_r) / (noize_dots_r - 1);
	vector<double> assessments_fhss;
	vector<double> assessments_un;
	assessments_fhss.resize(noize_dots_r);
	veroiatnosti_fhss.resize(noize_dots_r);
	assessments_un.resize(noize_dots_r);
	veroiatnosti_un.resize(noize_dots_r);
	int runs = 20; //кол-во прогонов для одной велечины шума
	for (int i = 0; i < noize_dots_r; i++)
	{
		double noize_r = noize_min_r + i * noize_step_r;
		for (int j = 0; j < runs; j++)
		{
			/////////////////////////////////////////////////Генерация новых сигналов для итерации исследования			
			Signals_Gen(bits_size, delay_size, noize_r);
			/////////////////////////////////////////////////
			//////////////////// Подготовка входных сигналов
			signal_buf signal1;
			signal_buf signal2;
			for (int i = 0; i < ImSignal1.size(); i++)
			{
				complex<float> buf; buf = ImSignal1[i].real() + comjd * ImSignal1[i].imag();
				signal1.push_back(buf);
			}
			for (int i = 0; i < ImSignal2.size(); i++)
			{
				complex<float> buf; buf = ImSignal2[i].real() + comjd * ImSignal2[i].imag();
				signal2.push_back(buf);
			}
			//////////////////// fast_convolution
			fast_convolution(signal1, sp.fir_s, sp.FHSS_Signals_initial_fl, GPU_FD);
			fast_convolution(signal2, sp.fir_s, sp.FHSS_Signals_fl, GPU_FD);
			ResearchRrr.clear();
			ResearchRrr2D.clear();
			ResearchRrr2D.resize(sp.operating_frequencies.size());
#pragma omp parallel for
			for (int i = 0; i < sp.operating_frequencies.size(); i++)
			{
				vector<float>buffer;
				sp.Uncertainty_omp(buffer, sp.FHSS_Signals_initial_fl[i], sp.FHSS_Signals_fl[i], _k);
				ResearchRrr2D[i] = buffer;
			}
			ResearchRrr.resize(ResearchRrr2D[0].size());
			for (int j = 0; j < sp.operating_frequencies.size(); j++)
			{
				for (int i = 0; i < ResearchRrr.size(); i++)
				{
					ResearchRrr[i] += ResearchRrr2D[j][i];
				}
			}
			if (ResearchRrr.size() != NULL)
			{
				double buff_ResearchRrr = 0;
				for (int i = 0; i < ResearchRrr.size(); i++)
				{
					if (ResearchRrr[i] > buff_ResearchRrr)
					{
						buff_ResearchRrr = ResearchRrr[i]; delay_lama = i;
					}
				}
			}
			delay_lama_r = delay_lama;
			delay_lama = int((double)delay_lama / sp.bit_time);
			sp.FHSS_Signals_initial_fl.clear();
			sp.FHSS_Signals_fl.clear();
			ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
			/////////////////////////////////////////////////
			double expected_delay = delay_size * sp.bit_time;
			assessments_fhss[i] += pow((expected_delay - delay_lama_r), 2);
			/////////////////////////////////////////////////
			double deltadelay = abs(expected_delay - delay_lama_r);
			if (deltadelay < (double(sp.bit_time) / 2)) veroiatnosti_fhss[i].second += 1;
			//////////////////////////////////////////////////
			sp.Uncertainty_omp(ResearchRrr, ImSignal1, ImSignal2, _k);
			if (ResearchRrr.size() != NULL)
			{
				double buff_ResearchRrr = 0;
				for (int i = 0; i < ResearchRrr.size(); i++)
				{
					if (ResearchRrr[i] > buff_ResearchRrr)
					{
						buff_ResearchRrr = ResearchRrr[i]; delay_lama = i;
					}
				}
			}
			delay_lama_r = delay_lama;
			delay_lama = int((double)delay_lama / sp.bit_time);
			deltadelay = abs(expected_delay - delay_lama_r);
			if (deltadelay < (double(sp.bit_time) / 2)) veroiatnosti_un[i].second += 1;
			//////////////////////////////////////////////////
		}
		assessments_fhss[i] /= runs;
		veroiatnosti_fhss[i].second /= runs;
		veroiatnosti_fhss[i].first = noize_r;
		veroiatnosti_un[i].second /= runs;
		veroiatnosti_un[i].first = noize_r;
	}
	//double test_pic = 1.2345;
	//string str = convertToStrPng<double>(&test_pic);
	//ViewerDraw(assessments_fhss, noize_min_r, noize_max_r, viewer3, "assessments.png");
	vectorDoubleToFile(veroiatnosti_fhss, "veroiatnosti.txt");
	vector<vector<double>> veroiatnosti_help; veroiatnosti_help.resize(2);
	for (int i = 0; i < veroiatnosti_fhss.size(); i++)
	{
		veroiatnosti_help[0].push_back(veroiatnosti_fhss[i].second);
		veroiatnosti_help[1].push_back(veroiatnosti_un[i].second);
	}
	ViewerDraw(veroiatnosti_help, noize_min_r, noize_max_r, viewer3, "veroiatnosti.png");
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(0);
}

//Исследвование производительности
void CSignalsDlg::OnBnClickedButton7()
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	vector<vector<double>> times;
	times.resize(2); //0-noMod 1-ompMod 2-ippMod

	int size_start = 100;
	int size_stop = size_start + 1000;
	int size_step = (size_stop - size_start) / 10;
	int povtor = 5;
	for (int i = size_start; i < size_stop; i += size_step)
	{
		u_int64 buffer_time = 0;
		for (int j = 0; j < povtor; j++)
		{
			Signals_Gen(i, i / 5, noize_lvl);
			ResearchRrr.clear();
			auto start = steady_clock::now();
			sp.Uncertainty(ResearchRrr, ImSignal1, ImSignal2, _k);
			auto end = steady_clock::now();
			auto elapsed = duration_cast<milliseconds>(end - start);
			buffer_time += elapsed.count();
		}
		buffer_time /= povtor;
		times[0].push_back(buffer_time);
	}
	for (int i = size_start; i < size_stop; i += size_step)
	{
		u_int64 buffer_time = 0;
		for (int j = 0; j < povtor; j++)
		{
			Signals_Gen(i, i / 5, noize_lvl);
			ResearchRrr.clear();
			auto start = steady_clock::now();
			sp.Uncertainty_ipp(ResearchRrr, ImSignal1, ImSignal2, _k);
			auto end = steady_clock::now();
			auto elapsed = duration_cast<milliseconds>(end - start);
			buffer_time += elapsed.count();
		}
		buffer_time /= povtor;
		times[1].push_back(buffer_time);
	}
	auto str_time = fh.get_time_str();
	str_time = "performance_" + str_time + ".png";
	str_time = LOGS_PATH + str_time;
	ViewerDraw(times, size_start, size_stop, viewer3, str_time);
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
	/////////////////////////////
	/////////////////////////////
	/////////////////////////////
	//UpdateData(TRUE);
	//SetCursor(LoadCursor(nullptr, IDC_WAIT));
	//updateSP();
	//vector<vector<double>> times;
	//times.resize(2); //0-noMod 1-ompMod 2-ippMod

	//int k_start = pow(2,0);
	//int k_stop = pow(2, 7);
	//int povtor = 5;
	//for (int i = k_start; i <= k_stop; i *= 2)
	//{
	//	u_int64 buffer_time = 0;
	//	for (int j = 0; j < povtor; j++)
	//	{
	//		Signals_Gen(bits_size, delay_size, noize_lvl);
	//		ResearchRrr.clear();
	//		auto start = steady_clock::now();
	//		sp.Uncertainty(ResearchRrr, ImSignal1, ImSignal2, i);
	//		auto end = steady_clock::now();
	//		auto elapsed = duration_cast<milliseconds>(end - start);
	//		buffer_time += elapsed.count();
	//	}
	//	buffer_time /= povtor;
	//	times[0].push_back(buffer_time);
	//}
	//for (int i = k_start; i <= k_stop; i *= 2)
	//{
	//	u_int64 buffer_time = 0;
	//	for (int j = 0; j < povtor; j++)
	//	{
	//		Signals_Gen(bits_size, delay_size, noize_lvl);
	//		ResearchRrr.clear();
	//		auto start = steady_clock::now();
	//		sp.Uncertainty_ipp(ResearchRrr, ImSignal1, ImSignal2, i);
	//		auto end = steady_clock::now();
	//		auto elapsed = duration_cast<milliseconds>(end - start);
	//		buffer_time += elapsed.count();
	//	}
	//	buffer_time /= povtor;
	//	times[1].push_back(buffer_time);
	//}
	//auto str_time = fh.get_time_str();
	//str_time = "performance_k_" + str_time + ".png";
	//str_time = LOGS_PATH + str_time;
	//ViewerDraw(times, 0, 7, viewer3, str_time);
	//SetCursor(LoadCursor(nullptr, IDC_ARROW));
	//UpdateData(FALSE);
}

//Исследование определения ВВЗ функцией неопределённости по каналам.
void CSignalsDlg::OnBnClickedButton8()
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();

	/*vector<int> delay_error_v;
	vector<double> peak_intensity_v;*/

	veroiatnosti_fhss.clear();
	//veroiatnosti_un.clear();
	double noize_min_r = -30;
	double noize_max_r = -25;
	int noize_dots_r = 10;
	double noize_step_r = (noize_max_r - noize_min_r) / (noize_dots_r - 1);
	veroiatnosti_fhss.resize(noize_dots_r);
	//veroiatnosti_un.resize(noize_dots_r);

	for (int i = 0; i < noize_dots_r; i++)
	{
		double noize_r = noize_min_r + i * noize_step_r;
		int try_size = 5;
		for (int j = 0; j < try_size; j++)
		{
			Signals_Gen(bits_size, bits_size / 10, noize_r);
			int found_delay;
			double pi = sp.Uncertainty_ipp_jtids(delay_size, ImSignal1, ImSignal2, _k, ResearchRrr, found_delay, delay_lama);
			double expected_delay = delay_size * sp.bit_time;
			double deltadelay = abs(expected_delay - found_delay);
			if ((deltadelay < (double(sp.bit_time) / 2)) && pi > 0.008)
			{
				veroiatnosti_fhss[i].second += 1;
			}
			/*if ((deltadelay < (double(sp.bit_time) / 2)))
			{
				veroiatnosti_un[i].second += 1;
			}*/
			/*peak_intensity_v.push_back(pi);
			auto file_name = fh.get_time_str();
			file_name = "Study2020_" + file_name;
			file_name += "_noizelvl_" + std::to_string((int)noize_lvl);
			file_name += "signalSize_" + std::to_string(ImSignal2.size());
			file_name += "k" + std::to_string(_k);
			file_name = LOGS_PATH + file_name;
			fh.write_vector_to_file(peak_intensity_v, file_name);*/
		}
		veroiatnosti_fhss[i].second /= try_size;
		veroiatnosti_fhss[i].first = noize_r;
		/*veroiatnosti_un[i].second /= try_size;
		veroiatnosti_un[i].first = noize_r;*/
		auto file_name = fh.get_time_str();
		file_name = "veroiatnosti_fhss_" + file_name;
		file_name += "_noizelvl_" + std::to_string((int)noize_r);
		file_name += "_signalSize_" + std::to_string(ImSignal2.size());
		file_name += "_k_" + std::to_string(_k);
		file_name += "_try_size_" + std::to_string(try_size);
		file_name += "_vector_size_" + std::to_string(try_size);
		file_name += ".txt";
		file_name = LOGS_PATH + file_name;
		vectorDoubleToFile(veroiatnosti_fhss, file_name);
	}
	vectorDoubleToFile(veroiatnosti_fhss, "veroiatnostiPro.txt");
	vector<vector<double>> veroiatnosti_help; veroiatnosti_help.resize(2);
	for (int i = 0; i < veroiatnosti_fhss.size(); i++)
	{
		veroiatnosti_help[0].push_back(veroiatnosti_fhss[i].second);
		veroiatnosti_help[1].push_back(veroiatnosti_un[i].second);
	}
	ViewerDraw(veroiatnosti_help, noize_min_r, noize_max_r, viewer3, "veroiatnosti.png");
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(0);
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
}



//VIcourse
void CSignalsDlg::TrueViewerDraw(vector<vector<double>>& data, double Xmin, double Xmax, CChartViewer& viewer_num, string PathPic, bool podpisi)
{
	if (data.empty())return;
	// In this example, we simply use random data for the 3 data series.
	vector<DoubleArray> Arr_dataReal; Arr_dataReal.resize(data.size());
	for (int i = 0; i < data.size(); i++)
	{
		Arr_dataReal[i] = vectorToArray(data[i]);
	}

	vector<double>Datatime;

	double OXstep = (Xmax - Xmin) / (data[0].size() - 1);
	for (double i = Xmin; i <= Xmax; i += OXstep)Datatime.push_back(i);
	DoubleArray timeStamps = vectorToArray(Datatime);

	// Create a XYChart object of size 600 x 400 pixels
	XYChart* c = new XYChart(850, 240);

	// Add a title box using grey (0x555555) 20pt Arial font
	//c->addTitle("", "arial.ttf", 20, 0x555555);

	// Set the plotarea at (70, 70) and of size 500 x 300 pixels, with transparent background and
	// border and light grey (0xcccccc) horizontal grid lines
	c->setPlotArea(70, 50, 700, 120, Chart::Transparent, -1, Chart::Transparent, 0xcccccc);

	// Add a legend box with horizontal layout above the plot area at (70, 35). Use 12pt Arial font,
	// transparent background and border, and line style legend icon.
	LegendBox* b = c->addLegend(20, 5, false, "arial.ttf", 12);
	b->setBackground(Chart::Transparent, Chart::Transparent);
	b->setLineStyleKey();

	// Set axis label font to 12pt Arial
	c->xAxis()->setLabelStyle("arial.ttf", 12);
	c->yAxis()->setLabelStyle("arial.ttf", 12);

	// Set the x and y axis stems to transparent, and the x-axis tick color to grey (0xaaaaaa)
	c->xAxis()->setColors(Chart::TextColor, Chart::TextColor, Chart::TextColor, 0xaaaaaa);
	c->yAxis()->setColors(Chart::TextColor);

	// Set the major/minor tick lengths for the x-axis to 10 and 0.
	c->xAxis()->setTickLength(10, 0);

	// For the automatic axis labels, set the minimum spacing to 80/40 pixels for the x/y axis.
	c->xAxis()->setTickDensity(80);
	c->yAxis()->setTickDensity(40);

	// Add a title to the y axis using dark grey (0x555555) 14pt Arial font
	c->yAxis()->setTitle("", "arial.ttf", 14, 0x555555);

	// Add a line layer to the chart with 3-pixel line width
	LineLayer* layer = c->addLineLayer();
	layer->setLineWidth(3);
	//
	if (podpisi) layer->setDataLabelFormat("{value|2} ");
	// Add 3 data series to the line layer
	for (int i = 0; i < Arr_dataReal.size(); i++)
	{
		stringstream ss;
		ss << "Data " << i;
		string name = ss.str();
		layer->addDataSet(Arr_dataReal[i], -1, name.c_str());
	}
	// The x-coordinates for the line layer
	layer->setXData(timeStamps);
	viewer_num.setChart(c);
	const char* chPathPic = PathPic.c_str();
	c->makeChart(chPathPic);
	delete c;
}
void CSignalsDlg::OnBnClickedButton10()
{	
	UpdateData(1);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	MySignals.Init(sampling, 1119e6, bitrate, bits_size, noize_lvl, 2, 0.3);
	MySignals.GetSignals();
	sp.addNoize(MySignals.CSignal1, noize_lvl);
	sp.addNoize(MySignals.CSignal2, noize_lvl);
	ImSignal1 = MySignals.CSignal1;
	ImSignal2 = MySignals.CSignal2;
	ImSpectr1.clear();
	ImSpectr2.clear();
	sp.FAST_FUR(ImSignal1, ImSpectr1, -1);
	sp.FAST_FUR(ImSignal2, ImSpectr2, -1);
	sp.spVertex(ImSpectr1);
	sp.spVertex(ImSpectr2);
	OnBnClickedCheck1();
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(0);
}

void CSignalsDlg::OnBnClickedButton11()
{
	sp.nonlinear_filtering(ImSignal1, 1119e6, sampling, bitrate);
	sp.nonlinear_filtering(ImSignal2, 1119e6, sampling, bitrate);
	ImSpectr1.clear();
	ImSpectr2.clear();
	sp.FAST_FUR(ImSignal1, ImSpectr1, -1);
	sp.FAST_FUR(ImSignal2, ImSpectr2, -1);
	sp.spVertex(ImSpectr1);
	sp.spVertex(ImSpectr2);
	OnBnClickedCheck1();
}
