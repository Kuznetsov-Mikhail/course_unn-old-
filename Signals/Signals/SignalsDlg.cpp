
// SignalsDlg.cpp: файл реализации
//

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
	, bitrate(5000000)
	, sampling(350000000)
	, bits_size(1700)
	, delay_size(200)
	, delay_lama(0)
	, Signals_or_Spectrs(FALSE)
	, noize_lvl(100)
	, f_dop(30000)
	, alfa(2e-05)
	, Dopler_On(FALSE)
	, scramble(TRUE)
	, parallel(TRUE)
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
	DDX_Text(pDX, IDC_EDIT8, alfa);
	DDX_Check(pDX, IDC_CHECK2, Dopler_On);
	DDX_Check(pDX, IDC_CHECK3, scramble);
	DDX_Check(pDX, IDC_CHECK4, parallel);
	DDX_Control(pDX, IDC_PROGRESS1, Prog_bar);
	DDX_Control(pDX, IDC_PROGRESS2, ProgBarRes);
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
	Prog_bar.SetRange(0, 50);
	Prog_bar.SetPos(0);
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
	int bits=1;
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
void CSignalsDlg::ViewerDraw(vector<double>& data, double Xmin, double Xmax, CChartViewer& viewer_num, string PathPic)
{
	// In this example, we simply use random data for the 3 data series.
	DoubleArray Arr_dataReal = vectorToArray(data);
	vector<double>Datatime;

	double OXstep = (Xmax - Xmin) / (data.size()-1);
	for (double i = Xmin; i <= Xmax; i+= OXstep)Datatime.push_back(i);
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
	const char* chPathPic = PathPic.c_str();
	c->makeChart(chPathPic);
	delete c;
}
void CSignalsDlg::updateSP() //Обновление переменных класса sp - класса обработки генерируемых сигналов
{
	sp.BrV = bitrate;
	sp.sampling = sampling;
}

void CSignalsDlg::OnBnClickedCancel()
{
	// TODO: добавьте свой код обработчика уведомлений
	CDialogEx::OnCancel();
}


void CSignalsDlg::OnBnClickedButton1() //Генерация сигналов
{
	UpdateData(1);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	ResearchRrr.resize(0);
	ImSpectr1.clear();
	ImSpectr2.clear();

	sp.Link16_Signals_Generator(ImSignal1, ImSignal2, bits_size, delay_size, scramble);
	ImSignal2Rch.clear();
	ImSignal2Rch = ImSignal2;

	if (Dopler_On)
	{
		sp.Dopler_scaling(ImSignal2, alfa);
		sp.Dopler_shift(ImSignal2, f_dop);
	}
	sp.addNoize(ImSignal2, noize_lvl);
	sp.FAST_FUR(ImSignal1, ImSpectr1, -1);
	sp.FAST_FUR(ImSignal2, ImSpectr2, -1);
	sp.spVertex(ImSpectr1);
	sp.spVertex(ImSpectr2);
	OnBnClickedCheck1();
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(0);
}


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
void CSignalsDlg::OnBnClickedButton3()//Функция Uncertainty
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	ResearchRrr.clear();
	sp.Uncertainty(ResearchRrr, ImSignal1, ImSignal2, 4);

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
void CSignalsDlg::OnBnClickedButton5() //Функция Uncertainty OMP
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	updateSP();
	ResearchRrr.clear();
	sp.Uncertainty_omp(ResearchRrr, ImSignal1, ImSignal2, 8);

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

void CSignalsDlg::OnBnClickedButton4() //Генератор фильтров
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
		Prog_bar.SetPos(i + 1);
	}
	Prog_bar.SetPos(0);
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
}
void CSignalsDlg::OnBnClickedButton2() //ФН для ППРЧ сигналов
{
	UpdateData(TRUE);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	Prog_bar.SetRange(0, 50);
	Prog_bar.SetPos(0);
	updateSP();
	if (!parallel)
	{		
		vector <complex<double>> ImSpectr1_local;
		vector <complex<double>> ImSpectr2_local;

		sp.FAST_FUR(ImSignal1, ImSpectr1_local, -1);
		sp.FAST_FUR(ImSignal2, ImSpectr2_local, -1);

		sp.spVertex(ImSpectr1_local);
		sp.spVertex(ImSpectr2_local);

		sp.FHSS_Signals_initial.clear();
		sp.FHSS_Signals_initial.resize(sp.operating_frequencies.size());
		sp.FHSS_Signals.clear();
		sp.FHSS_Signals.resize(sp.operating_frequencies.size());
		for (int i = 0; i < sp.operating_frequencies.size(); i++)
		{
			int start = ((double)ImSpectr2_local.size() / sampling) * ((sp.operating_frequencies[i] - 912.5) * pow(10, 6) - (sp.BrV / 4));
			int stop = ((double)ImSpectr2_local.size() / sampling) * ((sp.operating_frequencies[i] - 912.5) * pow(10, 6) + (sp.BrV / 4));
			sp.FHSS_Signals[i].resize(ImSpectr2_local.size());
			for (int j = start; j < stop; j++)
			{
				sp.FHSS_Signals[i][j] = ImSpectr2_local[j];
			}
			sp.spVertex(sp.FHSS_Signals[i]);
			vector <complex<double>> signal_buffer = sp.FHSS_Signals[i];
			sp.FAST_FUR(signal_buffer, sp.FHSS_Signals[i], 1);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			start = ((double)ImSpectr1_local.size() / sampling) * ((sp.operating_frequencies[i] - 912.5) * pow(10, 6) - (sp.BrV / 4));
			stop = ((double)ImSpectr1_local.size() / sampling) * ((sp.operating_frequencies[i] - 912.5) * pow(10, 6) + (sp.BrV / 4));
			sp.FHSS_Signals_initial[i].resize(ImSpectr2_local.size());
			for (int j = start; j < stop; j++)
			{
				sp.FHSS_Signals_initial[i][j] = ImSpectr2_local[j];
			}
			sp.spVertex(sp.FHSS_Signals_initial[i]);
			signal_buffer = sp.FHSS_Signals_initial[i];
			sp.FAST_FUR(signal_buffer, sp.FHSS_Signals_initial[i], 1);
		}
		ResearchRrr.clear();
		//ResearchRrr2D.clear();
		for (int i = 0; i < sp.operating_frequencies.size(); i++)
		{
			vector<double>buffer;
			sp.Uncertainty(buffer, sp.FHSS_Signals_initial[i], sp.FHSS_Signals[i], 4);
			//ResearchRrr2D.push_back(buffer);
			ResearchRrr.resize(buffer.size());
			for (int j = 0; j < ResearchRrr.size(); j++) ResearchRrr[j] += buffer[j];
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
		sp.FHSS_Signals_initial.clear();
		sp.FHSS_Signals.clear();
		ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
	}
	else
	{
		//////////////////// Подготовка входных сигналов
		signal_buf signal1;
		signal_buf signal2;
		for (int i = 0; i < ImSignal1.size(); i++)
		{
			complex<float> buf; buf = ImSignal1[i].real() + Comj * ImSignal1[i].imag();
			signal1.push_back(buf);
		}
		for (int i = 0; i < ImSignal2.size(); i++)
		{
			complex<float> buf; buf = ImSignal2[i].real() + Comj * ImSignal2[i].imag();
			signal2.push_back(buf);
		}
		//////////////////// fast_convolution
		fast_convolution(signal1, sp.fir_s, sp.FHSS_Signals_initial_fl, GPU_FD);
		Prog_bar.SetPos(25);
		fast_convolution(signal2, sp.fir_s, sp.FHSS_Signals_fl, GPU_FD);
		Prog_bar.SetPos(50);
		ResearchRrr.clear();
		ResearchRrr2D.clear();
		ResearchRrr2D.resize(sp.operating_frequencies.size());
		#pragma omp parallel for
		for (int i = 0; i < sp.operating_frequencies.size(); i++)
		{
			vector<float>buffer;
			sp.Uncertainty_omp(buffer, sp.FHSS_Signals_initial_fl[i], sp.FHSS_Signals_fl[i], 4);
			ResearchRrr2D[i]= buffer;
			//Prog_bar.SetPos(i+1);
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
		delay_lama = int((double)delay_lama / sp.bit_time);
		sp.FHSS_Signals_initial_fl.clear();
		sp.FHSS_Signals_fl.clear();
		ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
	}
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(FALSE);
}

void CSignalsDlg::OnBnClickedButton6()//Исследование
{
	UpdateData(1);
	SetCursor(LoadCursor(nullptr, IDC_WAIT));
	veroiatnosti.clear();
	double noize_min_r = -20;
	double noize_max_r = 10;
	int noize_dots_r = 5;
	double noize_step_r = (noize_max_r - noize_min_r) / noize_dots_r;
	vector<double> assessments;
	assessments.resize(noize_dots_r+1);
	veroiatnosti.resize(noize_dots_r + 1);
	int runs = 10; //кол-во прогонов для одной велечины шума

	Prog_bar.SetRange(0, runs);
	Prog_bar.SetPos(0);
	ProgBarRes.SetRange(0, noize_dots_r);
	ProgBarRes.SetPos(0);
	for (int i = 0; i < noize_dots_r+1; i++)
	{
		double noize_r = noize_min_r + i * noize_step_r;
		for (int j = 0; j < runs; j++)
		{
			/////////////////////////////////////////////////Генерация новых сигналов для итерации исследования			
			updateSP();
			ResearchRrr.resize(0);
			ImSpectr1.clear();
			ImSpectr2.clear();

			sp.Link16_Signals_Generator(ImSignal1, ImSignal2, bits_size, delay_size, scramble);
			ImSignal2Rch.clear();
			ImSignal2Rch = ImSignal2;

			if (Dopler_On)
			{
				sp.Dopler_scaling(ImSignal2, alfa);
				sp.Dopler_shift(ImSignal2, f_dop);
			}
			sp.addNoize(ImSignal2, noize_r);
			sp.FAST_FUR(ImSignal1, ImSpectr1, -1);
			sp.FAST_FUR(ImSignal2, ImSpectr2, -1);
			sp.spVertex(ImSpectr1);
			sp.spVertex(ImSpectr2);
			/////////////////////////////////////////////////
			//////////////////// Подготовка входных сигналов
			signal_buf signal1;
			signal_buf signal2;
			for (int i = 0; i < ImSignal1.size(); i++)
			{
				complex<float> buf; buf = ImSignal1[i].real() + Comj * ImSignal1[i].imag();
				signal1.push_back(buf);
			}
			for (int i = 0; i < ImSignal2.size(); i++)
			{
				complex<float> buf; buf = ImSignal2[i].real() + Comj * ImSignal2[i].imag();
				signal2.push_back(buf);
			}
			//////////////////// fast_convolution
			fast_convolution(signal1, sp.fir_s, sp.FHSS_Signals_initial_fl, GPU_FD);
			Prog_bar.SetPos(25);
			fast_convolution(signal2, sp.fir_s, sp.FHSS_Signals_fl, GPU_FD);
			Prog_bar.SetPos(50);
			ResearchRrr.clear();
			ResearchRrr2D.clear();
			ResearchRrr2D.resize(sp.operating_frequencies.size());
#pragma omp parallel for
			for (int i = 0; i < sp.operating_frequencies.size(); i++)
			{
				vector<float>buffer;
				sp.Uncertainty_omp(buffer, sp.FHSS_Signals_initial_fl[i], sp.FHSS_Signals_fl[i], 4);
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
			delay_lama_r= delay_lama;
			delay_lama = int((double)delay_lama / sp.bit_time);
			sp.FHSS_Signals_initial_fl.clear();
			sp.FHSS_Signals_fl.clear();
			ViewerDraw(ResearchRrr, ResearchRrr.size(), viewer3);
			/////////////////////////////////////////////////
			double expected_delay = delay_size * sp.bit_time;
			assessments[i] += pow((expected_delay - delay_lama_r), 2);
			/////////////////////////////////////////////////
			double deltadelay = abs(expected_delay - delay_lama_r);
			if (deltadelay < (double(sp.bit_time) / 2)) veroiatnosti[i].second += 1;
			Prog_bar.SetPos(j+1);
		}
		assessments[i] /= runs;
		veroiatnosti[i].second /= runs;
		veroiatnosti[i].first = noize_r;
		ProgBarRes.SetPos(i+1);
	}
	//double test_pic = 1.2345;
	//string str = convertToStrPng<double>(&test_pic);
	ViewerDraw(assessments, noize_min_r, noize_max_r, viewer3, "assessments.png");
	vectorDoubleToFile(veroiatnosti, "veroiatnosti.txt");
	SetCursor(LoadCursor(nullptr, IDC_ARROW));
	UpdateData(0);
}
