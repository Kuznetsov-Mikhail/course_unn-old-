#include <vector>
#include <complex>
#include <chrono>
#include <omp.h>
#include "fast_convolution.h"
//#include "cubic.h"
#define M_PI 3.1415926535
#define Comj complex<double>(0,1)


using namespace std;
using namespace std::chrono;

class Signals_Processing
{

public:
	Signals_Processing();
	virtual ~Signals_Processing();

	double Buffaza; // ���������� ��� �������� ����
	double V; //������� �������
	double delataW; // ������������ �������
	double sampling;// ������� �������������
	double BrV; //bitrate
	int bit_time;
	//������� ���������� ���� �� +-2M_PI
	void NormalPhaza(double& Phaza);
	//������� ���������� �������� �� ����� ������� ������� (no complex! no FHSS!)
	void SignalFill(vector<double>& mass, vector<bool> data, bool type);
	//������� ���������� �������� �� ����� ������� ������� (complex! no FHSS!)
	void SignalFill(vector<complex<double>>& mass, vector<bool> data, bool type);
	//������� ���������� ������ ��� �������� �� ����� ������� ������� (��� SignalFill)
	void DataFill(vector <bool>& Data1, vector <bool>& Data2, int Signal_size, int Data_size, int delay);
	//���������� ����
	//NoizeV - �������� � ��
	void addNoize(vector<double>& mass, double NoizeV);
	void addNoize(vector<complex<double>>& mass, double NoizeV);
	//������� �������������
	int Correlation(vector<double>& mass, vector<double> Signal1, vector<double> Siganl2);
	int Correlation(vector<double>& mass, vector<complex<double>> Signal1, vector<complex<double>> Siganl2);
	//������� ���������������
	void Uncertainty(vector<double>& mass, vector<complex<double>> Signal1, vector<complex<double>> Siganl2, int ksum);
	void Uncertainty(vector<float>& mass, signal_buf &Signal1, signal_buf &Siganl2, int ksum);
	void Uncertainty_omp(vector<double>& mass, vector<complex<double>> Signal1, vector<complex<double>> Siganl2, int ksum);
	void Uncertainty_omp(vector<float>& mass, signal_buf &Signal1, signal_buf &Signal2, int ksum);
	/**
	������� ������������� ������
	*param PhiDopler - ������� ������
	*param samp - ������� ������������� �������
	*/
	void Dopler_shift(vector<complex<double>>& mass, double PhiDopler);
	/**
	������� ������������� ���������������
	*param Signal
	*param speed - �������� �������
	*/
	void Dopler_scaling(vector <complex<double>>& Signal, double koeff);
	void InterSpline(vector<double> &Signal, vector<double>& NewSignal, double step);//������������ ��������
	void Linear_interpolation(vector<double>& Old_Data, vector<double>& New_Data, double step);//�������� ������������
	void Cubic_Inter_spline(vector<double>& Old_Data, vector<double>& New_Data, double step);
	void fur(vector <complex<double>>& data, int is);//������ �����
	void FAST_FUR(vector <complex<double>> Signal, vector <complex<double>>& Spectr,bool is);
	void spVertex(vector <complex<double>>& Spectr);
	double Max(vector<double> Mass);
	double Min(vector<double> Mass);	

	double spline(vector<double>&massx, vector<double>&massy, vector<double>&mit, int i, double h, double x) //���������������� �������
	{
		double s = ((massx[i + 1] - x) * (massx[i + 1] - x) * (2 * (x - massx[i]) + h) * massy[i]) / (h * h * h) + ((x - massx[i]) *
			(x - massx[i]) * (2 * (massx[i + 1] - x) + h) * massy[i + 1]) / (h * h * h) + ((massx[i + 1] - x) * (massx[i + 1] - x) *
			(x - massx[i]) * mit[i]) / (h * h) + ((x - massx[i]) * (x - massx[i]) * (x - massx[i + 1]) * mit[i + 1]) / (h * h);
		return s;
	}
	int step2(int sizein);//���������� ��������� �������� pow(2,n)
	/**
	������� ��������� ���� �������� ����+MSK(�������� �������)
	������� ������� - 51 ������� JTIDS
	����� ��������� ���� - 100
	������� ������������� 350���
	*param Signal
	*param speed - �������� �������
	*/
	void FHSS(vector <complex<double>> & Signal1, vector <complex<double>> & Signal2, int signal1Size, int delaySize);
	void FHSS(int signal1Size, int delaySize); //��������� 50 FHSS ��������
	/**
	������� ��������� ���� �������� ����+MSK(�������� �������) ������������ � ������� ����� JTIDS!!!
	������� ������� - 51 ������� JTIDS
	����� ��������� ���� - 100
	������� ������������� 350���
	*param Signal
	*param speed - �������� �������
	*/
	void Link16_Signals_Generator(vector <complex<double>>& Signal1, vector <complex<double>>& Signal2, int signal1Size, int delaySize,bool scramble);
	void Simple_Signals_Generator(vector <complex<double>>& Signal1, vector <complex<double>>& Signal2, int signal1Size, int delaySize);
	/**
	*���������� ��� ������� ��������������� ��� �������� � ������� ��������
	*/
	vector <vector<complex<double>>> FHSS_Signals;//50 �������� �����������
	vector <vector<complex<double>>> FHSS_Signals_initial;//50 �������� �������
	bank_buf FHSS_Signals_fl;
	bank_buf FHSS_Signals_initial_fl;
	bank_buf fir_s;

	int pseudo_bit(int start, int end);//����� ������ ���������� ������ �������
	const vector<int> operating_frequencies //������� JTIDS � ���
	{
		969,972,975,978,981,984,987,990,993,999,1002,1005,1008,1053,1056,1059,
		1062,1065,1113,1116,1119,1122,1125,1128,1131,1134,1137,1140,1143,1146,1149,1152,1155,
		1158,1161,1164,1167,1170,1173,1176,1179,1182,1185,1188,1191,1194,1197,1200,1203,1206
	};
	double peak_intensity(vector<double> mas);	//������������ ��������� �� ������������� ��������
};

