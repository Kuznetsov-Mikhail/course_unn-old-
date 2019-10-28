#ifdef FROMDLL
	#define FCLIBAPI extern "C" __declspec(dllexport)
#else
	#define FCLIBAPI extern "C" __declspec(dllimport)
#endif

// Filtration algorithms
enum Filtering_Algorithm
{
	CPU_TD,
	CPU_FD,
	GPU_TD,
	GPU_FD
};

// Alias for vector < complex < float > > 
typedef std::vector < std::complex < float > > signal_buf;
// Alias for vector < vector < complex < float > > >
typedef std::vector < signal_buf > bank_buf;

/** Fast_convolustion function
 *  @param signal	Input signal
 *  @param filters	Bank of filters
 *  @param outputs  Output signals
 *  @param alg		Algorithm
 */
FCLIBAPI void __stdcall fast_convolution(	const signal_buf &signal,
											const bank_buf &filters,
											bank_buf &outputs,
											Filtering_Algorithm alg);

 /**
* ��������� ������ �� �����
* @param file_name ��� �����
* @param filter    ������, � ������� ����� �������� ������������ �������
* @return          ������� ������ ��������
*/
FCLIBAPI bool __stdcall LoadFilter(    const std::string &file_name, std::vector <std::complex <float> > &filter);

/** ������������� ��������� ���-������
  * freq1				������ ������� ����� (������� ������� - ��, ���, ��� � �.�. �����, ����� ��� ������� ���� � ����� ��������)
  * freq2				������� ������� �����
  * trans_zone_width	������ ���������� ����
  * sampling_rate		������� �������������
  * fir					���������� �������������� ���-�������
*/
FCLIBAPI void __stdcall CreateFirFilter(		double freq1, 
												double freq2, 
												double trans_zone_width, 
												double sampling_rate, 
												std::vector <std::complex <float> > &fir);