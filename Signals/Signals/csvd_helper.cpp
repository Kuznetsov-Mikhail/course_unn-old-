#pragma once
#include "pch.h"
#include "Signals_Processing.h"
#include "cubic.h"
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Signals_Processing::trans_matr(const vector<vector<complex<double>>>& v1,
	const vector<vector<complex<double>>>& v2, vector<vector<complex<double>>>& v3)
{
	vector<vector<complex<double>>> buffer;
	buffer.clear();
	buffer.resize(v1.size());
	for (int i = 0; i < buffer.size(); i++)
	{
		buffer[i].resize(v2[0].size());
		for (int j = 0; j < buffer[i].size(); j++)
		{
			for (int r = 0; r < v2.size(); r++)
				buffer[i][j] += v1[i][r] * (v2[r][j]);
		}
	}
	v3 = buffer;
}
void Signals_Processing::transpose_conj(vector<vector<complex<double>>>& v1)
{
	vector<vector<complex<double>>> buffer; buffer.resize(v1[0].size());
	for (int i = 0; i < buffer.size(); i++)
	{
		buffer[i].resize(v1.size());
		for (int j = 0; j < buffer[i].size(); j++)
		{
			buffer[i][j] = conj(v1[j][i]);
		}
	}
	v1 = buffer;
}
void Signals_Processing::vec_to_2dvec(const vector<complex<double>>& v1, vector<vector<complex<double>>>& v2)
{
	v2.clear();
	int N = sqrt(v1.size());
	v2.resize(N);
	for (int i = 0; i < N; i++)
	{
		v2[i].resize(N);
		for (int j = 0; j < N; j++)
			v2[i][j] = v1[i * N + j];
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
int Signals_Processing::CSVD(vector<complex<double>> A, int M, int N, int NU, int NV, vector<double>& S,
	vector<complex<double>>& U, vector<complex<double>>& V)
{
	int i, j, k, k1, kk, kc, kc1, l, l1, ll;
	int nmax;
	int IP = 0, NP, N1;
	double Q_re, Q_im, R_re, R_im, Z, W, tmp, tmp1, EPS, CS, SN, F, H;
	double X, Y, G;
	//double* B, * C, * T;
	vector <double> B, C, T;
	double ETA = 1.2E-7, TOL = 2.4E-16;
	int index, index1;

	nmax = N;
	if (M > nmax) nmax = M;
	//B = new double[nmax];
	//C = new double[nmax];
	//T = new double[nmax];
	B.resize(nmax);
	C.resize(nmax);
	T.resize(nmax);

	for (i = 0; i < nmax; i++)
	{
		B[i] = 0;
		C[i] = 0;
		T[i] = 0;
	}

	NP = N + IP;
	N1 = N + 1;
	//	ponivenie porqdka
	C[0] = 0;
	k = 1;
m10:
	k1 = k + 1;

	kc = k - 1;
	kc1 = k1 - 1;
	//	iskl`~enie |lementow a(I,K),I=K+1,...,M
	Z = 0;
	for (i = kc; i < M; i++) { index = i * N + kc; Z += A[index].real() * A[index].real() + A[index].imag() * A[index].imag(); }
	B[kc] = 0;
	if (Z <= TOL) goto m70;
	Z = (double)sqrt(Z);
	B[kc] = Z;
	index = kc * N + kc;
	W = (double)sqrt(A[index].real() * A[index].real() + A[index].imag() * A[index].imag());
	Q_re = 1; Q_im = 0;
	if (fabs(W) >= TOL * TOL) { Q_re = A[index].real() / W; Q_im = A[index].imag() / W; }
	A[index] = std::complex<double>((double)(Q_re * (Z + W)), (double)(Q_im * (Z + W)));

	if (k == NP) goto m70;

	for (j = kc1; j < NP; j++)
	{
		Q_re = 0; Q_im = 0;
		for (i = kc; i < M; i++)
		{
			index = i * N + j;
			index1 = i * N + kc;
			Q_re += A[index].real() * A[index1].real() + A[index].imag() * A[index1].imag();
			Q_im += A[index].imag() * A[index1].real() - A[index].real() * A[index1].imag();
		}
		Q_re /= (Z * (Z + W));
		Q_im /= (Z * (Z + W));

		for (i = kc; i < M; i++)
		{
			index = i * N + j;
			index1 = i * N + kc;
			tmp = Q_re * A[index1].real() - Q_im * A[index1].imag();
			tmp1 = Q_im * A[index1].real() + Q_re * A[index1].imag();
			A[index] -= std::complex<double>((double)tmp, (double)tmp1);
		}
	}
	//	preobrazowanie fazy
	index = kc * N + kc;
	tmp = (double)sqrt(A[index].real() * A[index].real() + A[index].imag() * A[index].imag());
	Q_re = -A[index].real() / tmp;
	Q_im = A[index].imag() / tmp;

	for (j = kc1; j < NP; j++)
	{
		index = kc * N + j;
		tmp = Q_re * A[index].real() - Q_im * A[index].imag();
		tmp1 = Q_im * A[index].real() + Q_re * A[index].imag();
		A[index] = std::complex<double>((double)tmp, (double)tmp1);
	}
	//	iskl`~enie |lementow a(K,J),J=K+2,...,N
m70:
	if (k == N) goto m140;
	Z = 0;
	for (j = kc1; j < N; j++) { index = kc * N + j; Z += A[index].real() * A[index].real() + A[index].imag() * A[index].imag(); }

	C[kc1] = 0;
	if (Z <= TOL) goto m130;
	Z = (double)sqrt(Z);
	C[kc1] = Z;
	index = kc * N + kc1;
	W = (double)sqrt(A[index].real() * A[index].real() + A[index].imag() * A[index].imag());
	Q_re = 1; Q_im = 0;
	if (fabs(W) > TOL * TOL) { Q_re = A[kc * N + kc1].real() / W; Q_im = A[kc * N + kc1].imag() / W; }
	A[index] -= std::complex<double>((double)(Q_re * (Z + W)), (double)(Q_im * (Z + W)));

	for (i = kc1; i < M; i++)
	{
		Q_re = 0; Q_im = 0;
		for (j = kc1; j < N; j++)
		{
			index = i * N + j;
			index1 = kc * N + j;
			Q_re += A[index].real() * A[index1].real() + A[index].imag() * A[index1].imag();
			Q_im += A[index].imag() * A[index1].real() - A[index].real() * A[index1].imag();
		}
		Q_re /= (Z * (Z + W));
		Q_im /= (Z * (Z + W));

		for (j = kc1; j < N; j++)
		{
			index = i * N + j;
			index1 = kc * N + j;
			tmp = Q_re * A[index1].real() - Q_im * A[index1].imag();
			tmp1 = Q_im * A[index1].real() + Q_re * A[index1].imag();
			A[index] -= std::complex<double>((double)tmp, (double)tmp1);
		}
	}

	//	preobrazowanie fazy
	index = kc * N + kc1;
	tmp = (double)sqrt(A[index].real() * A[index].real() + A[index].imag() * A[index].imag());
	Q_re = -A[index].real() / tmp;
	Q_im = A[index].imag() / tmp;

	for (i = kc1; i < M; i++)
	{
		index = i * N + kc1;
		tmp = Q_re * A[index].real() - Q_im * A[index].imag();
		tmp1 = Q_im * A[index].real() + Q_re * A[index].imag();
		A[index] = std::complex<double>((double)tmp, (double)tmp1);
	}
m130:
	k = k1;
	goto m10;
	//	dopusk dlq prenebrevimo malyh |lementow
m140:
	EPS = 0;
	for (k = 0; k < N; k++)
	{
		S[k] = (double)B[k];
		T[k] = C[k];
		if ((S[k] + T[k]) > EPS) EPS = S[k] + T[k];
	}
	EPS = (double)(EPS * ETA);
	//	inicializaciq wy~islenij U i V
	if (NU == 0) goto m180;

	for (j = 0; j < NU; j++)
	{
		for (i = 0; i < M; i++)
		{
			U[i * M + j] = 0;
		}
		U[j * M + j] = std::complex<double>(1, 0);
	}
m180:
	if (NV == 0) goto m210;
	for (j = 0; j < NV; j++)
	{
		for (i = 0; i < N; i++)
		{
			V[i * N + j] = 0;
		}
		V[j * N + j] = std::complex<double>(1, 0);
	}

	//	QR-diagonalizaciq
m210:
	for (kk = 1; kk <= N; kk++)
	{
		k = N1 - kk;
		//	prowerka na ras}eplenie
	m220:
		for (ll = 1; ll <= k; ll++)
		{
			l = k + 1 - ll;
			if ((double)fabs(T[l - 1]) <= EPS) goto m290;
			if ((double)fabs(S[l - 2]) <= EPS) goto m240; // l-1 ???
		}
		//	sokra}enie B(L)
	m240:
		CS = 0;
		SN = 1;
		l1 = l - 1;
		for (i = l; i <= k; i++)
		{
			F = SN * T[i - 1];
			T[i - 1] = CS * T[i - 1];
			if ((double)fabs(F) <= EPS) goto m290;
			H = S[i - 1];
			W = (double)sqrt(F * F + H * H);
			S[i - 1] = (double)W;
			CS = H / W;
			SN = -F / W;
			if (NU == 0) goto m260;
			for (j = 0; j < N; j++)
			{
				index = j * M + l1 - 1;
				index1 = j * M + i - 1;
				X = U[index].real();
				Y = U[index1].real();
				U[index] = std::complex<double>((double)(X * CS + Y * SN), 0);
				U[index1] = std::complex<double>((double)(Y * CS - X * SN), 0);
			}
		m260:
			if (NP == N) continue;
			for (j = N1; j <= NP; j++)
			{
				index = (l1 - 1) * N + j - 1;
				index1 = (i - 1) * N + j - 1;
				Q_re = A[index].real();
				Q_im = A[index].imag();
				R_re = A[index1].real();
				R_im = A[index1].imag();
				A[index] = std::complex<double>((double)(Q_re * CS + R_re * SN), (double)(Q_im * CS + R_im * SN));
				A[index1] = std::complex<double>((double)(R_re * CS - Q_re * SN), (double)(R_im * CS - Q_im * SN));
			}
		}

		//	prowerka shodimosti
	m290:
		W = S[k - 1];
		if (l == k) goto m360;
		//	sdwig na~ala koordinat
		X = S[l - 1];
		Y = S[k - 2];
		G = T[k - 2];
		H = T[k - 1];
		F = ((Y - W) * (Y + W) + (G - H) * (G + H)) / (2 * H * Y);
		G = (double)sqrt(F * F + 1);
		if (F < 0) G = -G;
		F = ((X - W) * (X + W) + (Y / (F + G) - H) * H) / X;
		//	QR {ag
		CS = 1;
		SN = 1;
		l1 = l + 1;
		for (i = l1; i <= k; i++)
		{
			G = T[i - 1];
			Y = S[i - 1];
			H = SN * G;
			G = CS * G;
			W = (double)sqrt(H * H + F * F);
			T[i - 2] = W;
			CS = F / W;
			SN = H / W;
			F = X * CS + G * SN;
			G = G * CS - X * SN;
			H = Y * SN;
			Y = Y * CS;
			if (NV == 0) goto m310;
			for (j = 0; j < N; j++)
			{
				index = j * N + i - 1;
				X = V[index - 1].real();
				W = V[index].real();
				V[index - 1] = std::complex<double>((double)(X * CS + W * SN), 0);
				V[index] = std::complex<double>((double)(W * CS - X * SN), 0);
			}

		m310:
			W = (double)sqrt(H * H + F * F);
			S[i - 2] = (double)W;
			CS = F / W;
			SN = H / W;
			F = CS * G + SN * Y;
			X = CS * Y - SN * G;
			if (NU == 0) goto m330;
			for (j = 0; j < N; j++)
			{
				index = j * M + i - 1;
				Y = U[index - 1].real();
				W = U[index].real();
				U[index - 1] = std::complex<double>((double)(Y * CS + W * SN), 0);
				U[index] = std::complex<double>((double)(W * CS - Y * SN), 0);
			}

		m330:
			if (NP == N) continue;
			for (j = N1; j <= NP; j++)
			{
				index = (i - 2) * N + j - 1;
				index1 = (i - 1) * N + j - 1;
				Q_re = A[index].real();
				Q_im = A[index].imag();
				R_re = A[index1].real();
				R_im = A[index1].imag();
				A[index] = std::complex<double>((double)(Q_re * CS + R_re * SN), (double)(Q_im * CS + R_im * SN));
				A[index1] = std::complex<double>((double)(R_re * CS - Q_re * SN), (double)(R_im * CS - Q_im * SN));
			}
		}
		T[l - 1] = 0;
		T[k - 1] = F;
		S[k - 1] = (double)X;
		goto m220;
		//	shodimostx
	m360:
		if (W >= 0) continue;
		S[k - 1] = (double)-W;
		if (NV == 0) continue;
		for (j = 0; j < N; j++)
		{
			index = j * N + k - 1;
			V[index] = -V[index];
		}
	}

	//	uporqdo~enie singulqrnyh ~isel
	for (k = 0; k < N; k++)
	{
		G = -1;
		j = k;
		for (i = k; i < N; i++)
		{
			if (S[i] <= G) continue;
			G = S[i];
			j = i;
		}
		if (j == k) continue;
		S[j] = S[k];
		S[k] = (double)G;
		if (NV == 0) goto m410;
		for (i = 0; i < N; i++)
		{
			index = i * N + j;
			index1 = i * N + k;
			std::swap(V[index], V[index1]);
		}
	m410:
		if (NU == 0) goto m430;
		for (i = 0; i < N; i++)
		{
			index = i * M + j;
			index1 = i * M + k;
			std::swap(U[index], U[index1]);
		}
	m430:
		if (NP == N) continue;
		for (i = N1 - 1; i < NP; i++)
		{
			index = j * N + i;
			index1 = k * N + i;
			std::swap(A[index], A[index1]);
		}
	}
	//	obratnoe preobrazowanie
	if (NU == 0) goto m510;
	for (kk = 1; kk <= N; kk++)
	{
		k = N1 - kk;
		kc = k - 1;
		if (B[kc] == 0) continue;
		index = kc * N + kc;
		tmp = (double)sqrt(A[index].real() * A[index].real() + A[index].imag() * A[index].imag());
		Q_re = -A[index].real() / tmp;
		Q_im = -A[index].imag() / tmp;
		for (j = 0; j < NU; j++)
		{
			index = kc * M + j;
			tmp = Q_re * U[index].real() - Q_im * U[index].imag();
			tmp1 = Q_im * U[index].real() + Q_re * U[index].imag();
			U[index] = std::complex<double>((double)tmp, (double)tmp1);
		}
		for (j = 0; j < NU; j++)
		{
			Q_re = 0; Q_im = 0;
			for (i = kc; i < M; i++)
			{
				index = i * M + j;
				index1 = i * N + kc;
				Q_re += (U[index].real() * A[index1].real() + U[index].imag() * A[index1].imag());
				Q_im += (U[index].imag() * A[index1].real() - U[index].real() * A[index1].imag());
			}
			index = kc * N + kc;
			tmp = B[kc] * (double)sqrt(A[index].real() * A[index].real() + A[index].imag() * A[index].imag());
			Q_re /= tmp;
			Q_im /= tmp;
			for (i = kc; i < M; i++)
			{
				index = i * M + j;
				index1 = i * N + kc;
				tmp = Q_re * A[index1].real() - Q_im * A[index1].imag();
				tmp1 = Q_im * A[index1].real() + Q_re * A[index1].imag();
				U[index] -= std::complex<double>((double)tmp, (double)tmp1);
			}
		}
	}
m510:
	if (NV == 0) goto m570;
	if (N < 1) goto m570;
	for (kk = 2; kk <= N; kk++)
	{
		k = N1 - kk;
		kc = k - 1;
		k1 = k + 1;
		kc1 = k1 - 1;
		if (C[kc1] == 0) continue;
		index = kc * N + kc1;
		tmp = (double)sqrt(A[index].real() * A[index].real() + A[index].imag() * A[index].imag());
		Q_re = -A[index].real() / tmp;
		Q_im = A[index].imag() / tmp;
		for (j = 0; j < NV; j++)
		{
			index = kc1 * N + j;
			tmp = Q_re * V[index].real() - Q_im * V[index].imag();
			tmp1 = Q_im * V[index].real() + Q_re * V[index].imag();
			V[index] = std::complex<double>((double)tmp, (double)tmp1);
		}
		for (j = 0; j < NV; j++)
		{
			Q_re = 0; Q_im = 0;
			for (i = kc1; i < N; i++)
			{
				index = i * N + j;
				index1 = kc * N + i;
				Q_re += (A[index1].real() * V[index].real() - A[index1].imag() * V[index].imag());
				Q_im += (A[index1].imag() * V[index].real() + A[index1].real() * V[index].imag());
			}
			index = kc * N + kc1;
			tmp = C[kc1] * (double)sqrt(A[index].real() * A[index].real() + A[index].imag() * A[index].imag());
			Q_re /= tmp;
			Q_im /= tmp;

			for (i = kc1; i < N; i++)
			{
				index = i * N + j;
				index1 = kc * N + i;
				tmp = Q_re * A[index1].real() + Q_im * A[index1].imag();
				tmp1 = Q_im * A[index1].real() - Q_re * A[index1].imag();
				V[index] -= std::complex<double>((double)tmp, (double)tmp1);
			}

		}
	}
m570: IP = 0;

	//delete[] B;
	//delete[] C;
	//delete[] T;
	B.clear();
	C.clear();
	T.clear();

	return IP;
}
int Signals_Processing::nonlinear_filtering(vector<complex<double>>& signal, double f0, double sampling, double bitrate)
{
	int win_size = sampling / bitrate;
	if (win_size < 0)return -1;
	win_size = step2(win_size);
	vector<complex<double>> template_signal;
	template_signal.resize(win_size);
	Buffaza = 0;
	for (int i = 0; i < template_signal.size(); i++)
	{
		Buffaza += (2 * M_PI * f0 / sampling);
		NormalPhaza(Buffaza);
		template_signal[i] = cos(Buffaza) + comjd * sin(Buffaza);
	}
	//fur(template_signal, -1);
	//for (int i = 0; i < template_signal.size(); i++)
	//	template_signal[i] = pow(abs(template_signal[i]), 2);
	//fur(template_signal, 1);
	////////////////////////////////////////////////////////////////////
	vector<complex<double>> A; A.resize(pow(template_signal.size(), 2));
	for (int i = 0; i < template_signal.size(); i++)
	{
		for (int j = 0; j < template_signal.size(); j++)
		{
			A[i * template_signal.size() + j] = template_signal[abs(i - j)];
		}
	}
	vector<double> S; S.resize(win_size);
	vector<complex<double>> U, V; U.resize(win_size * win_size); V.resize(win_size * win_size);
	int error = CSVD(A, win_size, win_size, win_size, win_size, S, U, V);
	if (error == -1)return -1;
	////////////////////////////////////////////////////////////////////
	vector<vector<complex<double>>> SS, UU, VV, AA;
	vec_to_2dvec(U, UU);
	vec_to_2dvec(V, VV);
	transpose_conj(UU);
	transpose_conj(VV);
	for (int i = 0; i < S.size(); i++) if (abs(S[i] > 0.001))S[i] = 1. / S[i];
	SS.resize(win_size);
	for (int i = 0; i < win_size; i++)
	{
		SS[i].resize(win_size);
		SS[i][i] = S[i];
	}
	trans_matr(VV, SS, AA);
	trans_matr(AA, UU, AA);
	////////////////////////////////////////////////////////////////////
	for (int i = 0; i < signal.size(); i++)
	{
		if (i < signal.size() - win_size)
		{
			vector<vector<complex<double>>> buffer, out; buffer.resize(win_size);
			for (int j = i; j < i + win_size; j++)
			{
				buffer[j - i].resize(1);
				buffer[j - i][0] = signal[j];
			}
			transpose_conj(buffer);
			trans_matr(buffer, AA, out);
			transpose_conj(buffer);
			trans_matr(out, buffer, out);
			signal[i] = out[0][0];
		}
		else signal[i] = 0;
	}
}