#include <iostream>
#include <vector>
#include <tuple>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>

#include <functional>

using namespace std;

///////////////////////////////////////////////////////////
// ������ ������� ��� �������� ������� ������
auto uExact = [](double t, double x) -> double {
	return exp(-(M_PI / 3.0) * (M_PI / 3.0) * t) * sin(M_PI * x / 3.0);
};


// ���������� ������
double absErr = 0.0;

// ��� K(u)

function<double(double, double)> uExactU;
///////////////////////////////////////////////////////////


// ����� �������
void print_vec(const vector<vector<double>>& res)
{

	for (size_t i = 0; i < res.size(); i++) {
		for (size_t j = 0; j < res[i].size(); j++)
			cout << res[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}

void print_vec(const vector<double>& res)
{
	for (size_t i = 0; i < res.size(); i++)
		cout << fixed << res[i] << "\t";
	cout << endl;
}


// ����� �������� ���������� ������� � �������
// vec - ��������� �������
// eExact(t, x)
// t - ������� ��������� ����
// L - ����� �������
double dist(const vector<double>& vec, function<double(double, double)> uExact,
	double t, double L) {
	size_t n = vec.size() - 1;
	double res = 0.0;
	for (size_t i = 0; i <= n; ++i) {
		double x = L * i / n;
		res = max(res, abs(vec[i] - uExact(t, x)));
	}
	return res;
}

// ����� ��������
template<typename T>
vector<T> ThomasAlgorithm(
	const vector<T>& a, const vector<T>& b, const vector<T>& c,
	const vector<T>& d
) {
	size_t n = b.size();
	vector<T> alpha, beta;
	alpha.reserve(n - 1);
	beta.reserve(n - 1);

	// i = 1
	alpha.push_back(-c[0] / b[0]);
	beta.push_back(d[0] / b[0]);
	// i = 2..n-1
	for (size_t i = 1; i < n - 1; ++i) {
		T denom = -b[i] - a[i - 1] * alpha.back();
		alpha.push_back(c[i] / denom);
		beta.push_back((-d[i] + a[i - 1] * beta.back()) / denom);
	}
	vector<T> X;
	X.reserve(n);
	X.push_back((-d[n - 1] + a[n - 2] * beta[n - 2]) / (-b[n - 1] - a[n - 2] * alpha[n - 2]));
	for (size_t i = n - 2; i >= 0 && i < n; --i)
		X.push_back(alpha[i] * X.back() + beta[i]);

	reverse(X.begin(), X.end());
	return X;
}

// ������ � ���� �������
template<typename T>
ostream& operator<<(ostream& out, vector<T> v) {
	for (size_t i = 0; i < v.size() - 1; ++i)
		out << v[i] << "    ";
	out << v.back();
	return out;
}

// ������� ��������� ����������������
// crho - ������������ ������������� �������� ������������ � �������� ���������
// K(x) - ����������� ����������������
// uInit(x) - ��������� ������� (��� t = 0)
// uL(t) - ����������� �� ����� ����� (��� x = 0)
// PR(t) - ����� �� ������ ����� (��� x = L)
// sigma - �������� ����� (0 <= sigma <= 1)
// ���� condL = false, �� ����� ������� ���� u(0, t) = f1(t),
// ����� - ���� -K(0) du/dx (0, t) = f1(t)
// ���� condR = false, �� ������ ������� ���� u(L, t) = f2(t),
// ����� - ���� K(L) du/dx (L, t) = f2(t)
template<bool condL = false, bool condR = false>
void solver(double crho, function<double(double)> K, function<double(double)> uInit,
	double T, double L, double sigma, double tau, double h,
	function<double(double)> fL, function<double(double)> fR, ostream& out) {
	////////////////////////////
	// �������� ����������������
	ofstream cons("cons.txt");
	////////////////////////////


	size_t n = round(L / h);
	vector<double> xx;
	xx.reserve(n + 1);
	for (size_t i = 0; i <= n; ++i)
		xx.push_back(i * h);

	vector<double> u;
	u.reserve(n + 1);
	// �������������� ������� ���� ��������� ��������
	for (size_t i = 0; i <= n; ++i)
		u.push_back(uInit(xx[i]));

	vector<double> a;
	a.reserve(n);
	for (size_t i = 0; i < n; ++i)
		// ��������� ������ �������
		a.push_back((K(xx[i + 1]) + K(xx[i])) / 2);

	// ����� � C ����� �������� �� -1, ������ ��� � ��� � ������ �������� ���
	// ������������� � ���� �������������
	vector<double> A, B, C;
	A.reserve(n);
	B.reserve(n);
	C.reserve(n + 1);

	if (!condL) {
		C.push_back(1.0);
		B.push_back(0.0);
	}
	else {
		// TODO �� ��������
		double kappa = sigma * a[0] / h /
			(crho * h / (2 * tau) + sigma * a[0] / h);

		C.push_back(-1.0);
		B.push_back(kappa);
	}

	for (size_t i = 1; i <= n - 1; ++i) {
		A.push_back(sigma * a[i - 1] / h);
		B.push_back(sigma * a[i] / h);
		C.push_back(-1 * (sigma * (a[i - 1] + a[i]) / h + crho * h / tau));
	}

	if (!condR) {
		C.push_back(1.0);
		A.push_back(0.0);
	}
	else {
		double kappa = sigma * a[n - 1] / h /
			(crho * h / (2 * tau) + sigma * a[n - 1] / h);

		A.push_back(kappa);
		C.push_back(-1.0);
	}

	//cout << n << "    " << A.size() << "    " << C.size() << "    " << B.size() << endl;

	vector<double> w;
	w.reserve(n);

	// ��������� w ��� �������� ����
	for (size_t i = 0; i < n; ++i)
		w.push_back(a[i] * (u[i + 1] - u[i]) / h);

	// ���������� � ���� ������� ����
	out << u << endl;

	vector<double> F(n + 1, 0);

	size_t m = round(T / tau);
	for (size_t j = 1; j <= m; ++j) {
		double t = j * tau;

		if (!condL) {
			F[0] = fL(t);
		}
		else {
			double mu = (crho * u[0] * h / (2 * tau) + sigma * fL(t) + (1 - sigma) * (fL(t - tau) + w[0])) /
				(crho * h / (2 * tau) + sigma * a[0] / h);

			F[0] = -mu;
		}

		// ������ ����� ����� ��������� �� -1
		for (size_t i = 1; i < n; ++i) {
			F[i] = crho * h / tau * u[i] + (1 - sigma) * (w[i] - w[i - 1]);
			F[i] *= -1;
		}

		if (!condR) {
			F[n] = fR(t);
		}
		else {
			double mu = (crho * u[n] * h / (2 * tau) + sigma * fR(t) + (1 - sigma) * (fR(t - tau) - w[n - 1])) /
				(crho * h / (2 * tau) + sigma * a[n - 1] / h);
			F[n] = -mu;
		}

		// � ��������� � ���� ������������ ���� � ������� A, C, B
		//u = ThomasAlgorithm(A, C, B, F);

		////////////////////////////////////////////////
		// �������� ������ ���������� (����������������)
		vector<double> uNew = ThomasAlgorithm(A, C, B, F);
		auto sumOld = 0.5 * u[0] + 0.5 * u[n] + accumulate(u.begin() + 1, u.end() - 1, 0.0);
		auto sumNew = 0.5 * uNew[0] + 0.5 * uNew[n] + accumulate(uNew.begin() + 1, uNew.end() - 1, 0.0);
		//cout << sumOld << "    " << sumNew << "    " << crho * h * abs(sumNew - sumOld) << endl;
		cons << abs(sumNew - sumOld) << endl;
		u = uNew;
		////////////////////////////////////////////////

		// ������������� w ��� �������� ����
		for (size_t i = 0; i < n; ++i)
			w[i] = a[i] * (u[i + 1] - u[i]) / h;

		// ���������� � ����
		out << u << endl;
	}
	////////////////////////////
	// �������� ����������������
	cons.close();
	////////////////////////////

}




// ������� ��� 5 �������
void quest5() {
	double L = 1.0;
	double T = 1.0;
	//double T = 1000.0;

	double crho = 1.0;
	auto K = [](double x) -> double {return 1.0; };
	double sigma = 0.5;

	/*auto uL = [](double t) -> double {return 5.0 - t / 10.0; };
	//auto uR = [](double t) -> double {return 5.0 + t; };
	auto uR = [T](double t) -> double {return 5.0 - t; };

	//auto uLR = [](double t) -> double {return 5.0; };
	auto uLR = [T](double t) -> double {return 5.0 + t * (T - t); };*/

	auto PLR = [T](double t) -> double {return 0.0; };

	//auto uInit = [](double x) -> double {return 5.0; };
	//auto uInit = [L](double x) -> double {return 5.0 + x * (L - x); };
	auto uInit = [L](double x) -> double {return 5.0 - x * (L - x); };

	// ���� ������������ (h^2 > tau)
	//double h = 0.1;
	//double h = 0.5;
	//double tau = 0.009;

	// ��� ������������ (h^2 <= 2/3 * tau, sigma = )
	double h = 0.1;
	double tau = 1;

	ofstream out("test.txt");
	// ���������� ��������� ��� ���������� ��������
	out << L << "    " << T << endl;

	/*solver<false, false>(crho, K, uInit, T, L,
		sigma, tau, h,
		//uL, uR, out
		uLR, uLR, out
	);*/

	solver<true, true>(crho, K, uInit, T, L,
		sigma, tau, h,
		//uL, uR, out
		PLR, PLR, out
		);

	out.close();

}

void quest5_2() {
	double L = 10.0;
	double T = 10.0;

	double crho = 1.0;
	auto K = [](double x) -> double { return 1.0; };
	double sigma = 0.0;

	//auto uInit = [L](double x) -> double { return 5.0 - x * (L - x); };
	auto uInit = [L](double x) -> double { return 5.0 + x * (L - x); };
	auto uLR = [](double t) -> double { return 5.0; };

	// ������������ ������
	double h = 0.25;
	double tau = 15 * h * h / 20;

	// ���������� ������
	/*double h = 0.1;
	double tau = 0.005;*/

	ofstream out("test.txt");
	// ���������� ��������� ��� ���������� ��������
	out << L << "    " << T << endl;

	solver<false, false>(crho, K, uInit, T, L,
		sigma, tau, h,
		uLR, uLR, out
		//PLR, PLR, out
		);

	out.close();

}


std::fstream& GotoLine(std::fstream& file, unsigned int num) {
	file.seekg(std::ios::beg);
	for (int i = 0; i < num - 1; ++i) {
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	return file;
}

// ���������� �������
// sigma - ���
// t - ��������� ����, � �������� ����� ������
// q - �������� �����
void order_zvezda() {

	auto K = [](double x) -> double { return 1.0; };
	auto uLR = [](double t) -> double { return 0.0; };
	auto uInit = [](double x) -> double { return sin(M_PI * x / 3.0); };
	double L = 3.0;
	double T = 1.0;
	double t_zvezda = T / 10.0;
	double err;
	double err_max = 0.0;
	double h = 0.1;
	double tau = 0.01;
	double q = 1.0 / 2.0;
	double sigma = 0.0;
	double value;
	double cur = 0;
	int w = 14;
	vector <vector<double>> res;
	vector<double> str;

	cout << left << setw(w) << "tau" << setw(w) << "h" << setw(w) << "AbsErr(tau)" << setw(w) << "Delta" << setw(w) << "log(Delta)" << endl;


	for (size_t i = 0; i < 6; i++) {
		res.clear();

		//�������� �������
		ofstream out("order.txt");
		out.clear();
		solver<false, false>(1.0, K, uInit, T, L,
			sigma, tau, h,
			uLR, uLR, out
			);
		out.close();


		fstream file("order.txt"); // ��������� ���� ��� ������
		if (file.is_open()) {
			//���������� ����� � ����� (��������� �����)
			/*int count = std::count(std::istreambuf_iterator<char>(file),
				std::istreambuf_iterator<char>(), '\n');*/

			for (size_t j = 1; j <= round(t_zvezda / tau); j++) {
				/*cout << "��������� ����" << endl;
				cout << "���������� �����: " << round(t_zvezda / tau) << endl;*/
				// ��������� � ������� ����
				//cout << round(t_zvezda / tau) << endl;
				GotoLine(file, j);
				for (size_t jj = 0; jj <= round(L / h); jj++) {
					file >> value;
					str.push_back(value);
				}
				/*cout << "����" << endl;
				print_vec(str);*/
				res.push_back(str);
				str.clear();
			}
		}
		file.close();
		/*cout << "������ �������� ���������" << endl;
		print_vec(res);*/

		//double err = dist(res, uExactU, t_zvezda, L);
		size_t k = res[0].size() - 1;
		for (size_t jj = 0; jj < round(t_zvezda / tau); jj++) {
			err = 0.0;
			for (size_t i = 0; i <= k; ++i) {
				double x = L * i / k;
				err = max(err, abs(res[jj][i] - uExact(jj * tau, x)));
				//cout << "����������� " << res[jj][i] - uExact(jj * tau, x) << endl;
			}
			err_max = max(err_max, err);
		}
		if (i == 0)
			cout << left << setw(w) << tau << setw(w) << h << setw(w) << err_max << setw(w) << "---" << setw(w) << "---" << endl;
		else
			cout << left << setw(w) << tau << setw(w) << h << setw(w) << err_max << setw(w) << (err_max / cur) << setw(w) << log2(err_max / cur) / log2(q) << endl;
		cur = err_max;
		err_max = 0.0;

		h = q * h;
		tau = q * q * tau;

		/*h = q * h;
		tau = q * tau;*/
	}
	//print_vec(res);
}


// ���������� �������
// sigma - ���
// t - ��������� ����, � �������� ����� ������
// q - �������� �����
void order() {

	auto K = [](double x) -> double { return 1.0; };
	auto uLR = [](double t) -> double { return 0.0; };
	auto uInit = [](double x) -> double { return sin(M_PI * x / 3.0); };
	double L = 3.0;
	double T = 1.0;
	double err;
	double h = 0.5;
	double tau = 0.05;
	double q = 1.0 / 2.0;
	double sigma = 0;
	double value;
	double cur = 0;
	int w = 14;
	vector<double> res;

	cout << left << setw(w) << "tau" << setw(w) << "h" << setw(w) << "AbsErr(tau)" << setw(w) << "Delta" << setw(w) << "log(Delta)" << endl;


	for (size_t i = 0; i < 6; i++) {
		double t = round(T / tau/2);
		//cout << t << endl;
		res.clear();
		//�������� �������
		ofstream out("order.txt");
		out.clear();
		solver<false, false>(1.0, K, uInit, T, L,
			sigma, tau, h,
			uLR, uLR, out
			);
		out.close();


		fstream file("order.txt"); // ��������� ���� ��� ������
		if (file.is_open()) {
			//���������� ����� � ����� (��������� �����)
			/*int count = std::count(std::istreambuf_iterator<char>(file),
				std::istreambuf_iterator<char>(), '\n');*/

			GotoLine(file, t);
			for (size_t j = 0; j <= round(L / h); j++) {
				file >> value;
				res.push_back(value);
			}
			//print_vec(res);
		}
		file.close();
		/*cout << "������ �������� ���������" << endl;
		print_vec(res);*/

		size_t k = res.size() - 1;
		err = 0.0;
		for (size_t i = 0; i <= k; ++i) {
			double x = L * i / k;
			err = max(err, abs(res[i] - uExact((t-1)  * tau, x)));
			/*cout << res[i] << endl;
			cout << uExact((t-1) * tau, x) << endl;*/
		}
		if (i == 0)
			cout << left << setw(w) << tau << setw(w) << h << setw(w) << err << setw(w) << "---" << setw(w) << "---" << endl;
		else
			cout << left << setw(w) << tau << setw(w) << h << setw(w) << err << setw(w) << (err / cur) << setw(w) << log2(err / cur) / log2(q) << endl;
		cur = err;

		h = q * h;
		tau = q * q * tau;

		/*h = q * h;
		tau = q * tau;*/
	}
	//print_vec(res);
}

void test25() {

	double L = 1.0;
	double T = 10.0;
	double a = 3.0;
	double A = 5.0;
	auto K = [a](double x) -> double {return a; };
	auto uInit = [A, L](double x) -> double {return A * x / L; };
	auto uL = [](double t) -> double {return 0; };
	auto uR = [A](double t) -> double {return A * exp(-t); };

	ofstream out("test.txt");
	// ���������� ��������� ��� ���������� ��������
	out << L << "    " << T << endl;

	solver<false, false>(1.0, K, uInit, T, L,
		0.5, 0.01, 0.01,
		uL, uR, out
		);

	out.close();
}

// ���� ��� �������� ���������������� (��� ����� ����������������)
void testCons() {
	double L = 20.0;
	double T = 1.0;
	double a = 3.0;
	auto K = [a](double x) -> double {return a; };
	//auto uInit = [A, L](double x) -> double {return A * x / L; };
	auto uInit = [L](double x) -> double {return 10.0 + x * (L - x); };
	// ������������������ �����
	auto Plr = [](double t) -> double { return 0; };

	ofstream out("testCons.txt");
	// ���������� ��������� ��� ���������� ��������
	out << L << "    " << T << endl;

	// �������� �� ��������� � ����������� �� ������ ���������� �����?
	double sigma = 1.0;

	solver<true, true>(1.0, K, uInit, T, L,
		sigma, 0.005, 0.55,
		Plr, Plr, out
		);

	out.close();

}

void test1() {
	double crho = 2.0;
	double u0 = 15.0;
	double L = 1.0;
	double T = 10.0;

	//auto K = [](double x) -> double {return x; };
	// � ������ ������ � ������ ������ ������� ������ ����������� ���������
	// ������������ L/2
	//auto K = [L](double x) -> double {return sin(x / L * M_PI); };
	//auto K = [L](double x) -> double {return sin(x) + cos(x); };

	//auto uInit = [A, L](double x) -> double {return u0 + x * (L-x); };
	auto uLR = [u0](double t) -> double {return u0; };
	//auto uLR = [u0](double t) -> double {return u0 + t; };

	// ������ 2.5
	double a = 3.0;
	double A = 5.0;
	auto K = [a](double x) -> double {return a; };
	auto uInit = [A, L](double x) -> double {return A * x / L; };
	auto uL = [](double t) -> double {return 0; };
	auto uR = [A](double t) -> double {return A * exp(-t); };

	ofstream out("test.txt");
	// ���������� ��������� ��� ���������� ��������
	out << L << "    " << T << endl;

	/*solver<false, false>(crho, K, uInit, T, L,
		0.5, 0.01, 0.01,
		uLR, uLR, out
	);*/

	solver<false, false>(crho, K, uInit, T, L,
		0.5, 0.01, 0.01,
		uL, uR, out
		);

	/*solver<true, true>(crho, K, uInit, T, L,
		0.5, 0.01, 0.01,
		[](double t) -> double {return 0; }, [](double t) -> double {return 0; }, out
	);*/

	/*solver<true, true>(crho, K, uInit, T, L,
		0.5, 0.01, 0.01,
		[](double t) -> double {return 0; }, [](double t) -> double {return 0; }, out
	);*/

	out.close();
}

// ������� ��������� ���������������� � �������� ������������ �� ����� �������
// � ������� �� ������ � ������ K(u)
// crho - ������������ ������������� �������� ������������ � �������� ���������
// K(U) - ����������� ����������������
// uInit(x) - ��������� ������� (��� t = 0)
// M - ���������� �������� �� ����������� ������ ������� �������� (M > 0)
// ���� condL = false, �� ����� ������� ���� u(0, t) = f1(t),
// ����� - ���� -K(0) du/dx (0, t) = f1(t)
// ���� condR = false, �� ������ ������� ���� u(L, t) = f2(t),
// ����� - ���� K(L) du/dx (L, t) = f2(t)
template<bool condL = false, bool condR = false>
void solverU(double crho, function<double(double)> K, function<double(double)> uInit,
	double T, double L, size_t M, double tau, double h,
	function<double(double)> fL, function<double(double)> fR, ostream& out) {
	///////////////////////////
	// ��� ������� ������������
	ofstream outErr("error.txt");
	// ���������� ��� ��� �����������
	// ������������ ����� �� ��� X
	outErr << tau << endl;
	///////////////////////////

	size_t n = round(L / h);
	vector<double> xx;
	xx.reserve(n + 1);
	for (size_t i = 0; i <= n; ++i)
		xx.push_back(i * h);

	vector<double> u;
	u.reserve(n + 1);
	// �������������� ������� ���� ��������� ��������
	for (size_t i = 0; i <= n; ++i)
		u.push_back(uInit(xx[i]));

	// ���������� � ���� ������� ����
	out << u << endl;

	// �������� ������ ��� ������������ ����
	vector<double> A(n, 0);
	vector<double> C(n + 1, 0);
	vector<double> B(n, 0);
	vector<double> F(n + 1, 0);

	vector<double> kk(n + 1, 0);
	vector<double> aa(n, 0);

	size_t m = round(T / tau);
	for (size_t j = 1; j <= m; ++j) {
		double t = j * tau;

		vector<double> uNew = u;

		for (size_t k = 1; k <= n - 1; ++k)
			F[k] = crho * h / tau * u[k];

		for (size_t s = 0; s < M; ++s) {

			// ������� ����
			if (!condL) {
				C[0] = 1.0;
				B[0] = 0.0;
				F[0] = fL(t);
			}
			else {
				// TODO �� ��������
				double kappa = aa[0] / h /
					(crho * h / (2 * tau) + aa[0] / h);

				double mu = (crho * h / (2 * tau) * u[0] + fL(t)) /
					(crho * h / (2 * tau) + aa[0] / h);

				B[0] = -kappa;
				C[0] = 1.0;
				F[0] = mu;
			}

			kk[0] = K(uNew[0]);
			for (size_t k = 1; k <= n; ++k) {
				kk[k] = K(uNew[k]);
				aa[k - 1] = (kk[k] + kk[k - 1]) / 2.0;
			}

			if (!condR) {
				A[n - 1] = 0.0;
				C[n] = 1.0;
				F[n] = fR(t);
			}
			else {
				double kappa = aa[n - 1] / h /
					(crho * h / (2 * tau) + aa[n - 1] / h);

				double mu = (crho * h / (2 * tau) * u[n] + fR(t)) /
					(crho * h / (2 * tau) + aa[n - 1] / h);

				A[n - 1] = -kappa;
				C[n] = 1.0;
				F[n] = mu;
			}

			for (size_t k = 1; k <= n - 1; ++k) {
				A[k - 1] = -aa[k - 1] / h;
				C[k] = crho * h / tau + (aa[k] + aa[k - 1]) / h;
				B[k] = -aa[k] / h;

				F[k] = crho * h / tau * u[k];
			}

			// ������ ����
			uNew = ThomasAlgorithm(A, C, B, F);
		}

		u = uNew;

		// ���������� � ����
		out << u << endl;

		// ������ �����������:
		double err = dist(u, uExactU, t, L);
		// TODO ����� ���������� ����� � ����
		outErr << err << endl;
	}

	///////////////////////////
	// ��� ������� ������������
	outErr.close();
	///////////////////////////
}

void testU() {
	ofstream out("testU.txt");

	double sigma = 2.0;
	double kappa0 = 0.5;
	double c = 5.0;

	double L = 10.0;
	double h = 0.2;
	//double h = 0.05;
	double T = 1.0;
	double tau = 2e-4;
	//double tau = 1e-4;

	auto uInit = [](double x) -> double {return 0.0; };
	double u0 = pow(sigma * c * c / kappa0, 1.0 / sigma);
	auto uL = [u0, sigma](double t) -> double {return u0 * pow(t, 1.0 / sigma); };
	auto PRight = [](double t) -> double { return 0.0; };

	auto K = [kappa0, sigma](double u) -> double {return kappa0 * pow(u, sigma); };

	uExactU = [sigma, c, kappa0](double t, double x) -> double {
		if (x >= c * t)
			return 0.0;
		else
			return pow(sigma * c / kappa0 * (c * t - x), 1.0 / sigma);
	};

	// TODO ���������� ��������� L � T

	solverU<false, true>(1.0, K, uInit, T, L, 3, tau, h, uL, PRight, out);

	out.close();


}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "�� ��� �� �����!" << endl;

	//test1();
	/*test25();*/
	quest5();
	//quest5_2();

	/*testCons();
	testU();*/

	//order_zvezda();
	//order();

	cout << "�� ��� � �� �����!" << endl;
}