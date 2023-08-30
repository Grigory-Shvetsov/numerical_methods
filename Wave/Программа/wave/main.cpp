#include <iostream>
#include <vector>
#include <tuple>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include <functional>

using namespace std;


// ������ � ���� �������
template<typename T>
ostream& operator<<(ostream& out, vector<T> v) {
	for (size_t i = 0; i < v.size() - 1; ++i)
		out << v[i] << "    ";
	out << v.back();
	return out;
}

// a2 - ������� ������������ a (a^2), ������ ���� ������ ����
// f(x) = u(0, x)
// g(x) = u_t (0, x)
// phi(t) = u(t, 0)
// psi(t) = u(t, L)
// fxx - ������ ����������� ������� f
// L - �����
// T - ����� ��������������
// h - ��� �� ������������
// tau - ��� �� �������
// out - ���� ���������� ���������
void solver(double a2, function<double(double)> f, function<double(double)> g,
	function<double(double)> phi, function<double(double)> psi,
	function<double(double)> fxx, double a, double b, double T, double h, double tau,
	ostream& out) {

	size_t n = round(abs(b - a) / h);
	vector<double> xx;
	xx.reserve(n + 1);
	for (size_t i = 0; i <= n; ++i)
		xx.push_back(a + i * h);

	vector<double> u0;
	u0.reserve(n + 1);
	// ��������� ������� ����
	for (size_t i = 0; i <= n; ++i)
		u0.push_back(f(xx[i]));

	// ���������� ���
	out << u0 << endl;

	double t = tau;

	vector<double> u1;
	u1.reserve(n + 1);
	// ������� ������ ����
	u1.push_back(phi(t));
	for (size_t i = 1; i <= n - 1; ++i)
		u1.push_back(u0[i] + tau * g(xx[i]) + a2 * tau * tau / 2.0 * fxx(xx[i]));
	u1.push_back(psi(t));

	// ����� ���������� ���
	out << u1 << endl;

	// �������� ������ ��� ��������� ����
	vector<double> u2(n + 1, 0);

	size_t m = round(T / tau);

	// �����������, ������� ����� �����������
	double c = a2 * tau * tau / h / h;

	if (c > 1)
		cout << "����� ����� ���� ������������!!!" << endl;

	// ������ ���� ��� ��������, ������� �������� � 2
	for (size_t j = 2; j <= m; ++j) {
		t = j * tau;

		// ��� ���������� ��������� �������
		u2[0] = phi(t);
		u2[n] = psi(t);

		for (size_t i = 1; i <= n - 1; ++i) {
			u2[i] = 2.0 * u1[i] - u0[i] + c * (u1[i - 1] - 2.0 * u1[i] + u1[i + 1]);
		}

		// ���������� � ����� ����
		out << u2 << endl;
		u0 = u1;
		u1 = u2;
	}
}

// �� �� �����, �� � ������� ����������� �� �� ���������� ������
// TODO �� ��������, �� ������
void solver(double a2, function<double(double)> f, function<double(double)> g,
	function<double(double)> phi, function<double(double)> psi, double a, double b, double T, double h, double tau,
	ostream& out) {

	// ���������� �����������
	auto fxx = [f, h](double x) -> double {
		return (f(x - h) - 2.0 * f(x) + f(x + h)) / h / h;
	};
	solver(a2, f, g, phi, psi, fxx, a, b, T, h, tau, out);

}




void test1() {
	double a2 = 1.0;
	double T = 10.0;
	double a = 0.0;
	double b = 1.0;

	auto f = [](double x) -> double { return sin(M_PI * x); };
	auto g = [](double x) -> double { return 0.0; };
	auto phi = [](double x) -> double { return 0.0; };
	auto psi = [](double x) -> double { return 0.0; };

	auto fxx = [](double x) -> double { return -M_PI * M_PI * sin(M_PI * x); };
	double h = 0.02;
	double tau = 0.01;

	ofstream out("test1.txt");
	// ���������� ��������� ��� ���������� ��������
	out << abs(b - a) << "    " << T << endl;

	solver(a2, f, g, phi, psi, fxx, a, b, T, h, tau, out);

	out.close();
}

void test2() {
	double a2 = 1.0;
	double T = 10.0;
	double a = 0.0;
	double b = 1.0;

	auto f = [](double x) -> double { return x * (1 - x); };
	auto g = [](double x) -> double { return 0.0; };
	auto phi = [](double x) -> double { return 0.0; };
	auto psi = [](double x) -> double { return 0.0; };

	auto fxx = [](double x) -> double { return -2.0; };
	double h = 0.02;
	double tau = 0.01;

	//double h = 0.005;
	//double tau = 0.0005;

	ofstream out("test2.txt");
	// ���������� ��������� ��� ���������� ��������
	out << abs(b - a) << "    " << T << endl;

	solver(a2, f, g, phi, psi, fxx, a, b, T, h, tau, out);

	out.close();
}

void task1() {
	double a2 = 1.0;
	double T = 1.0;
	double a = -2;
	double b = 2;
	auto f = [](double x) -> double { if (x >= -1.0 / 3.0 && x <= 1.0 / 3.0) return 1.0;
	else return 0.0; };
	auto g = [](double x) -> double { return 0.0; };
	auto phi = [](double x) -> double { return 0.0;};
	auto psi = [](double x) -> double { return 0.0; };
	auto fxx = [](double x) -> double { return 0.0; };

	// gamma=0.1
	double h = 0.1;
	double tau = 0.01;

	// gamma=0.5
	/*double h = 0.02;
	double tau = 0.01;*/

	//// gamma=0.75	
	/*double h = 0.15;
	double tau = 0.1125;*/

	//// gamma=1.0
	/*double h = 0.01;
	double tau = 0.01;*/
	ofstream out1("task1_analytical.txt");
	// ���������� ��������� ��� ���������� ��������
	out1 << abs(b - a) << "    " << T << endl;
	solver(a2, f, g, phi, psi,fxx, a, b, T, h, tau, out1);
	ofstream out2("task1_numerical.txt");
	// ���������� ��������� ��� ���������� ��������
	out2 << abs(b - a) << "    " << T << endl;
	solver(a2, f, g, phi, psi, a, b, T, h, tau, out2);
}

void task2() {
	double a2 = 1.0;
	double T = 1.0;
	double a = -1;
	double b = 1;
	auto f = [](double x) -> double {return 0.0;};
	auto g = [](double x) -> double { if (x >= -1.0 / 2.0 && x <= 1.0 / 2.0) return 1-2*abs(x);
	else return 0.0; };
	auto phi = [](double x) -> double { return 0.0;};
	auto psi = [](double x) -> double { return 0.0; };
	auto fxx = [](double x) -> double { return 0.0; };

	// gamma=0.1
	/*double h = 0.1;
	double tau = 0.01;*/

	// gamma=0.5
	/*double h = 0.02;
	double tau = 0.01;*/

	//// gamma=0.75	
	/*double h = 0.15;
	double tau = 0.1125;*/

	//// gamma=1.0
	double h = 0.01;
	double tau = 0.01;
	ofstream out("task2.txt");
	// ���������� ��������� ��� ���������� ��������
	out << abs(b - a) << "    " << T << endl;
	solver(a2, f, g, phi, psi, fxx, a, b, T, h, tau, out);
	
}

void task3() {
	double a2 = 1.0;
	double T = 1.0;
	double a = 0;
	double b = 4*M_PI;
	auto f = [](double x) -> double {return 0.0;};
	auto g = [](double x) -> double {return 0.0; };
	auto phi = [](double x) -> double { return sin(x);};
	auto psi = [](double x) -> double { return 0.0; };
	auto fxx = [](double x) -> double { return 0.0; };

	// gamma=0.1
	/*double h = 0.1;
	double tau = 0.01;*/

	// gamma=0.5
	/*double h = 0.02;
	double tau = 0.01;*/

	//// gamma=0.75	
	double h = 0.15;
	double tau = 0.1125;

	//// gamma=1.0
	/*double h = 0.01;
	double tau = 0.01;*/
	ofstream out("task3.txt");
	// ���������� ��������� ��� ���������� ��������
	out << abs(b - a) << "    " << T << endl;
	solver(a2, f, g, phi, psi, fxx, a, b, T, h, tau, out);

}

void test_exact() {
	double a2 = 1.0;
	double T = 10.0;
	double a = 0.0;
	double b = M_PI;

	auto f = [](double x) -> double { return x; };
	auto g = [](double x) -> double { return x; };
	auto phi = [](double x) -> double { return 0.0; };
	auto psi = [](double x) -> double { return M_PI; };

	auto fxx = [](double x) -> double { return -M_PI * M_PI * sin(M_PI * x); };
	//double h = 0.01;
	//double tau = 0.05;
	double h = 0.02;
	double tau = 0.01;

	ofstream out("test1.txt");
	// ���������� ��������� ��� ���������� ��������
	out << abs(b - a) << "    " << T << endl;

	solver(a2, f, g, phi, psi, fxx, a, b, T, h, tau, out);

	out.close();
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "�� ��� �� �����!" << endl;

	/*test1();

	test2();*/

	task1();

	//task2();

	//task3();

	//test_exact();
	cout << "�� ��� � �� �����!" << endl;
}
