#include<iostream>
#include <functional>
#include <vector>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>


using namespace std;
using matrix = vector<vector<double>>;



// Вывод вектора
void print_vec(const matrix& res)
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

// метод прогонки
template<typename T>
vector<T> ThomasAlgorithm(
	const vector<T>& a, const vector<T>& b, const vector<T>& c,
	const vector<T>& d) {
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


// запись в файл матрицы
template<typename T>
ostream& operator<<(ostream& out, const vector<vector<T>>& v) {
	for (size_t i = 0; i < v.size(); i++) {
		for (size_t j = 0; j < v[0].size(); ++j)
			out << v[i][j] << "    ";
		out << endl;
	}
	return out;
}

//      x2
//      /\          f2
//    L2||----------------------
//      ||                     |
//      ||                     |
//  f3  ||                     |  f4
//      ||                     |
//      ||                     |
//      ||                     |
//      =============================>x1
//      0                      L1
//                  f1
// typeN - false - условия первого рода, true - второго, fN - функции, задающие условия
// eps - точность
// phi - правая часть уравнения: d^2(u)/dx^2 + d^2(u)/dy^2 = - phi(x, y)
template<bool type1 = false, bool type2 = false, bool type3 = false, bool type4 = false>
void solver(
	double L1, double L2, double eps, double h1, double h2, double tau,
	function<double(double, double)> phi,
	function<double(double)> f1, function<double(double)> f2,
	function<double(double)> f3, function<double(double)> f4,
	ostream& out)
{
	static_assert(!(type1 && type2 && type3 && type4), "Некорректно поставленная задача Неймана");

	size_t n1 = round(L1 / h1);
	size_t n2 = round(L2 / h2);

	// Инициализируем начальный слой и выделяем память под следующий и промежуточный
	matrix u0(n1 + 1, vector<double>(n2 + 1, 0));
	matrix u1(n1 + 1, vector<double>(n2 + 1, 0));
	matrix u2(n1 + 1, vector<double>(n2 + 1, 0));

	if (!type1) {
		// граничное условие 1 рода
		for (size_t i = 0; i <= n1; ++i)
			u0[i][0] = f1(i * h1);
	}
	else {
		// граничное условие 2 рода, просто записываем нули
		// TODO исправить (отдельно рассматривать узлы в углах)
		for (size_t i = 0; i <= n1; ++i)
			u0[i][0] = 0;
	}

	if (!type2) {
		// граничное условие 1 рода
		for (size_t i = 0; i <= n1; ++i)
			u0[i][n2] = f2(i * h1);
	}
	else {
		// граничное условие 2 рода, просто записываем нули
		// TODO исправить (отдельно рассматривать узлы в углах)
		for (size_t i = 0; i <= n1; ++i)
			u0[i][n2] = 0;
	}

	if (!type3) {
		// граничное условие 1 рода
		for (size_t j = 0; j <= n2; ++j)
			u0[0][j] = f3(j * h2);
	}
	else {
		// граничное условие 2 рода, просто записываем нули
		// TODO исправить (отдельно рассматривать узлы в углах)
		for (size_t j = 0; j <= n2; ++j)
			u0[0][j] = 0;
	}

	if (!type4) {
		// граничное условие 1 рода
		for (size_t j = 0; j <= n2; ++j)
			u0[n1][j] = f4(j * h2);
	}
	else {
		// граничное условие 2 рода, просто записываем нули
		// TODO исправить (отдельно рассматривать узлы в углах)
		for (size_t j = 0; j <= n2; ++j)
			u0[n1][j] = 0;
	}

	// Внутри области заполняем нулями (уже заполнено при объявлении)

	// Инициализируем матрицы СЛАУ (они трехдиагональные, нужно хранить шесть векторов)
	vector<double> A1(n1, 1.0 / h1 / h1), B1(n1 + 1, -2.0 * (1.0 / h1 / h1 + 1.0 / tau)), C1(n1, 1.0 / h1 / h1);
	vector<double> A2(n2, 1.0 / h2 / h2), B2(n2 + 1, -2.0 * (1.0 / h2 / h2 + 1.0 / tau)), C2(n2, 1.0 / h2 / h2);
	if (!type1) {
		B1[0] = 1.0;
		C1[0] = 0.0;
	}
	else {
		throw "ГУ 2 рода не реализованы";
	}
	if (!type2) {
		B1[n1] = 1.0;
		A1[n1 - 1] = 0.0;
	}
	else {
		throw "ГУ 2 рода не реализованы";
	}

	if (!type3) {
		B2[0] = 1.0;
		C2[0] = 0.0;
	}
	else {
		throw "ГУ 2 рода не реализованы";
	}
	if (!type4) {
		B2[n2] = 1.0;
		A2[n2 - 1] = 0.0;
	}
	else {
		throw "ГУ 2 рода не реализованы";
	}

	// Выделяем память под правые столбцы СЛАУ
	// TODO мог ошибиться с количеством
	vector<double> F1(n1 + 1), F2(n2 + 1);

	vector<double> u_vec;

	//////
	// кол-во итераций
	size_t k = 0;
	double oldNorm, newNorm = 0.0, nu = 0.0;
	//////

	while (k < 2 || newNorm > eps * (1 - nu)) {

		// TODO возможно, тут нужно отдельно рассматривать разные условия
		// и в зависимости от этого менять границы изменения j
		for (size_t j = 1; j <= n2 - 1; ++j) {
			// считаем F^k
			if (!type3) {
				F1[0] = f3(j * h2);
			}
			else {
				throw "ГУ 2 рода не реализованы";
			}
			for (size_t i = 1; i <= n1 - 1; ++i)
				F1[i] = -(1 / h2 / h2 * u0[i][j + 1] +
					2.0 * (1.0 / tau - 1.0 / h2 / h2) * u0[i][j] +
					1 / h2 / h2 * u0[i][j - 1] + phi(i * h1, j * h2));

			if (!type4) {
				F1[n1] = f4(j * h2);
			}
			else {
				throw "ГУ 2 рода не реализованы";
			}

			// прогонка по x
			u_vec = ThomasAlgorithm(A1, B1, C1, F1);

			// заполняем вычисленную строку
			for (size_t i = 0; i <= n1; ++i)
				u1[i][j] = u_vec[i];
		}

		// Непросчитанные строки j = 0, n2 копируем с прошлого слоя
		// TODO мб нужно действовать по-другому
		for (size_t i = 0; i <= n1; ++i) {
			u1[i][0] = u0[i][0];
			u1[i][n2] = u0[i][n2];
		}

		// TODO возможно, тут нужно отдельно рассматривать разные условия
		// и в зависимости от этого менять границы изменения j
		for (size_t i = 1; i <= n1 - 1; i++) {
			//считаем F^(k+1/2)
			if (!type1) {
				F2[0] = f1(i * h1);
			}
			else {
				throw "ГУ 2 рода не реализованы";
			}

			for (size_t j = 1; j <= n2 - 1; ++j)
				F2[j] = -(1 / h1 / h1 * u1[i + 1][j] +
					2.0 * (1.0 / tau - 1.0 / h1 / h1) * u1[i][j] +
					1 / h1 / h1 * u1[i - 1][j] + phi(i * h1, j * h2));

			if (!type2) {
				F2[n2] = f2(i * h1);
			}
			else {
				throw "ГУ 2 рода не реализованы";
			}

			// прогонка по y
			u_vec = ThomasAlgorithm(A2, B2, C2, F2);

			// заполняем вычисленный столбец
			for (size_t j = 0; j <= n2; ++j)
				u2[i][j] = u_vec[j];
		}

		// Непросчитанные столбцы i = 0, n1 копируем с прошлого слоя
		// TODO мб нужно действовать по-другому
		for (size_t j = 0; j <= n2; ++j) {
			u2[0][j] = u1[0][j];
			u2[n1][j] = u1[n1][j];
		}

		// Считаем норму разности (C-норму)
		double norm = 0.0;
		for (size_t i = 0; i <= n1; ++i)
			for (size_t j = 0; j <= n2; ++j)
				norm = max(norm, abs(u2[i][j] - u0[i][j]));

		// меняем критерий выхода
		k++;
		cout << "k = " << k << endl;
		if (k != 1) {
			oldNorm = newNorm;
			newNorm = norm;

			nu = newNorm - oldNorm;
		}
		else
			newNorm = norm;

		// Меняем значения
		u0 = u2;
	}

	out << u0;
}

// первый тест из методички (просто изменена площадка)
void test1() {
	double L1 = 1.0;
	double L2 = 2.0;
	double h1 = 0.001;
	double h2 = 0.001;
	double eps = 1e-3;
	double tau = 0.1;
	auto f = [](double x) { return 0.0; };
	auto phi = [](double x, double y) { return 0.0; };

	ofstream out("test1_1.txt");
	solver(L1, L2, eps, h1, h2, tau, phi, f, f, f, f, out);
	out.close();
}

// первый тест из методички с модификацией начального условия
void test1_2() {
	double L1 = 1.0;
	double L2 = 2.0;
	double h1 = 0.001;
	double h2 = 0.001;
	double eps = 1e-3;
	double tau = 0.1;
	auto f = [](double x) { return 1.0; };
	auto phi = [](double x, double y) { return 0.0; };

	ofstream out("test1_2.txt");
	solver(L1, L2, eps, h1, h2, tau, phi, f, f, f, f, out);
	out.close();
}

void my_test() {
	double L1 = 4.0;
	double L2 = 2.0;
	double h1 = 0.02;
	double h2 = 0.01;
	double eps = 1e-2;
	double tau = 0.1;
	auto f1 = [](double x) { return 6.0 - 3.0 / 4.0 * x; };
	auto f2 = [](double x) { return 3.0 - 3.0 / 4.0 * x; };
	auto f3 = [](double y) { return 6.0 - 3.0 / 2.0 * y; };
	auto f4 = [](double y) { return 3.0 - 3.0 / 2.0 * y; };
	auto phi = [](double x, double y) { return 0.0; };

	ofstream out("my_test.txt");
	solver(L1, L2, eps, h1, h2, tau, phi, f1, f2, f3, f4, out);
	out.close();
}

// решатель для подсчета порядка сходимости
//      x2
//      /\          f2
//    L2||----------------------
//      ||                     |
//      ||                     |
//  f3  ||                     |  f4
//      ||                     |
//      ||                     |
//      ||                     |
//      =============================>x1
//      0                      L1
//                  f1
// f1, f2, f3, f4 - функции, задающие условия 1 рода ()
// f0 - начальное распределение
// phi - правая часть уравнения: d^2(u)/dx^2 + d^2(u)/dy^2 = - phi(x, y)
// exact(t, x, y) - точное решение
// mod - на каждом m-м слое будет искатся ошибка
double solverOrder(
	double L1, double L2, double T, double h1, double h2, double tau,
	function<double(double, double)> phi,
	function<double(double)> f1, function<double(double)> f2,
	function<double(double)> f3, function<double(double)> f4,
	function<double(double, double)> f0,
	function<double(double, double, double)> exact,
	size_t mod = 1
	//ostream& out
)
{

	size_t n1 = round(L1 / h1);
	size_t n2 = round(L2 / h2);
	size_t m = round(T / tau);

	// Инициализируем начальный слой и выделяем память под следующий и промежуточный
	matrix u0(n1 + 1, vector<double>(n2 + 1, 0));
	matrix u1(n1 + 1, vector<double>(n2 + 1, 0));
	matrix u2(n1 + 1, vector<double>(n2 + 1, 0));

	for (size_t i = 0; i <= n1; ++i)
		for (size_t j = 0; j <= n2; ++j)
			u0[i][j] = f0(i * h1, j * h2);

	/*for (size_t i = 0; i <= n1; ++i)
		u0[i][0] = f1(i * h1);

	for (size_t i = 0; i <= n1; ++i)
		u0[i][n2] = f2(i * h1);

	for (size_t j = 0; j <= n2; ++j)
		u0[0][j] = f3(j * h2);

	for (size_t j = 0; j <= n2; ++j)
		u0[n1][j] = f4(j * h2);*/

		// Инициализируем матрицы СЛАУ (они трехдиагональные, нужно хранить шесть векторов)
	vector<double> A1(n1, 1.0 / h1 / h1), B1(n1 + 1, -2.0 * (1.0 / h1 / h1 + 1.0 / tau)), C1(n1, 1.0 / h1 / h1);
	vector<double> A2(n2, 1.0 / h2 / h2), B2(n2 + 1, -2.0 * (1.0 / h2 / h2 + 1.0 / tau)), C2(n2, 1.0 / h2 / h2);
	B1[0] = 1.0;
	C1[0] = 0.0;

	B1[n1] = 1.0;
	A1[n1 - 1] = 0.0;

	B2[0] = 1.0;
	C2[0] = 0.0;

	B2[n2] = 1.0;
	A2[n2 - 1] = 0.0;


	// Выделяем память под правые столбцы СЛАУ
	vector<double> F1(n1 + 1), F2(n2 + 1);

	vector<double> u_vec;

	double max_err = 0.0;

	for (size_t k = 1; k <= m; ++k) {
		for (size_t j = 1; j <= n2 - 1; ++j) {
			// считаем F^k
			F1[0] = f3(j * h2);
			for (size_t i = 1; i <= n1 - 1; ++i)
				F1[i] = -(1 / h2 / h2 * u0[i][j + 1] +
					2.0 * (1.0 / tau - 1.0 / h2 / h2) * u0[i][j] +
					1 / h2 / h2 * u0[i][j - 1] + phi(i * h1, j * h2));

			F1[n1] = f4(j * h2);

			// прогонка по x
			u_vec = ThomasAlgorithm(A1, B1, C1, F1);

			// заполняем вычисленную строку
			for (size_t i = 0; i <= n1; ++i)
				u1[i][j] = u_vec[i];
		}

		// Непросчитанные строки j = 0, n2 копируем с прошлого слоя
		// TODO мб нужно действовать по-другому
		for (size_t i = 0; i <= n1; ++i) {
			u1[i][0] = u0[i][0];
			u1[i][n2] = u0[i][n2];
		}

		for (size_t i = 1; i <= n1 - 1; i++) {
			//считаем F^(k+1/2)
			F2[0] = f1(i * h1);

			for (size_t j = 1; j <= n2 - 1; ++j)
				F2[j] = -(1 / h1 / h1 * u1[i + 1][j] +
					2.0 * (1.0 / tau - 1.0 / h1 / h1) * u1[i][j] +
					1 / h1 / h1 * u1[i - 1][j] + phi(i * h1, j * h2));

			F2[n2] = f2(i * h1);

			// прогонка по y
			u_vec = ThomasAlgorithm(A2, B2, C2, F2);

			// заполняем вычисленный столбец
			for (size_t j = 0; j <= n2; ++j)
				u2[i][j] = u_vec[j];
		}

		// Непросчитанные столбцы i = 0, n1 копируем с прошлого слоя
		for (size_t j = 0; j <= n2; ++j) {
			u2[0][j] = u1[0][j];
			u2[n1][j] = u1[n1][j];
		}

		// Считаем норму разности (C-норму)
		double norm = 0.0;
		for (size_t i = 0; i <= n1; ++i)
			for (size_t j = 0; j <= n2; ++j)
				norm = max(norm, abs(u2[i][j] - u0[i][j]));

		// Меняем значения
		u0 = u2;

		// считаем норму ошибки
		double err = 0.0;
		if (k % mod == 0)
			for (size_t i = 1; i < n1; ++i)
				for (size_t j = 1; j < n2; ++j)
					err = max(err, abs(u0[i][j] - exact(tau * k, i * h1, j * h2)));

		//max_err = err;
		max_err = max(max_err, err);
	}

	return max_err;
	//out << u0;
}

// Тест для оценки порядка метода
double order_test(size_t k = 1) {
	double L1 = 2.0;
	double L2 = 1.0;
	double h1 = 0.02 / pow(2, k - 1);
	double h2 = 0.01 / pow(2, k - 1);
	double T = 1.0;
	double tau = 0.1 / pow(2, k - 1);
	auto f = [](double x) { return 0.0; };
	auto f0 = [](double x, double y) { return sin(M_PI * x / 2.0) * sin(M_PI * y); };
	auto phi = [](double x, double y) { return 0.0; };
	auto exact = [](double t, double x, double y) {
		return exp(-5.0 * M_PI * M_PI / 4.0 * t) * sin(M_PI * x / 2.0) * sin(M_PI * y);
	};

	//ofstream out("test.txt");
	double max_err = solverOrder(L1, L2, T, h1, h2, tau, phi, f, f, f, f, f0, exact, k);

	return max_err;
	//out.close();
}

// Тест для поиска времени выхода на стационар
// TODO
void stationary_test(double tau) {
	cout << tau << endl;
	double L1 = 4.0;
	double L2 = 2.0;
	double h1 = 0.02;
	double h2 = 0.01;
	//double eps = 1e-2;
	double eps = 1e-2;
	auto f1 = [](double x) { return 6.0 - 3.0 / 4.0 * x; };
	auto f2 = [](double x) { return 3.0 - 3.0 / 4.0 * x; };
	auto f3 = [](double y) { return 6.0 - 3.0 / 2.0 * y; };
	auto f4 = [](double y) { return 3.0 - 3.0 / 2.0 * y; };
	auto phi = [](double x, double y) { return 0.0; };

	ofstream out("test.txt");
	solver(L1, L2, eps, h1, h2, tau, phi, f1, f2, f3, f4, out);
	out.close();
}


int main()
{
	//test();
	test1();
	/*test1_2();
	my_test();*/

	// таблицы
	//stationary_test(0.02);
	//stationary_test(0.02 / 2);
	//stationary_test(0.02 / 2 / 2);
	//stationary_test(0.02 / 2 / 2 / 2);
	//stationary_test(0.02 / 2 / 2 / 2 / 2);

	// порядок схемы
	// k = 1
	/*double absErr1 = order_test(1), absErr2;
	cout << 1 << "\t" << absErr1 << endl;
	for (size_t k = 2; k <= 6; ++k) {
		absErr2 = order_test(k);
		double delta = absErr1 / absErr2;
		double logDelta = log2(delta);
		cout << k << "\t" << absErr2 << "\t" << delta << "\t" << logDelta << endl;
		absErr1 = absErr2;
	}*/

}