#include <functional>
#include <fstream>
#include <iomanip>
#include <string>
#include<iostream>
#include<vector>

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;

// Обращение матрицы 2x2 (просто ввел формулы из Wolfram Mathematica)
template<typename T>
matrix<T> inv(const matrix<T>& m) {
	T det = m[0][0] * m[1][1] - m[1][0] * m[0][1];
	return matrix<T> {
		{  m[1][1] / det, -m[0][1] / det },
		{ -m[1][0] / det,  m[0][0] / det }
	};
}

// умножение матрицы и столбца
template<typename T>
vector<T> mult(const matrix<T>& A1, const vector<T>& A2) {
	size_t n = A2.size();
	vector<T> res;
	res.reserve(n); // оптимизируем память
	for (size_t i = 0; i < n; ++i) {
		T c = 0;
		for (size_t j = 0; j < n; ++j) {
			c += A1[i][j] * A2[j];
		}
		res.push_back(c);

	}
	return res;
}

double norm(double x) {
	return abs(x);
}

double norm(const vector<double>& vec) {
	double res = 0.0;
	for (double el : vec)
		res += abs(el);
	return res;
}

////////////////////////////////////////////


//Вывод вектора
void print_vec(const vector<vector<double>>& res, double h)
{

	for (size_t i = 0; i < res.size(); i++) {
		cout << i * h << "\t";
		for (size_t j = 0; j < res[i].size(); j++)
			cout << fixed << res[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}

void print_vec(const vector<vector<double>>& res)
{

	for (size_t i = 0; i < res.size(); i++) {
		for (size_t j = 0; j < res[i].size(); j++)
			cout << res[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}


void print_vec(const pair<vector<double>, vector<vector<double>>>& res)
{

	for (size_t i = 0; i < res.second.size(); i++) {
		cout << res.first[i] << "\t";
		for (size_t j = 0; j < res.second[i].size(); j++)
			cout << fixed << res.second[i][j] << "\t";
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


void inFile(const vector<vector<double>>& res, string str)
{
	ofstream out(str);
	for (size_t i = 0; i < res.size(); i++) {
		for (size_t j = 0; j < res[i].size(); j++)
			if (j == res[i].size() - 1)
				out << res[i][j];
			else
				out << res[i][j] << "\t";
		out << endl;
	}
	out.close();
}



//Вариант 1 Уравнение Ван-дер-Поля
// x1(0)=0.1, x2(0)=0.1 t=0...200
vector<double> variant1(double t, vector<double> x) {
	double delta = 0.3;
	double omega = 1.0;
	double alpha = 1.0;
	return vector<double>{x[1], 2 * delta * x[1] * (1.0 - alpha * x[0] * x[0]) - omega * omega * x[0] };
}

matrix<double> variant1Diff(double t, vector<double> x) {
	double delta = 0.3;
	double omega = 1.0;
	double alpha = 1.0;

	return matrix<double>{ {0.0, 1.0},
		{ -4.0 * alpha * delta * x[0] * x[1] - omega * omega, -2.0 * delta * alpha * x[0] * x[0] }};
}


//Вариант 20 Модель Лотки-Вольтерры
// x1(0)=1.0, x2(0)=4.0 t=0...150
vector<double> variant20(double t, vector<double> x) {
	double r1 = 0.4;
	double r2 = 0.1;
	double b11 = 0.05;
	double b12 = 0.1;
	double b21 = 0.08;
	double b22 = 0.003;
	return vector<double>{ r1* x[0] - b11 * x[0] * x[0] - b12 * x[0] * x[1], -r2 * x[1] - b22 * x[1] * x[1] + b21 * x[0] * x[1]  };
}



double mx1(double x = 0, double y = 0) {
	return y;
}

double mx2(double x = 0, double y = 0) {
	return 20 / 0.3 * x;
}

// возвращает функцию, описывающую математический маятник
function<vector<double>(double t, vector<double>)> mx(double k = 20.0, double m = 0.3) {
	double omega2 = k / m;
	return [omega2](double t, vector<double> x) -> vector<double> {
		return vector<double>{x[1], -omega2 * x[0]};
	};
}

// якобиан системы математического маятника
function<vector<vector<double>>(double t, vector<double>)> mxDiff(double k = 20.0, double m = 0.3) {
	double omega2 = k / m;
	return [omega2](double t, vector<double> x) -> vector<vector<double>> {
		return vector<vector<double>>{ {0.0, 1.0}, { -omega2, 0.0 }};
	};
}


double test1(double x, double y) {
	return 4 * y - 5 * x;
}

double test2(double x, double y) {
	return 4 * x;
}

vector<double> f1(double t, vector<double> x) {
	return vector<double>{ 2 * x[0] + x[1] * x[1] - 1, 6 * x[0] - x[1] * x[1] + 1 };
}


vector<double> f2(double t, vector<double> x) {
	return vector<double>{1 - x[0] * x[0] - x[1] * x[1], 2 * x[0] };
}


vector<double> f3(double t, vector<double> x) {
	return vector<double>{10.0 * (x[1] - x[0]), x[0] * (28.0 - x[2]) - x[1], x[0] * x[1] - 8 / 3 * x[2]};
}



// Явный метод Эйлера для маятника
vector<vector<double>> Euler_explicit_mx(double T = 10.0, double h = 0.01, double k = 20, double m = 0.3) {
	vector<vector<double>> res;
	vector<double> resn(2, 0), res0{ 1,1 };
	res.push_back(res0);
	for (size_t t = 0; t < int(T / h); t++) {
		resn[0] = res0[0] + res0[1] * h;
		resn[1] = res0[1] - h * (k / m) * res0[0];
		res0 = resn;
		res.push_back(resn);
	}

	return res;
}

//Неявный метод Эйлера для маятника
vector<vector<double>> Euler_implicit_mx(double T = 10.0, double h = 0.01, double k = 20, double m = 0.3) {
	vector<vector<double>> res;
	vector<double> resn(2, 0), res0{ 1,1 };
	double w2 = (k / m);
	res.push_back(res0);
	for (size_t t = 0; t < int(T / h); t++) {
		resn[0] = res0[0] * 1 / (1 + h * h * w2) + res0[1] * h / (1 + h * h * w2);
		resn[1] = res0[0] * (-1) * (h * w2) / (1 + h * h * w2) + res0[1] * 1 / (1 + h * h * w2);
		res0 = resn;
		res.push_back(resn);
	}

	return res;
}

// базовые операции над векторами

// сложение векторов
template<typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b) {
	vector<T> res;
	res.reserve(a.size());
	for (size_t i = 0; i < a.size(); i++)
		res.push_back(a[i] + b[i]);
	return res;
}

// вычитание векторов
template<typename T>
vector<T> operator-(const vector<T>& a, const vector<T>& b) {
	vector<T> res;
	res.reserve(a.size());
	for (size_t i = 0; i < a.size(); i++)
		res.push_back(a[i] - b[i]);
	return res;
}

// умножение вектора на число
template<typename T>
vector<T> operator*(const vector<T>& a, T b) {
	vector<T> res;
	res.reserve(a.size());
	for (size_t i = 0; i < a.size(); i++)
		res.push_back(a[i] * b);
	return res;
}

// умножение числа на вектор
template<typename T>
vector<T> operator*(T b, const vector<T>& a) {
	vector<T> res;
	res.reserve(a.size());
	for (size_t i = 0; i < a.size(); i++)
		res.push_back(a[i] * b);
	return res;
}



// Метод Рунге - Кутты второго порядка
template<typename T>
pair<vector<double>, vector<T>> RK2(function<T(double, T)> func,
	pair<double, double> interval, T initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<T> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);


	for (size_t i = 0; i < k; i++) {
		double t = tt.back();
		T x = xx.back();

		T k1 = func(t, x);
		T k2 = func(t + h / 2, x + (h / 2) * k1);

		T newX = x + h * k2;

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
	}

	return make_pair(tt, xx);

}

// Метод Рунге - Кутты четвертого порядка
// func - правая часть дифференциального уравнения, первый аргумент - время t
// interval - интервал интегрирования
// initVal - начальное значение
// h - шаг
// T может быть числом или вектором
template<typename T>
pair<vector<double>, vector<T>> RK4(function<T(double, T)> func,
	pair<double, double> interval, T initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<T> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);


	for (size_t i = 0; i < k; i++) {
		double t = tt.back();
		T x = xx.back();

		T k1 = func(t, x);
		T k2 = func(t + h / 2, x + (h / 2) * k1);
		T k3 = func(t + h / 2, x + (h / 2) * k2);
		T k4 = func(t + h, x + h * k3);

		T newX = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
	}

	return make_pair(tt, xx);

}


// Явный метод Эйлера
// func - правая часть дифференциального уравнения, первый аргумент - время t
// interval - интервал интегрирования
// initVal - начальное значение
// h - шаг
// T может быть числом или вектором
template<typename T>
pair<vector<double>, vector<T>>  Euler_explicit(function<T(double, T)> func,
	pair<double, double> interval, T initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<T> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);


	for (size_t i = 0; i < k; i++) {
		double t = tt.back();
		T x = xx.back();

		T newX = x + h * func(t, x);

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
	}

	return make_pair(tt, xx);
}

//// Метод Ньютона (без модификации)
//// diff - производная функции fun
//template<typename T, typename F, typename Diff>
//T Newton(F fun, Diff diff, T a, T b, T x0, T eps, size_t& iterCount, size_t p = 1) {
//	T xk = x0;
//	T coef;
//	vector<T> X;
//	X.push_back(xk);
//	do {
//		T newXk = xk - p * fun(xk) / diff(xk);
//		coef = abs(newXk - xk);
//		xk = newXk;
//		X.push_back(xk);
//		++iterCount;
//		//cout << newXk <<"\t" << abs(newXk - t7exakt) << endl;
//	} while (coef > eps);
//	degreeX(X, t7exakt);
//
//
//
//


//// Метод Ньютона (без модификации)
//// diff - производная функции fun
template<typename T, typename F, typename Diff>
//T Newton(F fun, Diff diff, T x0, T eps, size_t p = 1) {
T Newton(F fun, Diff diff, T x0, T eps) {
	T xk = x0;
	T coef;

	do {
		T newXk = xk - fun(xk) / diff(xk);
		coef = abs(newXk - xk);
		xk = newXk;
	} while (coef > eps);

	return xk;
}


vector<double> NewtonSystem(function<vector<double>(vector<double>)> f,
	function<matrix<double>(vector<double>)> jacobi,
	double eps,
	const vector<double>& x0)
{
	vector<double> xk = x0;
	double coef;

	size_t iterCount = 0;

	do {
		// F'^-1
		matrix<double> jacobiInv = inv(jacobi(xk));
		vector<double> F = f(xk);

		// F'^-1 * F
		vector<double> jacobiInvMultF = mult(jacobiInv, F);

		// x^k+1 = x^k - F'^-1 * F
		// x^k+1 - x^k = - F'^-1 * F
		// |x^k+1 - x^k| = |- F'^-1 * F| = |F'^-1 * F|
		// => coef = norm(jacobiInvMultF)
		vector<double> newXk;
		// Типо задел под n-мерный случай
		newXk.reserve(x0.size());
		for (size_t i = 0; i < x0.size(); ++i) {
			newXk.push_back(xk[i] - jacobiInvMultF[i]);
		}

		coef = norm(jacobiInvMultF);

		xk = newXk;

	} while (coef > eps);

	return xk;
}

// Неявный метод Эйлера для одного уравнений
pair<vector<double>, vector<double>> Euler1D_implicit(function<double(double, double)> func,
	function<double(double, double)> diff,
	pair<double, double> interval, double initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<double> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);


	for (size_t i = 0; i < k; i++) {
		double t = tt.back();
		double x = xx.back();

		// составляем решаемое уравнение g(x) = 0 и необходимая для метода
		// Ньютона производная
		function<double(double)> g = [func, t, h, x](double xk) -> double {
			return xk - x - h * func(t + h, xk);
		};
		// Якобиан g'(x) = 1 - h*f'(t+h, x), f' -  производная правой части по x
		function<double(double)> gDiff = [diff, t, h](double xk) -> double {
			return 1.0 - h * diff(t + h, xk);
		};
		double newX = Newton(g, gDiff, h, x);

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
	}

	return make_pair(tt, xx);
}


// Симметричная схема для одного уравнений
pair<vector<double>, vector<double>> Symm1D(function<double(double, double)> func,
	function<double(double, double)> diff,
	pair<double, double> interval, double initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<double> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);


	for (size_t i = 0; i < k; i++) {
		double t = tt.back();
		double x = xx.back();

		// коэффицинт для упрощения вычислений
		double c = x + (h / 2) * func(t, x);

		// составляем решаемое уравнение g(x) = 0 и необходимая для метода
		// Ньютона производная
		function<double(double)> g = [func, t, h, x, c](double xk) -> double {
			return xk - c - (h / 2) * func(t + h, xk);
		};
		// Якобиан g'(x) = 1 - (h/2)*f'(t+h, x), f' -  производная правой части по x
		function<double(double)> gDiff = [diff, t, h](double xk) -> double {
			return 1.0 - (h / 2) * diff(t + h, xk);
		};
		double newX = Newton(g, gDiff, h, x);

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
	}

	return make_pair(tt, xx);
}

//Неявный метод Эйлера для систем
template<typename T = vector<double>>
pair<vector<double>, vector<T>> EulerSystem_implicit(function<T(double, T)> func,
	function<matrix<double>(double, T)> diff,
	pair<double, double> interval, T initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<T> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);


	for (size_t i = 0; i < k; i++) {
		double t = tt.back();
		T x = xx.back();

		// составляем решаемое уравнение g(x) = 0 и необходимый для метода
		// Ньютона якобиан
		function<T(T)> g = [func, t, h, x](T xk) -> T {
			return xk - x - h * func(t + h, xk);
		};
		// Якобиан g'(x) = E - h*f'(t+h, x), f' -  якобиан правой части по x
		function<matrix<double>(T)> gDiff = [diff, t, h](T xk) -> matrix<double> {
			matrix<double> d = diff(t + h, xk);
			size_t n = d.size();
			for (size_t i = 0; i < n - 1; ++i) {
				d[i][i] = 1 - h * d[i][i];
				for (size_t j = i + 1; j < n; ++j) {
					d[i][j] *= -h;
					d[j][i] *= -h;
				}
			}
			d[n - 1][n - 1] = 1 - h * d[n - 1][n - 1];
			return d;
		};
		T newX = NewtonSystem(g, gDiff, h, x);

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
	}

	return make_pair(tt, xx);
}

// симметричная схема для систем
template<typename T = vector<double>>
pair<vector<double>, vector<T>> SymmSystem(function<T(double, T)> func,
	function<matrix<double>(double, T)> diff,
	pair<double, double> interval, T initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<T> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);


	for (size_t i = 0; i < k; i++) {
		double t = tt.back();
		T x = xx.back();

		//T newX = x + h * func(t, x);

		// коэффицинт для упрощения вычислений
		T c = x + (h / 2) * func(t, x);

		// составляем решаемое уравнение g(x) = 0 и необходимый для метода
		// Ньютона якобиан
		function<T(T)> g = [func, t, h, c](T xk) -> T {
			return xk - c - (h / 2) * func(t + h, xk);
		};
		// Якобиан g'(x) = E - (h/2)*f'(t+h, x), f' -  якобиан правой части по x
		function<matrix<double>(T)> gDiff = [diff, t, h](T xk) -> matrix<double> {
			matrix<double> d = diff(t + h, xk);
			size_t n = d.size();
			for (size_t i = 0; i < n - 1; ++i) {
				d[i][i] = 1 - (h / 2) * d[i][i];
				for (size_t j = i + 1; j < n; ++j) {
					d[i][j] *= -h / 2;
					d[j][i] *= -h / 2;
				}
			}
			d[n - 1][n - 1] = 1 - (h / 2) * d[n - 1][n - 1];
			return d;
		};
		T newX = NewtonSystem(g, gDiff, h, x);

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
	}

	return make_pair(tt, xx);
}

// TODO не уверен насчёт порядка (точность вроде бы хуже, чем у Рунге - Кутты)
// Метод Адамса 4 порядка
// func - правая часть дифференциального уравнения, первый аргумент - время t
// interval - интервал интегрирования
// initVal - начальное значение
// h - шаг
// T может быть числом или вектором
template<typename T>
pair<vector<double>, vector<T>> Adams(function<T(double, T)> func,
	pair<double, double> interval, T initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<T> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);

	vector<T> fn{ func(interval.first, initVal) };
	fn.reserve(4);

	// заполняем fn методом Рунге - Кутты 4 порядка (ещё 3 значения)
	for (size_t i = 0; i < 3; i++) {
		double t = tt.back();
		T x = xx.back();

		T k1 = fn.back();
		T k2 = func(t + h / 2, x + (h / 2) * k1);
		T k3 = func(t + h / 2, x + (h / 2) * k2);
		T k4 = func(t + h, x + h * k3);

		T newX = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

		// метод Рунге-Кутты второго порядка (модифицированный метод Эйлера)
		// для того, чтобы показать, что порядок метода зависит от внутреннего
		// метода, вычисляющего первые значения
		// Но оно никак не влияет
		/*T k1 = fn.back();
		T k2 = func(t + h / 2, x + (h / 2) * k1);
		T newX = x + h * k2;*/

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
		fn.push_back(func(t + h, newX));
	}

	// сам метод Адамса
	for (size_t i = 3; i < k; i++) {
		T x = xx.back();
		double newT = interval.first + (i + 1) * h;

		// расчёт нового значения
		T newX = x + (h / 24) * (55.0 * fn[3] + (-59.0) * fn[2] + 37.0 * fn[1] + (-9.0) * fn[0]);

		tt.push_back(newT);
		xx.push_back(newX);

		// меняем коэффициенты
		fn.erase(fn.begin());
		fn.push_back(func(newT, newX));
	}

	return make_pair(tt, xx);
}

// TODO не уверен насчёт порядка (точность вроде бы хуже, чем у Рунге - Кутты)
// Метод "прогноз-коррекция" 4 порядка
// func - правая часть дифференциального уравнения, первый аргумент - время t
// interval - интервал интегрирования
// initVal - начальное значение
// h - шаг
// T может быть числом или вектором
template<typename T>
pair<vector<double>, vector<T>> PEC(function<T(double, T)> func,
	pair<double, double> interval, T initVal, double h) {

	size_t k = (interval.second - interval.first) / h;
	vector<double> tt;
	tt.reserve(k + 1);
	vector<T> xx;
	xx.reserve(k + 1);

	tt.push_back(interval.first);
	xx.push_back(initVal);

	vector<T> fn{ func(interval.first, initVal) };
	fn.reserve(4);

	// заполняем fn методом Рунге - Кутты 4 порядка (ещё 3 значения)
	for (size_t i = 0; i < 3; i++) {
		double t = tt.back();
		T x = xx.back();

		T k1 = fn.back();
		T k2 = func(t + h / 2, x + (h / 2) * k1);
		T k3 = func(t + h / 2, x + (h / 2) * k2);
		T k4 = func(t + h, x + h * k3);

		T newX = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

		// метод Рунге-Кутты второго порядка (модифицированный метод Эйлера)
		// для того, чтобы показать, что порядок метода зависит от внутреннего
		// метода, вычисляющего первые значения
		// Но оно никак не влияет
		/*T k1 = fn.back();
		T k2 = func(t + h / 2, x + (h / 2) * k1);
		T newX = x + h * k2;*/

		tt.push_back(interval.first + (i + 1) * h);
		xx.push_back(newX);
		fn.push_back(func(t + h, newX));
	}

	// сам метод
	for (size_t i = 3; i < k; i++) {
		T x = xx.back();
		double newT = interval.first + (i + 1) * h;

		// расчёт нового значения
		T newX0 = x + (h / 24) * (55.0 * fn[3] + (-59.0) * fn[2] + 37.0 * fn[1] + (-9.0) * fn[0]);
		T newX1 = x + (h / 24) * (9.0 * func(newT, newX0) + 19.0 * fn[3] + (-5.0) * fn[2] + fn[1]);

		tt.push_back(newT);
		xx.push_back(newX1);

		// меняем коэффициенты
		fn.erase(fn.begin());
		fn.push_back(func(newT, newX1));
	}

	return make_pair(tt, xx);
}

// функция-пример с лабы
double exactTest(double t) {
	return exp(sin(t * t));
}

double fTest(double t, double x) {
	return 2.0 * t * cos(t * t) * x;
}

// наш пример
double labTest(double t) {
	return t * t * log(t) * sin(t / 5.0);
}

double labTestDiff(double t, double x) {
	return 0.2 * x / tan(t / 5.0) + 2.0 * x / t + t * sin(t / 5.0);
}
double labTestDiffDiff(double t, double x) {
	return 2.0 / t + 0.2 / tan(t / 5.0);
}

enum method {
	M_EULER_EXPL, M_EULER_IMPL, M_SYMM, M_RK2, M_RK4, M_ADAMS, M_PEC
};

// Вычисление порядка (с помощью эталонного решения)
// exact - точное решение
// func - производная exact (правая часть дифференциального уравнения)
// meth - исследуемый метод
// point - в какой точке исследуется порядок (SIZE_MAX - в последней)
// initH - начальный шаг
// count - количество делений шага (вычисления порядка)
void order(function<double(double)> exact, function<double(double, double)> func, function<double(double, double)> diff,
	pair<double, double> interval, double initVal, method meth, size_t point = SIZE_MAX, double initH = 0.1, size_t count = 6) {

	// ширина столбца при выводе
	size_t w = 14;


	double h = initH;
	double prev, cur;
	pair<vector<double>, vector<double>> res;

	//auto res = Euler_explicit<double>(fTest, interval, initVal, h);
	//auto res = RK2<double>(fTest, interval, initVal, h);
	//res = RK4<double>(fTest, interval, initVal, h);
	//auto res = Adams<double>(fTest, interval, initVal, h);
	//auto res = PEC<double>(fTest, interval, initVal, h);
	switch (meth) {
	case M_EULER_EXPL:
		res = Euler_explicit<double>(func, interval, initVal, h);
		break;
	case M_EULER_IMPL:
		res = Euler1D_implicit(func, diff, interval, initVal, h);
		break;
	case M_SYMM:
		res = Symm1D(func, diff, interval, initVal, h);
		break;
	case M_RK2:
		res = RK2<double>(func, interval, initVal, h);
		break;
	case M_RK4:
		res = RK4<double>(func, interval, initVal, h);
		break;
	case M_ADAMS:
		res = Adams<double>(func, interval, initVal, h);
		break;
	case M_PEC:
		res = PEC<double>(func, interval, initVal, h);
		break;
	default:
		throw "Метод не поддерживается";
		break;
	}

	//prev = abs(exact(res.first.back()) - res.second.back());
	if (point == SIZE_MAX)
		prev = abs(exact(res.first.back()) - res.second.back());
	else
		prev = abs(exact(res.first[point]) - res.second[point]);
	//cout << res.first[10] << " " << res.second[10] << endl;
	//cout << res.first[20] << " " << res.second[20] << endl;
	cout << left << setw(w) << "tau" << setw(w) << "AbsErr(tau)" << setw(w) << "Delta" << setw(w) << "log(Delta)" << endl;
	cout << setw(w) << h << setw(w) << prev << setw(w) << "-" << setw(w) << "-" << endl;


	h /= 2.0;

	for (size_t i = 0; i < count; i++) {
		//res = Euler_explicit<double>(fTest, interval, initVal, h);
		//res = RK2<double>(fTest, interval, initVal, h);
		//res = RK4<double>(fTest, interval, initVal, h);
		//res = Adams<double>(fTest, interval, initVal, h);
		//res = PEC<double>(fTest, interval, initVal, h);
		switch (meth) {
		case M_EULER_EXPL:
			res = Euler_explicit<double>(func, interval, initVal, h);
			break;
		case M_EULER_IMPL:
			res = Euler1D_implicit(func, diff, interval, initVal, h);
			break;
		case M_SYMM:
			res = Symm1D(func, diff, interval, initVal, h);
			break;
		case M_RK2:
			res = RK2<double>(func, interval, initVal, h);
			break;
		case M_RK4:
			res = RK4<double>(func, interval, initVal, h);
			break;
		case M_ADAMS:
			res = Adams<double>(func, interval, initVal, h);
			break;
		case M_PEC:
			res = PEC<double>(func, interval, initVal, h);
			break;
		default:
			throw "Метод не поддерживается";
			break;
		}

		//cur = abs(exact(res.first.back()) - res.second.back());
		if (point == SIZE_MAX)
			cur = abs(exact(res.first.back()) - res.second.back());
		else
			//cur = abs(exact(res.first[point * (i + 1)]) - res.second[point * (i + 1)]);
			cur = abs(exact(res.first[point * round(pow(2, i + 1))]) - res.second[point * round(pow(2, i + 1))]);
		cout << setw(w) << h << setw(w) << cur << setw(w) << (prev / cur) << setw(w) << log2(prev / cur) << endl;


		prev = cur;
		h /= 2.0;
	}

}

// TODO сделать
// Вычисление порядка (с помощью правила Эйткена)
// exact - точное решение
// func - производная exact (правая часть дифференциального уравнения)
// meth - исследуемый метод
// point - в какой точке исследуется порядок (SIZE_MAX - в последней)
// initH - начальный шаг
// count - количество делений шага (вычисления порядка)
void eitken(function<double(double, double)> func, function<double(double, double)> diff,
	pair<double, double> interval, double initVal, method meth, size_t point = SIZE_MAX, double initH = 0.1, size_t count = 6) {

	// ширина столбца при выводе
	size_t w = 14;


	double h = initH;
	double y1, y2, y3;
	pair<vector<double>, vector<double>> res;

	switch (meth) {
	case M_EULER_EXPL:
		res = Euler_explicit<double>(func, interval, initVal, h);
		break;
	case M_EULER_IMPL:
		res = Euler1D_implicit(func, diff, interval, initVal, h);
		break;
	case M_SYMM:
		res = Symm1D(func, diff, interval, initVal, h);
		break;
	case M_RK2:
		res = RK2<double>(func, interval, initVal, h);
		break;
	case M_RK4:
		res = RK4<double>(func, interval, initVal, h);
		break;
	case M_ADAMS:
		res = Adams<double>(func, interval, initVal, h);
		break;
	case M_PEC:
		res = PEC<double>(func, interval, initVal, h);
		break;
	default:
		throw "Метод не поддерживается";
		break;
	}

	//prev = abs(exact(res.first.back()) - res.second.back());
	if (point == SIZE_MAX)
		y1 = res.second.back();
	else
		y1 = res.second[point];

	cout << left << setw(w) << "tau" << setw(w) << "ErrEstim(tau)" << setw(w) << "pEstim" << endl;
	cout << setw(w) << h << setw(w) << "-" << setw(w) << "-" << endl;

	h /= 2.0;

	switch (meth) {
	case M_EULER_EXPL:
		res = Euler_explicit<double>(func, interval, initVal, h);
		break;
	case M_EULER_IMPL:
		res = Euler1D_implicit(func, diff, interval, initVal, h);
		break;
	case M_SYMM:
		res = Symm1D(func, diff, interval, initVal, h);
		break;
	case M_RK2:
		res = RK2<double>(func, interval, initVal, h);
		break;
	case M_RK4:
		res = RK4<double>(func, interval, initVal, h);
		break;
	case M_ADAMS:
		res = Adams<double>(func, interval, initVal, h);
		break;
	case M_PEC:
		res = PEC<double>(func, interval, initVal, h);
		break;
	default:
		throw "Метод не поддерживается";
		break;
	}

	if (point == SIZE_MAX)
		y2 = res.second.back();
	else
		y2 = res.second[point * 2];

	cout << setw(w) << h << setw(w) << abs(y2 - y1) << setw(w) << "-" << endl;

	h /= 2.0;

	for (size_t i = 1; i < count; i++) {
		switch (meth) {
		case M_EULER_EXPL:
			res = Euler_explicit<double>(func, interval, initVal, h);
			break;
		case M_EULER_IMPL:
			res = Euler1D_implicit(func, diff, interval, initVal, h);
			break;
		case M_SYMM:
			res = Symm1D(func, diff, interval, initVal, h);
			break;
		case M_RK2:
			res = RK2<double>(func, interval, initVal, h);
			break;
		case M_RK4:
			res = RK4<double>(func, interval, initVal, h);
			break;
		case M_ADAMS:
			res = Adams<double>(func, interval, initVal, h);
			break;
		case M_PEC:
			res = PEC<double>(func, interval, initVal, h);
			break;
		default:
			throw "Метод не поддерживается";
			break;
		}

		if (point == SIZE_MAX)
			y3 = res.second.back();
		else
			y3 = res.second[point * round(pow(2, i + 1))];
		cout << setw(w) << h << setw(w) << abs(y3 - y2) << setw(w) << log2(abs((y1 - y2) / (y2 - y3))) << endl;

		y1 = y2;
		y2 = y3;
		h /= 2.0;
	}

}

void pend_err(vector<double> x0 = { 1.0, 1.0 }, double k = 20, double m = 0.3) {
	double omega = sqrt(k / m);

	function<vector<double>(double t)> exakt = [x0, omega](double t) -> vector<double> {
		return vector<double>{x0[0] * cos(t* omega) + x0[1] / omega * sin(t * omega),
			x0[1] * cos(t* omega) - x0[0] * omega * sin(t * omega)
		};
	};

	auto pendulum = mx(k, m);
	auto pendDiff = mxDiff(k, m);

	auto interval = make_pair(0.0, 10.0);

	double tau1 = 0;
	double tau2 = 1.0;
	//double tau2 = 3.71039e-6;

	double tau;

	double eps = 1e-2;
	//double eps = 1e-4;
	//double eps = 1e-7;

	double err;

	while (true) {
		double newErr = 0.0;

		tau = (tau1 + tau2) / 2;

		auto res = Euler_explicit<vector<double>>(pendulum, interval, x0, tau);
		//auto res = EulerSystem_implicit<vector<double>>(pendulum, pendDiff, interval, x0, tau);
		//auto res = SymmSystem<vector<double>>(pendulum, pendDiff, interval, x0, tau);
		//auto res = RK2<vector<double>>(pendulum, interval, x0, tau);
		//auto res = RK4<vector<double>>(pendulum, interval, x0, tau);
		//auto res = Adams<vector<double>>(pendulum, interval, x0, tau);
		//auto res = PEC<vector<double>>(pendulum, interval, x0, tau);

		for (size_t i = 0; i < res.second.size(); ++i) {
			vector<double> val = exakt(res.first[i]);
			vector<double> approx = res.second[i];

			double norm = sqrt(pow(val[0] - approx[0], 2.0) + pow(val[1] - approx[1], 2.0));

			/*if (norm > eps) {
				err = DBL_MAX;
				break;
			}*/

			newErr = max(newErr, norm);

			//err = min(err, norm);
		}

		err = newErr;

		if (abs(err - eps) < eps * 1e-3)
			break;
		else if (err > eps)
			tau2 = tau;
		else
			tau1 = tau;

		cout << "[" << tau1 << ", " << tau2 << "]: " << err << endl;
	}

	cout << "Нужный шаг: " << tau << "; ошибка: " << err << endl;
}


// Шаг метода Рунге - Кутты (используется в РК с правилом Рунге)
template<typename T>
inline T RK4step(function<T(double, T)> func, double t, T x, double h) {
	T k1 = func(t, x);
	T k2 = func(t + h / 2, x + (h / 2) * k1);
	T k3 = func(t + h / 2, x + (h / 2) * k2);
	T k4 = func(t + h, x + h * k3);

	return x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

ostream& operator<<(ostream& out, vector<double>& vec) {
	for (double el : vec)
		out << el << "  ";
	return out;
}


// TODO делает лишний шаг
// Метод Рунге - Кутты с коррекцией шага
// func - правая часть дифференциального уравнения, первый аргумент - время t
// interval - интервал интегрирования
// initVal - начальное значение
// eps - требуемая точность
// initH - начальный шаг
// T может быть числом или вектором
template<typename T>
pair<vector<double>, vector<T>> RK4auto(function<T(double, T)> func,
	pair<double, double> interval, T initVal, double eps, double initH) {

	// порядок метода
	const double p = 4.0;
	// коэффициент при вычислении ошибки
	const double c = 1.0 / (pow(2.0, p) - 1.0);

	double t = interval.first;
	double h = initH;
	T x = initVal;

	vector<double> tt;
	tt.reserve((interval.second - interval.first) / initH);
	vector<T> xx;
	xx.reserve((interval.second - interval.first) / initH);
	tt.push_back(t);
	xx.push_back(initVal);


	while (t < interval.second) {
		T newX1 = RK4step(func, t, x, h);
		T newX2 = RK4step(func, t + h / 2, RK4step(func, t, x, h / 2), h / 2);

		double err = norm(c * (newX1 + (-1.0) * newX2));
		if (err < eps) {
			t += h;
			x = newX1;
			tt.push_back(t);
			xx.push_back(x);

			//if (err < eps / 6.0)
			if (err < eps / 6.0 && t + 2 * h <= interval.second - h * 1e-7)
				h *= 2.0;
		}
		else {
			h /= 2.0;
		}

		if (t + h > interval.second - h * 1e-7)
			h = interval.second - t;
	}

	return make_pair(tt, xx);
}

// Метод Рунге - Кутты с коррекцией шага и выводом в файл
// func - правая часть дифференциального уравнения, первый аргумент - время t
// interval - интервал интегрирования
// initVal - начальное значение
// eps - требуемая точность
// initH - начальный шаг
// exact - точное решение
// out - поток вывода
// T может быть числом или вектором
template<typename T>
void RK4auto(function<T(double, T)> func,
	pair<double, double> interval,
	T initVal,
	double eps,
	double initH,
	function<T(double)> exact,
	ostream& out) {

	// порядок метода
	const double p = 4.0;
	// коэффициент при вычислении ошибки
	const double c = 1.0 / (pow(2.0, p) - 1.0);
	//const double c = pow(2.0, p) / (pow(2.0, p) - 1.0);

	double t = interval.first;
	double h = initH;
	T x = initVal;


	while (t < interval.second) {
		T newX1 = RK4step(func, t, x, h);
		T newX2 = RK4step(func, t + h / 2, RK4step(func, t, x, h / 2), h / 2);

		double err = c * norm(newX1 + (-1.0) * newX2);
		if (err < eps) {
			out << t << "  " << x << "  " << norm(exact(t) + (-1.0) * x) << "  " << err << "  " << h << endl;

			t += h;
			//x = newX1;
			x = newX2;

			//if (err < eps / 6.0)
			if (err < eps / 6.0 && t + 2 * h <= interval.second - h * 1e-7)
				h *= 2.0;
		}
		else {
			h /= 2.0;
		}

		if (t + h > interval.second - h * 1e-7)
			h = interval.second - t;
	}
}


int main()
{
	setlocale(LC_ALL, "Russian");
	double eps1 = 1e-2;
	double eps2 = 1e-4;
	double eps3 = 1e-7;
	double h = 0.01;
	vector<double> initvalue{ 0,0 };
	vector<function<double(double, double)>> functions_mx{ mx1,mx2 };




	//Пример 1
	/*auto result_Euler_1 = Euler_explicit<vector<double>>(f1, pair<double, double>{ 0, 1.0 }, initvalue, h);
	inFile(result_Euler_1.second, "Эйлер.txt");
	print_vec(result_Euler_1);*/


	/*auto result_RK2_1 = RK2<vector<double>>(f1, pair<double, double>{0.0, 1.0}, initvalue, h);
	inFile(result_RK2_1.second, "Рунге-Кутта2.txt");
	print_vec(result_RK2_1);*/

	/*auto result_RK4_1 = RK4<vector<double>>(f1, pair<double, double>{0.0, 1.0}, initvalue, h);
	inFile(result_RK4_1.second, "Рунге-Кутта4.txt");
	print_vec(result_RK4_1);


	auto result_Adams_1 = Adams<vector<double>>(f1, pair<double, double>{0.0, 1.0}, initvalue, h);
	inFile(result_Adams_1.second, "Адамс.txt");
	print_vec(result_Adams_1);

	auto result_PEC_1 = PEC<vector<double>>(f1, pair<double, double>{0.0, 1.0}, initvalue, h);
	inFile(result_PEC_1.second, "Прогноз-коррекция.txt");
	print_vec(result_PEC_1);*/




	// Пример 2
	//auto result_Euler_2 = Euler_explicit<vector<double>>(f2, pair<double, double>{ 0, 1.0 }, initvalue, h);
	////print_vec(result_Euler2);
	//inFile(result_Euler2.second, "Эйлер.txt");

	/*auto result_RK2_2 = RK2<vector<double>>(f2, pair<double, double>{0.0, 1.0}, initvalue, h);
	inFile(result_RK2_2.second, "Рунге-Кутта2.txt");
	print_vec(result_RK2_2);*/

	//auto result_RK4_2 = RK4<vector<double>>(f2, pair<double, double>{0.0, 1.0}, initvalue, h);
	////print_vec(result_RK4_2);
	//inFile(result_RK4_2.second, "Рунге-Кутта.txt");

	//auto result_Adams_2 = Adams<vector<double>>(f2, pair<double, double>{0.0, 1.0}, initvalue, h);
	////print_vec(result_Adams_2);
	//inFile(result_Adams_2.second, "Адамс.txt");

	//auto result_PEC_2 = PEC<vector<double>>(f2, pair<double, double>{0.0, 1.0}, initvalue, h);
	////print_vec(result_PEC_2);
	//inFile(result_PEC_2.second, "Прогноз-коррекция.txt");



	//Пример 3

	//vector<double> initvalue3{ 4,4,4 };
	// 
	// 
	//auto result_Euler3 = Euler_explicit<vector<double>>(f3, pair<double, double>{ 0, 1.0 }, initvalue3, h);
	////print_vec(result_Euler3);
	//inFile(result_Euler3.second, "Эйлер.txt");

	/*auto result_RK2_3 = RK2<vector<double>>(f3, pair<double, double>{0.0, 1.0}, initvalue3, h);
	inFile(result_RK2_3.second, "Рунге-Кутта2.txt");
	print_vec(result_RK2_3);*/

	//auto result_RK4_3 = RK4<vector<double>>(f3, pair<double, double>{0.0, 1.0}, initvalue3, h);
	////print_vec(result_RK4_3);
	//inFile(result_RK4_3.second, "Рунге-Кутта.txt");

	//auto result_Adams_3 = Adams<vector<double>>(f3, pair<double, double>{0.0, 1.0}, initvalue3, h);
	////print_vec(result_Adams_3);
	//inFile(result_Adams_3.second, "Адамс.txt");

	//auto result_PEC_3 = PEC<vector<double>>(f3, pair<double, double>{0.0, 1.0}, initvalue3, h);
	////print_vec(result_PEC_3);
	//inFile(result_PEC_3.second, "Прогноз-коррекция.txt");

	//order(exactTest, fTest, pair<double, double>{ 0, 1 }, 1, 0.1, 12);




	/////////////////
	// Часть 2. Оценка порядка
	auto interv = make_pair(1.0, 5.0);
	auto f = labTestDiff;
	auto exact = labTest;
	auto diff = labTestDiffDiff;

	//order(exactTest, fTest, pair<double, double>{ 0, 1 }, 1, M_RK2, 10);
	//order(exactTest, fTest, pair<double, double>{ 0, 1 }, 1, M_RK2);

	//order(exact, f, interv, 0, M_EULER_EXPL);
	//order(exact, f, interv, 0.0, M_EULER_EXPL, 10);

	// сокращенный набор таблиц (только с помощью эталонного решения и в узле
	// близком к t0)
	cout << "--- Сокращенный набор таблиц ---" << endl;
	size_t point = 10;
	cout << "Явный метод Эйлера:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_EULER_EXPL, 10);
	cout << endl;
	cout << "Неявный метод Эйлера:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_EULER_IMPL, 10);
	// TODO не работает
	//cout << endl;
	//cout << "Симметричная схема:" << endl;
	//order(labTest, labTestDiff, diff, interv, 0.0, M_SYMM, 10);
	cout << endl;
	cout << "Метод Рунге-Кутты 2 порядка:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_RK2, 10);
	cout << endl;
	cout << "Метод Рунге-Кутты 4 порядка:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_RK4, 10);
	cout << endl;
	cout << "Метод Адамса-Башфорта:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_ADAMS, 10);
	cout << endl;
	cout << "Метод прогноза и коррекции:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_PEC, 10);

	/*cout << endl << endl << "------------------------------------------------------" << endl;
	cout << "--- Полный набор таблиц для РК4 и Адамса 4 порядка ---" << endl;
	cout << "------------------------------------------------------" << endl;
	// эталонное решение
	cout << "РК4, эталонное решение, t = 10 tau:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_RK4, 10);
	cout << endl;
	cout << "РК4, эталонное решение, t = T:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_RK4);
	// правило Эйткена
	cout << endl;
	cout << "РК4, правило Эйткена, t = 10 tau:" << endl;
	eitken(labTestDiff, diff, interv, 0.0, M_RK4, 10);
	cout << endl;
	cout << "РК4, правило Эйткена, t = T:" << endl;
	eitken(labTestDiff, diff, interv, 0.0, M_RK4);

	cout << endl;
	cout << "Адамс, эталонное решение, t = 10 tau:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_ADAMS, 10);
	cout << endl;
	cout << "Адамс, эталонное решение, t = T:" << endl;
	order(labTest, labTestDiff, diff, interv, 0.0, M_ADAMS);
	// правило Эйткена
	cout << endl;
	cout << "Адамс, правило Эйткена, t = 10 tau:" << endl;
	eitken(labTestDiff, diff, interv, 0.0, M_ADAMS, 10);
	cout << endl;
	cout << "Адамс, правило Эйткена, t = T:" << endl;
	eitken(labTestDiff, diff, interv, 0.0, M_ADAMS);*/


	/////////////////
	// Часть 2.5. Иллюстрация правила Рунге
	ofstream out("RK4auto.txt");
	double eps = 1e-6;
	//RK4auto<double>(fTest, pair<double, double>{ 0, 1 }, 1.0, eps, 0.1, exactTest, out);
	RK4auto<double>(labTestDiff, pair<double, double>{ 1, 5 }, 0.0, eps, 0.05, labTest, out);
	out.close();
	/////////////////



	/*
	/////////////////
	// Часть 3. Задача о маятнике
	// уравнения математического маятника
	auto pendulum = mx();
	auto pendDiff = mxDiff();

	vector<double> init{ 1.0, 1.0 };
	auto interval = make_pair(0.0, 10.0);
	double tau = 0.01;

	// явный метод Эйлера
	auto Euler_expl_mx = Euler_explicit<vector<double>>(pendulum, interval, init, tau);
	inFile(Euler_expl_mx.second, "Маятник_явнЭйлер.txt");
	// неявный метод Эйлера
	auto Euler_impl_mx = EulerSystem_implicit<vector<double>>(pendulum, pendDiff, interval, init, tau);
	inFile(Euler_impl_mx.second, "Маятник_неявнЭйлер.txt");
	// симметричная схема
	auto Symm_mx = SymmSystem<vector<double>>(pendulum, pendDiff, interval, init, tau);
	inFile(Symm_mx.second, "Маятник_симметричная.txt");
	// Рунге-Кутты 2 порядка
	auto RK2_mx = RK2<vector<double>>(pendulum, interval, init, tau);
	inFile(RK2_mx.second, "Маятник_РК2.txt");



	auto result_Euler_explicit_mx = Euler_explicit_mx();
	inFile(result_Euler_explicit_mx, "Эйлер_явный_mx.txt");

	auto result_Euler_implicit_mx = Euler_implicit_mx();
	inFile(result_Euler_implicit_mx, "Эйлер_неявный_mx.txt");
		/////////////////
		pend_err();
		*/




		/*
		/////////////////
		// Часть 4. Фазовый портрет
		vector<double> init1{0.1, 0.1};
		vector<double> init2{ -0.1, -0.1 };
		pair<double, double> inter4{ 0.0, 200.0 };
		auto RK4_var1_1 = RK4<vector<double>>(variant1, inter4, init1, 0.1);
		auto RK4_var1_2 = RK4<vector<double>>(variant1, inter4, init2, 0.1);
		//auto RK4_var1_3 = RK4<vector<double>>(variant1, inter4, vector{-1.7, 2.3}, 0.1);
		auto Symm_var1_1 = SymmSystem<vector<double>>(variant1, variant1Diff, inter4, init1, 0.1);
		auto Symm_var1_2 = SymmSystem<vector<double>>(variant1, variant1Diff, inter4, init2, 0.1);
		auto EulerImpl_var1_1 = EulerSystem_implicit<vector<double>>(variant1, variant1Diff, inter4, init1, 0.1);
		auto EulerImpl_var1_2 = EulerSystem_implicit<vector<double>>(variant1, variant1Diff, inter4, init2, 0.1);

		inFile(RK4_var1_1.second, "part4_rk4_1.txt");
		inFile(RK4_var1_2.second, "part4_rk4_2.txt");
		//inFile(RK4_var1_3.second, "part4_rk4_3.txt");
		inFile(Symm_var1_1.second, "part4_symm_1.txt");
		inFile(Symm_var1_2.second, "part4_symm_2.txt");
		inFile(EulerImpl_var1_1.second, "part4_impl_1.txt");
		inFile(EulerImpl_var1_2.second, "part4_impl_2.txt");
		/////////////////
		*/

}