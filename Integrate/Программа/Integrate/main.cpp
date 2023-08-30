#include <iostream>
#include <vector>
#include <array>
#include <tuple>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include <functional>

using namespace std;

// TODO:
// 4. Решение сингулярных уравнений

// запись в файл вектора
template<typename T>
ostream& operator<<(ostream& out, const vector<T>& v) {
	for (size_t i = 0; i < v.size(); i++)
		out << v[i] << "    ";
	return out;
}

// запись в файл сетки с решением
template<typename T>
ostream& operator<<(ostream& out, const pair<vector<T>, vector<T>>& p) {
	out << p.first << endl << p.second;
	return out;
}


/////////////////////////////////////////////////
// Метод Гаусса
/////

//считает сумму, чего надо вычесть
template<typename T>
T sum(size_t i, size_t N, const vector<vector<T>>& vec_A, const vector<T>& vec_X)
{
	T sum = 0;
	int k = 0;
	for (size_t j = i + 1; j < N; j++)
	{
		sum += vec_A[i][j] * vec_X[k];
		k++;
	}
	return sum;
}


//обратный ход
// сложность O(n^2 / 2)
template<typename T>
void reverse_steps(int N, const vector<vector<T>>& vec_A, const vector<T>& vec_B, vector<T>& vec_X)
{
	T x = vec_B[N - 1] / vec_A[N - 1][N - 1];
	vec_X.insert(vec_X.begin(), x);

	// нельзя менять на size_t
	for (int i = N - 2; i > -1; i--) {
		x = (vec_B[i] - sum(i, N, vec_A, vec_X)) / vec_A[i][i];
		vec_X.insert(vec_X.begin(), x);
	}
}

//перестановка строк (частичный выбор)
template<typename T>
void main_str(size_t i, vector<vector<T>>& vec_A, vector<T>& vec_B)
{
	size_t k = i;
	// Ищем максимум
	T max = abs(vec_A[i][i]);
	for (size_t j = i + 1; j < vec_B.size(); j++) {
		// TODO стоит проверять на ноль (но что использовать в качетсве нуля для чисел с плавающей запятой)
		if (abs(vec_A[j][i]) >= max) {
			max = abs(vec_A[j][i]);
			k = j;
		}
	}

	// тут поменял
	if (max <= numeric_limits<T>::epsilon() * 10)
		throw "Вырожденная матрица";

	// Меняем строки
	swap(vec_A[k], vec_A[i]);
	swap(vec_B[k], vec_B[i]);
}




//Метод Гаусса
template<typename T>
vector<T> gauss(vector<vector<T>>& A, vector<T>& B) {
	// Поверим, что матрица A квадратная и имеет такой же размер, что и у B
	size_t N = A.size();
	vector<T> X;
	// тут поменял
	T c_ik;
	for (size_t k = 0; k < N - 1; ++k) {
		main_str(k, A, B);//строку поднимаем с max-элементом
		for (size_t i = k + 1; i < N; ++i) {
			//for (size_t i = k; i < N; ++i) {
			c_ik = A[i][k] / A[k][k];

			B[i] -= c_ik * B[k];
			for (size_t j = k + 1; j < N; ++j)
				A[i][j] -= c_ik * A[k][j];
		}
	}

	// вроде бы это еще проверяет ранг (то, что ранг матрицы A равен рангу расширенной матрицы)
	if (abs(A[N - 1][N - 1]) <= numeric_limits<T>::epsilon() * 10)
		throw "Вырожденная матрица";

	reverse_steps(N, A, B, X); // обратный ход
	return X;
}
/////////////////////////////////////////////////

// Метод квадратур. Внутри используется квадратурная формула трапеций (имеет
// второй порядок) и метод Гаусса.
// K(x, s) - ядро интегрального уравнения
// f(x) - правая часть
// [a, b] - отрезок, на котором ищется решение
// N - количество "отрезков" (N+1 точка)
// возвращает пару "сетка"-"значения в узлах"
pair<vector<double>, vector<double>> quadratureMethod(
	function<double(double, double)> K, function<double(double)> f,
	double a, double b, size_t N
) {
	double h = (b - a) / N;

	// создаем сетку
	vector<double> xx;
	xx.reserve(N + 1);
	for (size_t i = 0; i <= N; ++i)
		xx.push_back(a + h * i);

	// правая часть СЛАУ
	vector<double> ff;
	ff.reserve(N + 1);
	for (size_t i = 0; i <= N; ++i)
		ff.push_back(f(xx[i]));

	// матрица СЛАУ
	vector<vector<double>> A;
	for (size_t i = 0; i <= N; ++i) {
		vector<double> Arow;
		Arow.reserve(N + 1);
		Arow.push_back(-h / 2 * K(xx[i], xx[0]));
		for (size_t j = 1; j <= N - 1; ++j) {
			Arow.push_back(-h * K(xx[i], xx[j]));
		}
		Arow.push_back(-h / 2 * K(xx[i], xx[N]));
		Arow[i] += 1.0;
		A.push_back(Arow);
	}

	// решаем получившуюся СЛАУ
	vector<double> res = gauss(A, ff);

	return make_pair(xx, res);
}

// отдельно можно задать lambda
pair<vector<double>, vector<double>> quadratureMethod(
	function<double(double, double)> K, double lambda, function<double(double)> f,
	double a, double b, size_t N
) {
	auto newK = [K, lambda](double x, double y) {
		return lambda * K(x, y);
	};

	return quadratureMethod(newK, f, a, b, N);
}

// Метод простой итерации. Внутри используется квадратурная формула трапеций
// (имеет второй порядок)
// K(x, s) - ядро интегрального уравнения
// f(x) - правая часть
// [a, b] - отрезок, на котором ищется решение
// N - количество "отрезков" (N+1 точка)
// eps - требуемая точность
// возвращает пару "сетка"-"значения в узлах"
pair<vector<double>, vector<double>> simpleIterationMethod(
	function<double(double, double)> K, function<double(double)> f,
	double a, double b, size_t N, double eps
) {
	double h = (b - a) / N;

	// создаем сетку
	vector<double> xx;
	xx.reserve(N + 1);
	for (size_t i = 0; i <= N; ++i)
		xx.push_back(a + h * i);

	// вектор правой части
	vector<double> ff;
	ff.reserve(N + 1);
	for (size_t i = 0; i <= N; ++i)
		ff.push_back(f(xx[i]));

	// вспомогательная матрица (зашита формула трапеций)
	vector<vector<double>> A;
	for (size_t i = 0; i <= N; ++i) {
		vector<double> Arow;
		Arow.reserve(N + 1);
		Arow.push_back(h / 2 * K(xx[i], xx[0]));
		for (size_t j = 1; j <= N - 1; ++j) {
			Arow.push_back(h * K(xx[i], xx[j]));
		}
		Arow.push_back(h / 2 * K(xx[i], xx[N]));
		A.push_back(Arow);
	}

	// в качестве начального приближения берем правую часть f
	vector<double> u0 = ff, u1(N + 1);

	double crit = eps + 1;
	// количество итераций (так, для справки)
	size_t k = 0;
	// собственно сам итерационный процесс
	while (crit >= eps) {

		u1 = ff;

		// перемножаем матрицу A на текущее приближение
		for (size_t i = 0; i <= N; ++i) {
			for (size_t j = 0; j <= N; ++j)
				u1[i] += A[i][j] * u0[j];
		}

		// считаем критерий выхода по ранице двух слоев (вроде бы он же по невязке)
		crit = 0.0;
		for (size_t i = 0; i <= N; ++i)
			crit = max(crit, abs(u0[i] - u1[i]));

		u0 = u1;

		++k;

		/////////////////////////////
		// для теста
		//auto u_exact = [](double x) { return cos(x) - 2.0 / M_PI * sin(x); };
		//double norm = 0.0;
		//for (size_t i = 0; i < xx.size(); ++i)
		//	norm = max(norm, abs(u_exact(xx[i]) - u0[i]));
		//cout << "k = " << k << ":    " << norm << endl;
	}

	cout << "Кол-во итераций: " << k << endl;

	return make_pair(xx, u1);
}

// отдельно можно задать lambda
pair<vector<double>, vector<double>> simpleIterationMethod(
	function<double(double, double)> K, double lambda, function<double(double)> f,
	double a, double b, size_t N, double eps
) {
	auto newK = [K, lambda](double x, double y) {
		return lambda * K(x, y);
	};

	return simpleIterationMethod(newK, f, a, b, N, eps);
}

// phi(x)
// psi(s)
pair<vector<double>, vector<double>> degenerateNucleiMethod(
	vector<function<double(double)>> phi, vector<function<double(double)>> psi,
	double lambda, function<double(double)> f, double a, double b, size_t N
) {
	size_t m = phi.size();
	if (m != psi.size())
		throw "размерности не совпадают!";

	double h = (b - a) / N;

	// создаем сетку
	vector<double> xx;
	xx.reserve(N + 1);
	for (size_t i = 0; i <= N; ++i)
		xx.push_back(a + h * i);

	vector<vector<double>> alpha(m, vector<double>(m));
	vector<double> betha(m);

	// функция для интегрирования (методом трапеций)
	auto integrate = [&xx, N, h](
		function<double(double)> func
		) {
			double res = 0.0;
			// граничные точки
			res += h / 2.0 * (func(xx[0]) + func(xx[N]));
			for (size_t i = 1; i < N; ++i)
				res += h * func(xx[i]);
			return res;
	};

	// считаем коэффициенты
	for (size_t i = 0; i < m; ++i) {
		auto psi_i = psi[i];
		auto bethaFunc = [&psi_i, &f](double x) {
			return psi_i(x) * f(x);
		};
		betha[i] = integrate(bethaFunc);
		for (size_t j = 0; j < m; ++j) {
			auto phi_j = phi[j];
			auto alphaFunc = [&psi_i, &phi_j](double x) {
				return psi_i(x) * phi_j(x);
			};

			alpha[i][j] = integrate(alphaFunc);
		}
	}

	// для матрицы СЛАУ переиспользуем матрицу alpha
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < m; ++j)
			alpha[i][j] *= -lambda;

	for (size_t i = 0; i < m; ++i)
		alpha[i][i] += 1.0;

	// решаем СЛАУ
	vector<double> C = gauss(alpha, betha);

	vector<double> u;
	u.reserve(N + 1);
	// проецируем решение на сетку
	for (size_t i = 0; i <= N; ++i) {
		double x = xx[i];
		double sum = 0.0;
		for (size_t i = 0; i < m; ++i)
			sum += C[i] * phi[i](x);
		u.push_back(f(x) + lambda * sum);
	}

	return make_pair(xx, u);
}

void test1() {
	auto f = [](double x) { return sin(M_PI * x); };
	auto K = [](double x, double t) { return 0.5; };

	double eps = 1e-6;

	//auto res = quadratureMethod(K, f, 0.0, 1.0, 6);
	auto res = simpleIterationMethod(K, f, 0.0, 1.0, 6, eps);

	ofstream out("test1.txt");
	out << res;
	out.close();
}

void test2() {
	auto f = [](double x) { return 1.0 + x * x; };
	auto K = [](double x, double t) { return x / (1 + t * t); };

	double eps = 1e-6;

	//auto res = quadratureMethod(K, f, 0.0, 1.0, 2);
	auto res = simpleIterationMethod(K, f, 0.0, 1.0, 2, eps);

	ofstream out("test2.txt");
	out << res;
	out.close();
}

// тест для проверки погрешности метода простой итерации в зависимости от итерации
void testSimpleIter() {
	auto K = [](double x, double s) { return s * sin(x) / 2.0 / M_PI; };
	auto f = [](double x) { return cos(x); };
	double eps = 1e-7;

	auto res = simpleIterationMethod(K, f, 0, M_PI, 50, eps);
}

// тест для проверки влияния количества членов на погрешность
// N - число узлов
void testDegen(size_t m, size_t N) {
	vector<function<double(double)>> phi, psi;
	phi.reserve(m); psi.reserve(m);
	phi.push_back(
		[](double x) { return x * x; }
	);
	psi.push_back(
		[](double s) {return 1.0; }
	);

	for (size_t i = 1; i < m; ++i) {
		phi.push_back(
			[i](double x) { return pow(x, 2 * (i + 1)); }
		);
		double fact = tgamma(i + 1);
		psi.push_back(
			[i, fact](double s) { return pow(s, 4 * i) / fact; }
		);
	}

	auto f = [](double x) { return x * x * x - exp(x * x) + 1.0; };

	auto res = degenerateNucleiMethod(phi, psi, 4.0, f, 0.0, 1.0, N);

	auto exact = [](double x) { return x * x * x; };

	double norm = 0.0;
	for (size_t i = 0; i < res.first.size(); ++i)
		norm = max(norm, abs(exact(res.first[i]) - res.second[i]));

	cout << "m = " << m << ":    " << norm << endl;
}

// TODO вбить примеры из методички
// пример 1 из методички
void task1_1() {
	auto K = [](double x, double t) {
		return 0.5 * (1 - x * cos(x * t));
	};
	//auto f = [](double x) { return 0.5 * (1.0 + sin(x)); };
	auto f = [](double x) { return x * x + sqrt(x); };
	//double a = 0.0, b = 1.0;
	double a = 0.1, b = 1.0;

	auto resQuad = quadratureMethod(K, f, a, b, 12);
	double eps = 1e-4;
	//auto resIter = simpleIterationMethod(K, f, a, b, 10, eps);

	ofstream out("task1.txt");
	out << resQuad;
	//out << resIter;
	out.close();
}

///////////////////////////////////////////////////////////
// сингулярные уравнения
///////////////////////////////////////////////////////////

using vec2d = array<double, 2>;

// скалярное произведение двух двумерных векторов
double operator*(const vec2d& l, const vec2d& r) {
	return l[0] * r[0] + l[1] * r[1];
}

// ядро интегрального уравнения
vec2d Q(vec2d& r, vec2d& rho) {
	double denom = 2.0 * M_PI * ((r[0] - rho[0]) * (r[0] - rho[0]) +
		(r[1] - rho[1]) * (r[1] - rho[1]));

	return vec2d{ (rho[1] - r[1]) / denom, (r[0] - rho[0]) / denom };
}

// i-я контрольная точка
vec2d k_i(size_t N, size_t i) {
	return vec2d{ cos(2 * M_PI * (i - 0.5) / N), sin(2 * M_PI * (i - 0.5) / N) };
}

// j-я расчетная точка
vec2d c_j(size_t N, size_t j) {
	return vec2d{ cos(2 * M_PI * (j - 1.0) / N), sin(2 * M_PI * (j - 1.0) / N) };
}

// номаль в контрольной точке k_i
vec2d n_i(size_t N, size_t i) {
	return k_i(N, i);
}

pair<vector<vec2d>, vector<double>> singular(
	size_t N,
	function<double(vec2d)> f
) {
	// создаем расчетные точки c_i, контрольные точки k_i (они же нормали)
	vector<vec2d> c, k;
	c.reserve(N); k.reserve(N);

	for (size_t i = 0; i < N; ++i) {
		c.push_back(c_j(N, i + 1)); // в формуле нумерация с 1
		k.push_back(k_i(N, i + 1));
	}

	// длина дуг
	double delta_l = 2 * M_PI / N;

	// создаем матрицу СЛАУ и правую часть

	vector<vector<double>> A;
	A.reserve(N + 1);

	vector<double> F;
	F.reserve(N + 1);

	for (size_t i = 0; i < N; ++i) {
		vector<double> Arow;
		Arow.reserve(N + 1);

		for (size_t j = 0; j < N; ++j)
			Arow.push_back(
				// n(k[i]) = k[i]
				k[i] * Q(k[i], c[j]) * delta_l
			);

		Arow.push_back(1.0);

		A.push_back(Arow);

		F.push_back(f(k[i]));
	}

	// заполняем последнюю строку
	vector<double> Arow(N + 1, delta_l);
	Arow[N] = 0.0;
	A.push_back(Arow);

	F.push_back(0.0);

	// решаем СЛАУ
	vector<double> g_R = gauss(A, F);

	//cout << "R = " << g_R[N] << endl;

	// вектор значений - это всё кроме последнего элемента R
	g_R.pop_back();
	return make_pair(c, g_R);
}

// метод без R (срабатывает исключение на вырожденную матрицу)
pair<vector<vec2d>, vector<double>> singular2(
	size_t N,
	function<double(vec2d)> f
) {
	// создаем расчетные точки c_i, контрольные точки k_i (они же нормали)
	vector<vec2d> c, k;
	c.reserve(N); k.reserve(N);

	for (size_t i = 0; i < N; ++i) {
		c.push_back(c_j(N, i + 1)); // в формуле нумерация с 1
		k.push_back(k_i(N, i + 1));
	}

	// длина дуг
	double delta_l = 2 * M_PI / N;

	// создаем матрицу СЛАУ и правую часть

	vector<vector<double>> A;
	A.reserve(N);

	vector<double> F;
	F.reserve(N);

	for (size_t i = 0; i < N; ++i) {
		vector<double> Arow;
		Arow.reserve(N);

		for (size_t j = 0; j < N; ++j)
			Arow.push_back(
				// n(k[i]) = k[i]
				k[i] * Q(k[i], c[j]) * delta_l
			);

		A.push_back(Arow);

		F.push_back(f(k[i]));
	}

	// решаем СЛАУ
	vector<double> g = gauss(A, F);

	cout << "R = " << 0 << endl;

	return make_pair(c, g);
}

// запись решения сингулярного уравнения в файл (вместо координаты оставляет угол)
ostream& operator<<(ostream& out, const pair<vector<vec2d>, vector<double>>& v) {
	for (size_t i = 0; i < v.first.size(); i++) {
		//double phi = atan2(v.first[i][1], v.first[i][0]) - M_PI;
		double phi = atan2(v.first[i][1], v.first[i][0]);

		if (phi >= 0)
			out << phi << "    ";
		else
			out << (phi + 2.0 * M_PI) << "    ";
	}
	out << endl;
	out << v.second;
	return out;
}

// правая часть для нечетных вариантов
function<double(vec2d)> f_odd(size_t N) {
	return [N](vec2d arg) {
		return sin((N + 1) / 2 * atan2(arg[1], arg[0]));
	};
}

// правая часть для четных вариантов
function<double(vec2d)> f_even(size_t N) {
	return [N](vec2d arg) {
		return cos(N / 2.0 * atan2(arg[1], arg[0]));
	};
}

void singular_test(size_t N = 100) {
	//auto f = f_odd(1);
	auto f = f_odd(13);
	//auto f = f_odd(7);
	//auto f = f_even(2);

	auto res = singular(N, f);
	//auto res = singular2(N, f);


	ofstream out("singular.txt");
	out << res;
	out.close();
}

// сравнение с точным решением
double singular_test_exact(size_t N) {
	//auto f = f_odd(1);
	size_t var = 13;
	//size_t var = 1;
	auto f = f_odd(var);
	//auto f = f_odd(7);
	//auto f = f_even(2);

	auto res = singular(N, f);
	//auto res = singular2(N, f);
	auto exact = [var](double phi) {
		return -2.0 * cos((var + 1) / 2 * phi);
	};

	ofstream out("singular.txt");
	out << res;
	out.close();

	double maxErr = 0.0;

	for (size_t i = 0; i < res.first.size(); ++i) {
		double phi = atan2(res.first[i][1], res.first[i][0]);
		double exact_val = exact(phi);
		double num_val = res.second[i];

		maxErr = max(maxErr, abs(exact_val - num_val));
	}

	return maxErr;
}


int main()
{
	setlocale(LC_ALL, "Russian");
	//test1();
	//test2();

	//for (size_t m = 1; m <= 12; ++m)
	//	testDegen(m, 200);

	task1_1();

	//testSimpleIter();

	//singular_test(1000);

	// тест на R
	/*for (size_t N = 10; N <= 100; N += 10) {
		cout << "N = " << N << ":" << endl;
		singular_test(N);
		cout << endl;
	}*/

	//cout << singular_test_exact(30);

	// тест на нужное кол-во узлов для точного решения
	//for (size_t N = 3; N <= 12; N++) {
	//	cout << "N = " << N << ": " << singular_test_exact(N) << endl;
	//}
}


exp(5*pow())