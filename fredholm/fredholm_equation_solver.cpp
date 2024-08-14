#include <functional>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>


using Coefficients = std::vector<double>;

Coefficients four_bode_rule(double a, double b) {
	Coefficients coeff(5, 0);

	double h = (b - a) / 4;

	coeff[0] = h * 14. / 45;
	coeff[1] = h * 64. / 45;
	coeff[2] = h * 24. / 45;
	coeff[3] = h * 64. / 45;
	coeff[4] = h * 14. / 45;

	return coeff;
}

Coefficients bode_rule(double a, double b, double h) {

	// N = 5 * r - (r - 1)
	// N = 5 * r - r + 1
	// N = 4 * r + 1
	// (N - 1) / 4

	int n = (b - a) / h;

	int r = (n - 1 + 3) / 4;

	if (r <= 0) {
		r = 1;
	}

	n = 4 * r + 1;

	h = (b - a) / (n - 1); 

	std::vector<double> coeff(n, 0);

	for (int i = 0; i < n - 1; i += 4) {
		coeff[i] += h * 14. / 45;
		coeff[i + 1] = h * 64. / 45;
		coeff[i + 2] = h * 24. / 45;
		coeff[i + 3] = coeff[i + 1];
		coeff[i + 4] = h * 14. / 45;
	}

	return coeff;
}

double integrate(std::function<double(double)> f, double a, double b, double h) {
	Coefficients coefficients = bode_rule(a, b, h);
	double step = (b - a) / (coefficients.size() - 1);
	double sum {0.};
	for (int i = 0; i < coefficients.size(); ++i) {
		sum += coefficients[i] * f(a + step * i);
	}
	return sum;
}

double adaptive_integrate(std::function<double(double)> f, double a, double b, double ea, double eo) {
	double l = a;
	double r = b;
	double I = 0;
	double E = 0;
	int k = 0;
	int kmax = 0;

	while (l < r) {
		double h = r - l;
		double I_h = integrate(f, l, r, h);
		double I_h_2 = integrate(f, l, (l + r) / 2, h / 2) + integrate(f, (l + r) / 2, r, h / 2);

		double d = std::abs(I_h - I_h_2) / 31;

		if (d < std::max(ea, eo * std::abs(I)) * h / (b - a)) {
			k = 0;
			l = r;
			r = b;
			I += I_h_2 + d;
			E += d;
		} else {
			if (k > kmax) {
				std::cout << "unusual point" << std::endl;
				I += I_h_2 + d;
				E += d;
				l = r;
				r = b;
				continue;
			}
			r = (l + r) / 2;
			k += 1;
		}
	}
	return I;
}

void solve(std::vector<double> & res, const std::vector<std::vector<double>>& mtx) {

	if (mtx.size() == 0) {
		throw std::runtime_error("solve: mtx.size() == 0");
	}
	if (mtx[0].size() == 0) {
		throw std::runtime_error("solve: mtx[0].size() == 0");
	}

	if (mtx.size() > mtx[0].size()) {
		throw std::runtime_error("solve: mtx.size() > mtx[0].size()");
	}

	std::vector<std::vector<double>> mtx_copy {mtx};

	size_t n = mtx_copy.size();
	size_t m = mtx_copy[0].size();

	if (m == 1) {
		throw std::runtime_error("solve: m == 1");
	}

	auto find_non_zero_row = [&mtx_copy](size_t i) -> size_t {
		// double eps = 1e-8;
		for (size_t j = i; j < mtx_copy.size(); ++j) {
			if (std::fabs(mtx_copy[j][i]) > 0) {
				return j;
			}
		}
		return std::numeric_limits<size_t>::max();
	};

	auto swap_rows = [&mtx_copy](size_t i, size_t j) {
		std::swap(mtx_copy[i], mtx_copy[j]);
	};

	auto normalize_row = [&mtx_copy](size_t i) {
		double diviser = mtx_copy[i][i];
		for (size_t j = i; j < mtx_copy[i].size(); ++j) {
			mtx_copy[i][j] /= diviser;
		}
	};

	auto add_row = [&mtx_copy](size_t dst, const size_t src, double coeff) {
		for (unsigned j = 0; j < mtx_copy[0].size(); ++j) {
			mtx_copy[dst][j] += mtx_copy[src][j] * coeff;
		}
	};

	for (size_t col_idx = 0; col_idx < m - 1; ++col_idx) {
		size_t non_zero_row = find_non_zero_row(col_idx);

		if (non_zero_row == std::numeric_limits<size_t>::max()) {
			continue;
		}

		swap_rows(col_idx, non_zero_row);
		normalize_row(col_idx);
		for (size_t row_idx = 0; row_idx < n; ++row_idx) {
			 if (row_idx == col_idx) {
			 	continue;
			 }
			 double coeff = -mtx_copy[row_idx][col_idx];
			 add_row(row_idx, col_idx, coeff);
		}
	}

	res.resize(m - 1);

	for (size_t row_idx = 0; row_idx < m - 1; ++row_idx) {
		res[row_idx] = mtx_copy[row_idx].back();
	}
}

std::function<double(double)> solve_fredholm_equation(std::function<double(double, double)> K, 
                                                      std::function<double(double)> f, double a, double b, double eps) {

	double h = (b - a) / 2;

	std::function<double(double)> prev_func = [](double x){return 0;};
	std::function<double(double)> next_func = [](double x){return 0;};

	while(true) {
		prev_func = next_func;
		h /= 2;
		Coefficients A = bode_rule(a, b, h);

		size_t n = A.size();
		h = (b - a) / (n - 1);

		std::vector<std::vector<double>> equations(n, std::vector<double>(n + 1));

		for (size_t i = 0; i < n; ++i) {
			double s_i = a + i * h;
			for (size_t j = 0; j < n; ++j) {
				double s_j = a + j * h;
				equations[i][j] = -A[j] * K(s_i, s_j);
			}
			equations[i][n] = f(s_i);
			equations[i][i] += 1;
		}

		std::vector<double> U;

		solve(U, equations);

		next_func = [A, f, U, K, a, b, h](double x) -> double {
			double res = f(x);
			size_t n = A.size();

			for (size_t i = 0; i < n; ++i) {
				double s_i = a + h * i;
				res += A[i] * K(x, s_i) * U[i];
			}
			return res;
		};

		auto diff_squared = [&next_func, &prev_func](double x) {
			return std::pow((next_func(x) - prev_func(x)), 2);
		};

		double L2 = std::sqrt(integrate(diff_squared, a, b, h));
		// double L2 = adaptive_integrate(diff_squared, a, b, 1e-3);

		std::cout << "n: " << n << std::endl;
		std::cout << "L2: "<< L2 << std::endl;

		if (L2 < eps) {
			break;
		}
	}

	return next_func;
}


int main(int argc, char *argv[]) {

	// integration 

	std::cout << std::fixed;
	std::cout << std::setprecision(16);

	double ea = 1e-3;
	double eo = 1e-3;

	// {

	// 	double ea = 1e-6;
	// 	double eo = 1e-6;

	// 	auto f = [](double x) {
	// 		return 1 / x;
	// 	};

	// 	std::cout << adaptive_integrate(f, 1, std::exp(1), ea, eo) << std::endl;
	// }

	// return 0;

	{
		std::cout << "\nSolve equation:\n";
		double a = 0;
		double b = 1;
		auto K = [](double x, double s) {
			return x * std::exp(s) / 2;
		};

		auto f = [](double x) {
			return std::exp(-x);
		};

		auto u = solve_fredholm_equation(K, f, a, b, ea);

		auto y = [](double x) {
			return std::exp(-x) + x;
		};

		auto diff_squared = [&u, &y](double x) {
			return std::pow((u(x) - y(x)), 2);
		};

		std::cout << "\nResult L2: "<< std::sqrt(adaptive_integrate(diff_squared, a, b, ea, eo)) << std::endl;
	}

	// return 0;

	// {
	// 	std::cout << "\n\nSolve equation:\n";
	// 	double a = 0;
	// 	double b = 1;

	// 	auto K = [](double x, double s) {
	// 		return x * s * s;
	// 	};

	// 	auto f = [](double x) {
	// 		return 1;
	// 	};

	// 	auto u = solve_fredholm_equation(K, f, a, b, ea);

	// 	auto y = [](double x) {
	// 		return 1 + (4 * x) / 9;
	// 	};

	// 	int N = 10;
	// 	for (int i = 0; i < N; ++i) {
	// 		double x = double(i) / N;
	// 		std::cout << "u=" << u(x) << ", y=" << y(x) << std::endl;
	// 	}

	// 	auto diff_squared = [&u, &y](double x) {
	// 		return std::pow((u(x) - y(x)), 2);
	// 	};

	// 	std::cout << "\nResult L2: "<< std::sqrt(adaptive_integrate(diff_squared, a, b, ea, eo)) << std::endl;
	// }

	// {
	// 	std::cout << "\n\nSolve equation:\n";
	// 	double a = 0;
	// 	double b = 1;

	// 	auto K = [](double x, double s) {
	// 		return x * s / 2;
	// 	};

	// 	auto f = [](double x) {
	// 		return 5 * x / 6;
	// 	};

	// 	auto u = solve_fredholm_equation(K, f, a, b, ea);

	// 	auto y = [](double x) {
	// 		return x;
	// 	};

	// 	auto diff_squared = [&u, &y](double x) {
	// 		return std::pow((u(x) - y(x)), 2);
	// 	};

	// 	std::cout << "\nResult L2: "<< std::sqrt(adaptive_integrate(diff_squared, a, b, ea, eo)) << std::endl;
	// }
	// // return 0;

	// {
	// 	std::cout << "\n\nSolve equation:\n";
	// 	double a = 0;
	// 	double b = std::atan(1) * 8;

	// 	auto K = [](double x, double s) {
	// 		return std::sin(x) * std::cos(s);
	// 	};

	// 	auto f = [](double x) {
	// 		return std::cos(2 * x);
	// 	};

	// 	auto y = [](double x) {
	// 		return std::cos(2 * x);
	// 	};

	// 	auto u = solve_fredholm_equation(K, f, a, b, ea);

	// 	// int N = 10;
	// 	// for (int i = 0; i < N; ++i) {
	// 	// 	double x =  b * i / N;
	// 	// 	std::cout << "u=" << u(x) << ", y=" << y(x) << std::endl;
	// 	// }

	// 	auto diff_squared = [&u, &y](double x) {
	// 		return std::pow((u(x) - y(x)), 2);
	// 	};

	// 	std::cout << "\nResult L2: "<< std::sqrt(adaptive_integrate(diff_squared, a, b, ea, eo)) << std::endl;
	// }

	// {
	// 	std::cout << "\n\nSolve equation:\n";
	// 	double a = 0;
	// 	double b = 1

	// 	double lambd = 1;

	// 	auto K = [lambd](double x, double s) {
	// 		return lambd * (2 x - s);
	// 	};

	// 	auto f = [](double x) {
	// 		return x / 6;
	// 	};

	// 	auto y = [lambd](double x) {
	// 		return ;
	// 	};

	// 	auto u = solve_fredholm_equation(K, f, a, b);

	// 	auto diff_squared = [&u, &y](double x) {
	// 		return std::pow((u(x) - y(x)), 2);
	// 	};

	// 	std::cout << "\nResult L2: "<< integrate(diff_squared, a, b, h) << std::endl;
	// }

	return 0;
}