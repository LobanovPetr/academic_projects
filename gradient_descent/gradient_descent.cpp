#include <iostream>
#include <cstring>
#include <string>
#include <initializer_list>
#include <cmath>
#include <functional>
#include <iomanip>
#include <sstream>
#include <limits>

class vector {
public:
	explicit vector(unsigned n);
	explicit vector(std::initializer_list<double>);
	vector(const vector&);
	vector& operator=(const vector&);
	double & operator[](unsigned i);
	double operator[](unsigned i) const;
	unsigned size() const;
	~vector();
private:
	unsigned n_;
	double * value_;
};

vector::vector(unsigned n) : n_(n) {
	if (n_ == 0) {
		throw std::range_error("Can't create vector with size 0");
	}
	value_ = new double[n];
}

vector::vector(std::initializer_list<double> list) {
	n_ = list.size();
	value_ = new double[n_];
	unsigned i = 0;
	for (auto v : list) {
		value_[i++] = v;
	}
}

vector::vector(const vector& obj) {
	n_ = obj.n_;
	value_ = new double [n_];
	memcpy(value_, obj.value_, sizeof(double) * n_);
}

vector& vector::operator=(const vector& vec) {
	if (this == &vec) {
		return *this;
	}
	delete [] value_;
	n_ = vec.n_;
	value_ = new double [n_];
	memcpy(value_, vec.value_, sizeof(double) * n_);
	return *this;
}

double& vector::operator[](unsigned i) {
	return value_[i];
}

double vector::operator[](unsigned i) const {
	return value_[i];
}

unsigned vector::size() const {
	return n_;
}

vector::~vector() {
	delete [] value_;
}

vector operator+(const vector& u, const vector& v) {
	if (u.size() != v.size()) {
		throw std::runtime_error("Different vector sizes: " + std::to_string(u.size()) + " and " + std::to_string(v.size()));
	}

	vector w(u);
	for (unsigned i = 0; i < w.size(); ++i) {
		w[i] += v[i];
	}
	return w;
}

vector operator-(const vector& u, const vector& v) {
	if (u.size() != v.size()) {
		throw std::runtime_error("Different vector sizes: " + std::to_string(u.size()) + " and " + std::to_string(v.size()));
	}

	vector w(u);
	for (unsigned i = 0; i < w.size(); ++i) {
		w[i] -= v[i];
	}
	return w;
}

vector operator*(double a, const vector& v) {
	vector w(v);
	for (unsigned i = 0; i < w.size(); ++i) {
		w[i] *= a;
	}
	return w;
}

class matrix {
public:
	matrix(unsigned n, unsigned m);
	matrix(unsigned n, unsigned m, double fill_value);
	matrix(const matrix& mtx);
	matrix(const vector& vec);
	const matrix& operator=(const matrix& mtx);
	double * operator[](unsigned i);
	const double* operator[](unsigned i) const;
	unsigned rows() const;
	unsigned cols() const;
	operator vector() const;
	~matrix();
private:
	unsigned n_;
	unsigned m_;
	double * values_;
};

matrix::matrix(unsigned n, unsigned m) : n_(n), m_(m) {
	if (n == 0 || m == 0) {
		throw std::range_error("Can't create matrix with zero sizes");
	}
	values_ = new double[n_ * m_];
}

matrix::matrix(unsigned n, unsigned m, double fill_value) : n_(n), m_(m) {
	if (n == 0 || m == 0) {
		throw std::range_error("Can't create matrix with zero size");
	}
	values_ = new double[n_ * m_];
	for (size_t i = 0; i < n_ * m_; ++i) {
		values_[i] = fill_value;
	}
}

matrix::matrix(const matrix& mtx) {
	n_ = mtx.n_;
	m_ = mtx.m_;
	values_ = new double[n_ * m_];
	for (size_t i = 0; i < n_ * m_; ++i) {
		values_[i] = mtx.values_[i];
	}
}

matrix::matrix(const vector& vec) {
	n_ = vec.size();
	m_ = 1;

	values_ = new double[n_ * m_];
	for (size_t i = 0; i < n_ * m_; ++i) {
		values_[i] = vec[i];
	}
}

const matrix& matrix::operator=(const matrix& mtx) {
	if (&mtx == this) {
		return mtx;
	}
	delete [] values_;
	n_ = mtx.n_;
	m_ = mtx.m_;

	values_ = new double [n_ * m_];
	for (size_t i = 0; i < n_ * m_; ++i) {
		values_[i] = mtx.values_[i];
	}
	return *this;
}

double * matrix::operator[](unsigned i) {
	return values_ + m_ * i;
}

const double* matrix::operator[](unsigned i) const {
	return values_ + m_ * i;
}

unsigned matrix::rows() const {
	return n_;
}

unsigned matrix::cols() const {
	return m_;
}

matrix::operator vector() const {
	if (m_ != 1) {
		throw std::runtime_error("Can't create vector from matrix m != 1");
	}
	vector vec(n_);

	for (unsigned i = 0; i < n_; ++i) {
		vec[i] = values_[i];
	}
	return vec;
}

matrix::~matrix() {
	delete [] values_;
}

matrix operator*(const matrix& a, const matrix& b) {
	if (a.cols() != b.rows()) {
		std::stringstream log;
		log << "Incorrect matrices sizes: ("<< a.rows() << "; " << a.cols() << ") and ("
											<< b.rows() << "; " << b.cols() << ")"; 
		throw std::range_error(log.str());
	}
	matrix c(a.rows(), b.cols());

	for (unsigned row = 0; row < a.rows(); ++row) {
		for (unsigned col = 0; col < b.cols(); ++col) {
			double res = 0.;
			for (unsigned i = 0; i < a.cols(); ++i) {
				res += a[row][i] * b[i][col];
			}
			c[row][col] = res;
		}
	}
	return c;
}

std::ostream& operator<<(std::ostream& os, const matrix& mtx) {
	if (mtx.rows() == 0 || mtx.cols() == 0) {
		os << "[]";
		return os;
	}
	os << "[";
	for (unsigned i = 0; i < mtx.rows() - 1; ++i) {
		os << "{";
		for (unsigned j = 0; j < mtx.cols() - 1; ++j) {
			os << mtx[i][j] << "; ";
		}
		os << mtx[i][mtx.cols() - 1] << "}\n";
	}
	os << "{";
	for (unsigned j = 0; j < mtx.cols() - 1; ++j) {
		os << mtx[mtx.rows() - 1][j] << "; ";
	}
	os << mtx[mtx.rows() - 1][mtx.cols() - 1] << "}]";
	return os;
}


matrix inverse(const matrix& a) {
	if (a.rows() != a.cols()) {
		throw std::runtime_error("Can't inverse no square matrix");
	}
	matrix acopy(a);
	matrix ainv(a.rows(), a.rows(), 0.);
	for (unsigned i = 0; i < a.rows(); ++i) {
		ainv[i][i] = 1.;
	}

	auto swap_rows = [](matrix& input, unsigned i, unsigned j){
		if (i == j) {
			return;
		}

		for (unsigned k = 0; k < input.cols(); ++k) {
			std::swap(input[i][k], input[j][k]);
		}
	};

	auto find_non_zero = [](const matrix& input, unsigned i) -> unsigned {
		// double eps = 1e-8;
		for (unsigned j = i; j < input.rows(); ++j) {
			if (std::fabs(input[j][i]) > 0) {
				return j;
			}
		}
		return std::numeric_limits<unsigned>::max();
	};

	auto multiply_row = [](matrix& input, unsigned i, double coeff) {
		for (unsigned j = 0; j < input.cols(); ++j) {
			input[i][j] *= coeff;
		}
	};

	auto add_row = [](matrix& input, unsigned dst, const unsigned src, double coeff) {
		// std::cout << "add "<< src << " to " << dst << std::endl;
		for (unsigned j = 0; j < input.cols(); ++j) {
			input[dst][j] += input[src][j] * coeff;
		}
	};

	for (unsigned col = 0; col < acopy.cols(); ++col) {
		unsigned row_idx = find_non_zero(acopy, col);
		if (row_idx == std::numeric_limits<unsigned>::max()) {
			throw std::runtime_error("Can't find non zero row");
		}

		// std::cout << "row idx: " << row_idx << std::endl;

		swap_rows(acopy, col, row_idx);
		swap_rows(ainv, col, row_idx);

		double val = acopy[col][col];

		// std::cout << "val: " << val << std::endl;

		multiply_row(acopy, col, 1/val);
		multiply_row(ainv, col, 1/val);

		for (unsigned row = 0; row < acopy.rows(); ++row) {
			if (row == col) {
				continue;
			}

			double coeff = -acopy[row][col];
			add_row(acopy, row, col, coeff);
			add_row(ainv, row, col, coeff);
		}

		// std::cout << "acopy:\n"<< acopy << std::endl << std::endl;
	}
	// std::cout << "acopy:\n"<< acopy << std::endl << std::endl;
	return ainv;
}

std::ostream& operator<<(std::ostream& os, const vector& vec) {
	if (vec.size() == 0) {
		os << "{}";
		return os;
	}
	os << "{";
	for (unsigned i = 0; i < vec.size() - 1; ++i) {
		os << vec[i] << "; ";
	}
	os << vec[vec.size() - 1] << "}";
	return os;
}


double derivative(std::function<double(const vector&)> f, unsigned n, const vector& x, double h) {
	vector v(x);
	double result { 0. };
	double xn = x[n];
	// f(x + h) - f(x - h) / 2h

	v[n] = xn + h;

	result += f(v);

	v[n] = xn - h;

	result -= f(v);

	return result / (2 * h);
}

double second_derivative_by_one_variable(std::function<double(const vector&)> f, unsigned n, const vector& x, double h) {
	// f(x - h) - 2 * f(x) + f(x + h) / h^2
	vector v(x);
	double result { 0. };
	double xn = x[n];

	v[n] = xn - h;
	result += f(v);

	v[n] = xn;
	result -= 2 * f(v);

	v[n] = xn + h;
	result += f(v);

	return result / (h * h);
}

double second_derivative_by_different_variables(std::function<double(const vector&)> f, unsigned i, unsigned j, const vector& x, double h) {
	double result {0.};
	vector v(x);

	double xi = x[i];
	double xj = x[j];

	v[i] = xi + h; // (x + h, y)
	v[j] = xj;
	result += f(v);

	v[i] = xi; // (x, y)
	v[j] = xj + h; // (x, y + h)

	result += f(v);

	v[i] = xi;
	v[j] = xj; // (x, y)
	result -= 2 * f(v);

	v[i] = xi + h; // (x + h, y)
	v[j] = xj - h; // (x + h, y - h)
	result -= f(v);

	v[i] = xi - h;
	v[j] = xj + h;
	result -= f(v);

	v[i] = xi;
	v[j] = xj - h;
	result += f(v);

	v[i] = xi - h;
	v[j] = xj;
	result += f(v);

	return result / (2 * h * h);
}

double second_derivative(std::function<double(const vector&)> f, unsigned i, unsigned j, const vector& x, double h) {
	if (i != j) {
		return second_derivative_by_different_variables(f, i, j, x, h);
	}
	return second_derivative_by_one_variable(f, i, x, h);
}

matrix hessian(std::function<double(const vector&)> f, const vector& x, double h) {
	matrix hess(x.size(), x.size());
	for (unsigned i = 0; i < x.size(); ++i) {
		for (unsigned j = 0; j < x.size(); ++j) {
			hess[i][j] = second_derivative(f, i, j, x, h);
		}
	}
	return hess;
}


vector gradient(std::function<double(const vector&)> f, const vector& x, double h) {
	vector grad(x.size());

	for (unsigned i = 0; i < x.size(); ++i) {
		grad[i] = derivative(f, i, x, h);
	}
	return grad;
}

double square_norm(const vector& v) {
	double result {0.};
	for (unsigned i = 0; i < v.size(); ++i) {
		result += v[i] * v[i];
	}
	return std::sqrt(result);
}

// xk+1 = xk - ak * grad(f, xk)
// ak = minarg (f(xk - a * grad(f, xk)))  g(a) = f(xk - a * grad(f, xk))
//
//

double minarg_step_spliting(std::function<double(double)> f, double eps) {
	double b = 2;
	double l = 0.5;

	double a = b;
	double init = f(0.);

	while (f(a) >= init - 1e-6) {
		a = l * a;
	}
	return a;
}

double minarg_gold_ratio(std::function<double(double)> f, double eps) {
	double b = 2;
	double a = 0;

	double r = (3 - std::sqrt(5)) / 2;
	while (b - a > eps) {
		double c = a + r * (b - a);
		double d = b - r * (b - a);

		if (f(c) < f(d)) {
			b = d;
		} else {
			a = c;
		}
	}
	return (a + b) / 2;
}



vector gradient_descent(std::function<double(const vector&)> f, const vector& init_point, std::function<double(std::function<double(double)>, double)> minarg, const double eps) {
	vector prev_x(init_point.size());
	double prev_f;

	vector next_x = init_point;
	double next_f = f(init_point);

	double x_stop;
	double f_stop;
	double grad_stop;

	int step = 0;
	std::cout << "START GRADIENT DESCENT:\n";
	std::cout << "INPUT:\n";
	std::cout << "\tx_0 = " << next_x << std::endl;
	std::cout << "\tf(x_0) = " << next_f << std::endl;
	std::cout << "\n\n";
	do {
		prev_x = next_x;
		prev_f = next_f;

		vector grad = gradient(f, prev_x, eps);
		auto one_arg_func = [f, &prev_x, &grad](double a){return f(prev_x - a * grad);};
		double a = minarg(one_arg_func, eps);

		next_x = prev_x - a * grad;
		next_f = f(next_x);
		vector next_grad = gradient(f, next_x, eps);

		std::cout << "STEP " << step + 1 << std::endl;
		std::cout << "\tx_" << step << " = " << prev_x << std::endl;
		std::cout << "\th_" << step << " = f'(x_" << step << ") = " << grad << std::endl;
		std::cout << "\ta_" << step << " = " << a << std::endl;

		std::cout << "\n";
		std::cout << "\tx_" << step + 1 << " = " << "x_" << step << " + a_" << step << " * h_" << step << std::endl;
		std::cout << "\tx_" << step + 1 << " = " << next_x << std::endl;
		std::cout << "\tf(x_" << step + 1 << ") = " << next_f << std::endl;
		std::cout << "\tf'(x_" << step + 1 << ") = " << next_grad << std::endl;

		x_stop = square_norm(next_x - prev_x);
		f_stop = std::fabs(next_f - prev_f);
		grad_stop = square_norm(next_grad);

		std::cout << "\n";
		std::cout << "\t|x_" << step + 1 << " - x_" << step << "| = " << x_stop << std::endl;
		std::cout << "\t|f(x_" << step + 1 << ") - f(x_" << step << ")| = " << f_stop << std::endl;
		std::cout << "\t|f'(x_" << step + 1 << ")| = " << grad_stop << std::endl;
		std::cout << "\n\n\n";
		step++;
	} while (x_stop > eps || f_stop > eps || grad_stop > eps);
	return next_x;
}

vector gradient_descent_newton(std::function<double(const vector&)> f, const vector& init_point, std::function<double(std::function<double(double)>, double)> minarg, const double eps) {
	vector prev_x(init_point.size());
	double prev_f;

	vector next_x = init_point;
	double next_f = f(init_point);

	double x_stop;
	double f_stop;
	double grad_stop;

	int step = 0;
	std::cout << "START NEWTON GRADIENT DESCENT:\n";
	std::cout << "INPUT:\n";
	std::cout << "\tx_0 = " << next_x << std::endl;
	std::cout << "\tf(x_0) = " << next_f << std::endl;
	std::cout << "\n\n";
	do {
		prev_x = next_x;
		prev_f = next_f;

		vector grad = gradient(f, prev_x, eps);
		matrix hess = hessian(f, prev_x, eps);
		matrix inv_hess(1,1);
		try {
			inv_hess = inverse(hess);
		} catch (std::exception& e) {
			std::cout << "ZERO HESSIAN.\n";
			return next_x;
		}
		vector h = inv_hess * (-1 * grad);
		auto one_arg_func = [f, &prev_x, &h](double a){return f(prev_x + a * h);};
		double a = minarg(one_arg_func, eps);

		next_x = prev_x + a * h;
		next_f = f(next_x);
		vector next_grad = gradient(f, next_x, eps);

		std::cout << "STEP " << step + 1 << std::endl;
		std::cout << "\tx_" << step << " = " << prev_x << std::endl;
		std::cout << "\tf'(x_" << step << ") = " << grad << std::endl;
		std::cout << "\tf''(x_" << step << ") =\n" << hess << std::endl << std::endl;
		std::cout << "\tf''(x_" << step << ")^(-1) = \n" << inv_hess << std::endl << std::endl;
		std::cout << "\th_" << step << " = " << h << std::endl;
		std::cout << "\ta_" << step << " = " << a << std::endl;

		std::cout << "\n";
		std::cout << "\tx_" << step + 1 << " = " << "x_" << step << " + a_" << step << " * h_" << step << std::endl;
		std::cout << "\tx_" << step + 1 << " = " << next_x << std::endl;
		std::cout << "\tf(x_" << step + 1 << ") = " << next_f << std::endl;
		std::cout << "\tf'(x_" << step + 1 << ") = " << next_grad << std::endl;

		x_stop = square_norm(next_x - prev_x);
		f_stop = std::fabs(next_f - prev_f);
		grad_stop = square_norm(next_grad);

		std::cout << "\n";
		std::cout << "\t|x_" << step + 1 << " - x_" << step << "| = " << x_stop << std::endl;
		std::cout << "\t|f(x_" << step + 1 << ") - f(x_" << step << ")| = " << f_stop << std::endl;
		std::cout << "\t|f'(x_" << step + 1 << ")| = " << grad_stop << std::endl;
		std::cout << "\n\n\n";
		step++;
	} while (x_stop > eps || f_stop > eps || grad_stop > eps);
	return next_x;
}


// two dim functions:

// (x0 - 1)^2 + (x1 + 1)^2
auto f0 = [](const vector& x){ 
	return (x[0] - 1) * (x[0] - 1) + (x[1] + 1) * (x[1] + 1);
};

// x0^2 + x1^4
auto f1 = [](const vector& x) {
	return x[0] * x[0] + x[1] * x[1] * x[1] * x[1];
};

// x1^2 - x2^2
auto f2 = [](const vector& x) {
	return x[0] * x[0] - x[1] * x[1];
};

auto g = [](const vector& x) {
	return x[0] * x[1];
};

auto f3 = [](const vector& x) {
	return x[0] * x[0] + x[1] * x[1];
};


// exp()
auto f4 = [](const vector& x) {
	return std::exp(f3(x));
};


auto f5 = [](const vector& x) {
	return -std::cos(f3(x));
};




// three dimension fucntions:

auto g0 = [](const vector& x) {
	return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
};


auto g1 = [](const vector& x) {
	return (x[0] - 1) * (x[0] - 1) + 10 * (x[1] - 2) * (x[1] - 2) + 0.1 * (x[2] -3) * (x[2] -3);
};

auto g2 = [](const vector& x) {
	return (x[0] - 1) * (x[0] - 1) * (x[0] - 1) * (x[0] - 1) + 10 * (x[1] - 2) * (x[1] - 2) + 0.1 * (x[2] -3) * (x[2] -3);
};

auto g3 = [](const vector& x) {
	return std::exp(g1(x));
};

int main(int args, char* argv[]) {

	std::cout << std::fixed;
	std::cout << std::setprecision(8);

	
	double eps = 1e-3;

	auto f = g3;
	std::cout << "STEP SPLITING ALGO:\n";
	vector p = gradient_descent(f, vector({0.9, 1.9, 2.9}), minarg_step_spliting, std::sqrt(eps));

	std::cout << "GOLD RATIO ALGO:\n";
	gradient_descent_newton(f, p, minarg_gold_ratio, eps);

	return 0;
}
