#include <iostream>
#include <functional>
#include <cmath>
#include <vector>

using vector = std::vector<double>;

namespace hw {

double norm(const vector& v) {
	double to_return;

	for (int i = 0; i < v.size(); ++i) {
		to_return += std::abs(v[i]);
	}
	return to_return;
}

std::ostream& operator<<(std::ostream& os, const vector& v) {
	os << "{ ";
	for (auto& a : v) os << a << " ";
	os << "}";
	return os;
}


std::function<vector(double, vector)> F;

double a = 0;
double b = 1;


double eps = 1e-5;

double X0 = a;
std::vector<double> Y0;


vector operator*(double h, const vector& v) {
	vector u(v);
	for (auto& a : u) a*= h;
	return u;
}

vector operator*(const vector& v, double h) {
	return h * v;
}

vector operator+(const vector& a, const vector& b) {
	vector u(a);
	for (int i = 0; i < a.size(); ++i) {
		u[i] += b[i];
	}
	return u;
}

vector operator/(const vector& v, double h) {
	return v * (1/ h);
}

vector operator-(const vector& a, const vector& b) {
	return a + (-1) * b;
}

void m0(void) {

	int p = 2;

	int counter = 0;

	auto k1 = [&counter](double x, const vector& y, double h) {
		counter++;
		return h * F(x, y);
	};

	auto k2 = [k1, &counter](double x, const vector& y, double h) {
		counter++;
		return h * F(x + h / 2, y + k1(x, y, h) / 2);
	};

	auto g = [k2](double x, const vector& y, double h) {
		return y + k2(x, y, h);
	};


	double h = b - a;

	double x = X0;
	vector y = Y0;


	double mul = 1. / (1 - std::pow(2, -p));


	while (x < b) {

		// std::cout << "x: " << x << std::endl;
		// std::cout << "y: " << y << std::endl;
		// std::cout << "h: " << h << std::endl;
		if (x + h > b) {
			h = b - x;
		}

		vector yh = g(x, y, h);

		vector yh2;

		{
			double h2 = h / 2;
			yh2 = g(x, y, h2);
			yh2 = g(x + h2, yh2, h2);
		}

		double ro = norm(yh - yh2) * mul;

		if (ro > (eps * h) / (b - x)) {
			h /= 2;
		} else {
			y = yh2;
			x = x + h;
			h *= 2;
		}
	}
	std::cout << "method 0." << std::endl;
	std::cout << "y: "<< y << std::endl;
	std::cout << "counter: "<< counter << std::endl;
}

void m1(void) {
	int p = 3;

	int counter = 0;

	auto k1 = [&counter](double x, vector y, double h) {
		counter++;
		return h * F(x, y);
	};

	auto k2 = [k1, &counter](double x, vector y, double h) {
		counter++;
		return h * F(x + h / 2, y + k1(x, y, h) / 2);
	};

	auto g = [k2](double x, vector y, double h) {
		return y + k2(x, y, h);
	};



	auto k2_ = [k1, &counter](double x, vector y, double h) {
		counter++;
		return h * F(x + h / 3, y + k1(x, y, h) / 3);
	};

	auto k3_ = [k2_, &counter](double x, vector y, double h) {
		counter++;
		return h * F(x + 2 * h / 3, y + 2 * k2_(x, y, h) / 3);
	};

	auto g_ = [k1, k3_, &counter](double x, vector y, double h) {
		return y + (k1(x, y, h) + 3 * k3_(x, y, h)) / 4;
	};

	double h = b - a;

	double x = X0;
	vector y = Y0;


	double mul = 1. / (1 - std::pow(2, -p));


	while (x < b) {

		// std::cout << "x: " << x << std::endl;
		// std::cout << "y: " << y << std::endl;
		// std::cout << "h: " << h << std::endl;
		if (x + h > b) {
			h = b - x;
		}

		vector yh = g(x, y, h);

		vector yh_ = g_(x, y, h);

		double ro = norm(yh - yh_) * mul;

		if (ro > (eps * h) / (b - x)) {
			h /= 2;
		} else {
			y = yh_;
			x = x + h;
			h *= 2;
		}
	}
	std::cout << "method 1." << std::endl;
	std::cout << "y: "<< y << std::endl;
	std::cout << "counter: "<< counter << std::endl;
}

void m2(void) {
	int p = 3;
	int counter = 0;

	auto k1_ = [&counter](double x, vector y, double h) {
		counter++;
		return h * F(x, y);
	};

	auto k2_ = [&counter](double x, vector y, double h, vector k1) {
		counter++;
		return h * F(x + h / 2, y + k1 / 2);
	};

	auto k3_ = [&counter](double x, vector y, double h, vector k2) {
		counter++;
		return h * F(x + h / 2, y + k2 / 2);
	};

	auto k4_ = [&counter](double x, vector y, double h, vector k3) {
		counter++;
		return h * F(x + h, y + k3);
	};

	auto g = [](double x, vector y, double h, vector k1, vector k2, vector k3, vector k4) {
		return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
	};

	auto E = [](vector k1, vector k2, vector k3, vector k4) {
		return 2 * (k1 - k2 - k3 + k4) / 3;
	};

	double h = b - a;

	double x = X0;
	vector y = Y0;



	while (x < b) {

		// std::cout << "x: " << x << std::endl;
		// std::cout << "y: " << y << std::endl;
		// std::cout << "h: " << h << std::endl;
		if (x + h > b) {
			h = b - x;
		}

		vector k1 = k1_(x, y, h);
		vector k2 = k2_(x, y, h, k1);
		vector k3 = k3_(x, y, h, k2);
		vector k4 = k4_(x, y, h, k3);

		vector yh = g(x, y, h, k1, k2, k3, k4);

		double ro = norm(E(k1, k2, k3, k4));

		if (ro > (eps * h) / (b - x)) {
			h /= 2;
		} else {
			y = yh;
			x = x + h;
			h *= 2;
		}
	}
	std::cout << "method 2." << std::endl;
	std::cout << "y: "<< y << std::endl;
	std::cout << "counter: "<< counter << std::endl;
}

void test0(void) {
	std::cout <<  std::endl <<  "-----------" << std::endl;
	std::cout << "test 0" << std::endl;
	F = [](double, vector){return vector{1.};};
	a = 0;
	b = 1;
	X0 = a;
	Y0 = {0};

	m0();
	m1();
	m2();

}

void test1(void) {
	std::cout <<  std::endl <<  "-----------" << std::endl;
	std::cout << "test 1" << std::endl;
	F = [](double x, vector y){return vector{x};};
	a = 0;
	b = 1;
	X0 = a;
	Y0 = {0};

	m0();
	m1();
	m2();

}

void test2(void) {
	std::cout <<  std::endl <<  "-----------" << std::endl;
	std::cout << "test 2" << std::endl;
	F = [](double x, vector y){return vector{x * x};};
	a = 0;
	b = 1;
	X0 = a;
	Y0 = {0};

	m0();
	m1();
	m2();
}

void test3(void) {
	std::cout <<  std::endl <<  "-----------" << std::endl;
	std::cout << "test 3" << std::endl;
	F = [](double x, vector y){return vector{y[0]};};

	a = 0;
	b = 1;

	X0 = a;
	Y0 = {1};

	m0();
	m1();
	m2();
}

void test4(void) {
	std::cout <<  std::endl <<  "-----------" << std::endl;
	std::cout << "test 4" << std::endl;

	F = [](double x, vector y) {
		return vector{y[1],y[0]};
	};

	a = 0;
	b = 1;

	X0 = a;
	Y0 = {1, 1};

	m0();
	m1();
	m2();
}
}


int main() {
	using namespace hw;

	std::cout << "eps: " << eps << std::endl;
	// test0();
	test1();
	test2();
	test3();
	test4();

	return 0;
}
