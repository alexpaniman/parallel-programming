#include <numeric>
#include <thread>
#include <future>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <iostream>


const size_t MAX_ITERATIONS = 32;

typedef double (*integration_func_type)(double x);

struct integration_scope {
    integration_func_type f;
    double a, b;
};

typedef double (*integration_method_func_type)(integration_scope scope, int number_of_points);

static double calculate_delta(double a, double b, int number_of_points) {
    return (b - a) / (number_of_points - 1);
}

double integrate_trapezoid(integration_scope scope, int number_of_points) {
    double sum = 0;

    double delta = calculate_delta(scope.a, scope.b, number_of_points);
    for (int i = 0; i < number_of_points - 1; ++ i) {
        sum += delta * (scope.f(scope.a + delta * i) + scope.f(scope.a + delta * (i + 1))) / 2;
    }

    return sum;
}

double integrate_simpson(integration_scope scope, int number_of_points) {
    double sum = 0;

    double delta = calculate_delta(scope.a, scope.b, number_of_points);
    for (int i = 0; i < number_of_points - 1; ++ i) {
        sum += delta * (scope.f(scope.a + delta * i)
                        + 4 * scope.f(scope.a + delta * (i + 0.5))
                        + scope.f(scope.a + delta * (i + 1))) / 6;
    }

    return sum;
}

double integrate_to_precision(integration_scope scope, integration_method_func_type method, double eps, const char* name = nullptr) {
    int number_of_points = 2;

    double previous_integral = method(scope, number_of_points);
    number_of_points *= 2;

    double current_integral = 0;

    double min_difference = INFINITY;
    for (size_t i = 0; i < MAX_ITERATIONS; ++ i) {
        current_integral = method(scope, number_of_points);

        double difference = std::abs(previous_integral - current_integral);
        if (name) {
            printf("%s fabs(%lf - %lf)\t=\t%lf\n", name, current_integral, previous_integral, difference);
        }

        min_difference = std::min(min_difference, difference);

        if (!std::isfinite(current_integral)) {
            fprintf(stderr, "I tried but could reach eps %lf since current method iteration gave %lf\n", eps, current_integral);
            return previous_integral;
        }

        if (difference < eps)
            return current_integral;

        previous_integral = current_integral;
        number_of_points *= 2;
    }

    fprintf(stderr, "Giving up, could not reached desired goal of %lf eps, got %lf eps at best\n", eps, min_difference);
    return NAN;
}



int main() {
    const std::size_t threads_num = 15;
    const double precision = 1e-9;

    integration_scope scope {
        [](double x) { return sin(1/x); },
        .00001, 1.
    };

    std::vector<double> results(threads_num);
    {
        std::vector<std::jthread> threads;

        for (std::size_t i = 0; i < threads_num; ++i) {
            double subscope_length = (scope.b - scope.a) / threads_num;
            integration_scope subscope { scope.f,
                scope.a + subscope_length * i,
                scope.a + subscope_length * (i + 1)
            };
            std::cout << i << "\t[" << subscope.a << ",\t" << subscope.b << "]\n";

            threads.emplace_back([=, &result = results[i]] {
                result = integrate_to_precision(subscope, integrate_trapezoid, precision);
            });
        }
    }

    double integral = std::accumulate(results.begin(), results.end(), 0.); 
    std::cout << integral << "\n";
}
