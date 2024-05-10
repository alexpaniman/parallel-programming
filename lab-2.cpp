#include <chrono>
#include <mutex>
#include <numeric>
#include <queue>
#include <thread>
#include <future>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <optional>
#include <sstream>


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

template <typename type>
class locked_stack {
public:
    void push(type value) {
        std::lock_guard<std::mutex> lock(mutex_);
        stack_.push_back(std::move(value));
    }

    std::optional<type> pop() {
        std::lock_guard<std::mutex> lock(mutex_);

        if (stack_.empty())
            return {};

        type value = stack_.back();
        stack_.pop_back();

        return value;
    }

private:
    std::vector<type> stack_;
    std::mutex mutex_;
};


struct integration_task {
    integration_scope scope;
    int number_of_points;
};

class integrator_thread {
public:
    integrator_thread(double precison, locked_stack<integration_task> &tasks, double &result):
        precision_(precison), tasks_(tasks), result_(result) {}


    void operator()() {
        while (std::optional<integration_task> maybe_task = tasks_.pop()) {
            auto task = *maybe_task;

            while (true) {
                double previous_integral = integrate_trapezoid(task.scope, task.number_of_points);
                task.number_of_points *= 2;

                double current_integral = integrate_trapezoid(task.scope, task.number_of_points);
                double difference = std::abs(previous_integral - current_integral);
                if (difference < precision_) {
                    result_ += current_integral;

                    std::stringstream ss;
                    ss << std::this_thread::get_id() << ":\t" << task.scope.a << "\t\t" << task.scope.b << "\n";

                    std::cerr << ss.str();
                    break;
                }

                // Didn't get good enough results yet...
                integration_task new_task = task;
                double middle_point = (task.scope.a + task.scope.b) / 2.;
                new_task.scope.a = middle_point;

                // Shrink current scope in half:
                task.scope.b = middle_point;
                tasks_.push(new_task); // Add half of the task into stack
            }
        }
    }


private:
    double precision_;

    locked_stack<integration_task> &tasks_;
    double &result_;
};




int main() {
    const std::size_t threads_num = 32;
    const double precision = 1e-05;

    integration_scope scope {
        [](double x) { return sin(1./x); },
        1e-5, 10000.
    };

    std::vector<double> results(threads_num, 0);

    locked_stack<integration_task> tasks;
    tasks.push({ scope, 2 });

    {
        std::vector<std::jthread> threads;
        for (std::size_t i = 0; i < threads_num; ++i)
            threads.emplace_back(integrator_thread(precision, tasks, results[i]));
    }

    double integral = std::accumulate(results.begin(), results.end(), 0.); 
    std::cout << integral << "\n";
}
