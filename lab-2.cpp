#include <chrono>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <queue>
#include <thread>
#include <future>
#include <cmath>
#include <vector>
#include <iostream>
#include <optional>
#include <sstream>


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



static void print_precise(std::ostream &os, double number) {
    os << std::setw(20) << std::fixed << std::setprecision(10) << number;
}


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

                    ss << std::this_thread::get_id() << ":\t[";
                    print_precise(ss, task.scope.a); ss << "\t\t";
                    print_precise(ss, task.scope.b); ss << "]\n";

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
    const double precision = 1e-10;

    integration_scope scope {
        [](double x) { return sin(1./x); },
        1e-5, 1.
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
    std::cout << "integral: " << integral << "\n";
}
