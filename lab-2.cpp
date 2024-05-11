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
#include <functional>


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


class integrator_pool {
public:
    integrator_pool(integration_scope main_scope, double precision, std::size_t num_threads):
        num_threads_(num_threads), precision_(precision), main_scope_(main_scope),
        results_(num_threads, 0) {}

    double run() {
        tasks_.push({ main_scope_.a, main_scope_.b, 2 });

        {
            std::vector<std::jthread> threads;
            for (std::size_t i = 0; i < num_threads_; ++i)
                threads.emplace_back(std::bind(&integrator_pool::run_single_thread, this, i));
        }

        return std::accumulate(results_.begin(), results_.end(), 0.); 
    }

private:
    std::size_t num_threads_;

    double precision_;
    integration_scope main_scope_;


    struct integration_task {
        double a, b;
        int num_points;
    };

    locked_stack<integration_task> tasks_;
    std::condition_variable finished_;

    std::vector<double> results_;


    void run_single_thread(int thread_index) {
        while (std::optional<integration_task> maybe_task = tasks_.pop()) {
            auto task = *maybe_task;

            while (true) {
                integration_scope scope { main_scope_.f, task.a, task.b };

                double integral_1x_points = integrate_trapezoid(scope, task.num_points);
                double integral_2x_points = integrate_trapezoid(scope, task.num_points *= 2);

                double difference = std::abs(integral_1x_points - integral_2x_points);
                if (difference < precision_) {
                    results_[thread_index] += integral_2x_points;

                    // ----------------- LOGGING -----------------
                    std::stringstream ss;

                    ss << std::this_thread::get_id() << ":\t[";
                    print_precise(ss, scope.a); ss << "\t\t";
                    print_precise(ss, scope.b); ss << "]\n";

                    std::cerr << ss.str();
                    // -------------------------------------------
                    break;
                }

                double middle_point = (scope.a + scope.b) / 2.;
                tasks_.push({ middle_point, task.b, task.num_points });

                // Shrink current scope in half:
                task.b = middle_point;
            }
        }
    }
};




int main() {
    const std::size_t threads_num = 32;
    const double precision = 1e-10;

    integration_scope scope {
        [](double x) { return sin(1./x); },
        1e-10, 1.
    };

    integrator_pool pool(scope, 1e-10, 32);

    std::cout << "integral: " << pool.run() << "\n";
}
