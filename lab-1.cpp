#include <fstream>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

const int ROOT = 0; // ранг главного процесса

const double t_max = 50;
const double x_max = 60;
const double t_step = 0.1;
const double x_step = 0.1;

double func(double t, double x);
double fi(double x);
double ksi(double t);

double func(double t, double x){
    return x*t;
}

double fi(double x){
    return x*x*x / 12;
}

double ksi(double t){
    return t*t*t / 12;
}

void fill_layer(int knot_t, int knot_x, std::vector<std::vector<double>> &u) {
    if (knot_x != x_max / x_step - 1) {
        // центральная явная трехточечная схема
        u[knot_t + 1][knot_x] = func(knot_t, knot_x) * t_step + 
                      1 / 2. * (u[knot_t][knot_x + 1] + 
                      u[knot_t][knot_x - 1]) -  
                      t_step / (2. * x_step) * (u[knot_t][knot_x + 1] - u[knot_t][knot_x - 1]);
    } 
    
    else {
        // явный левый уголок
        u[knot_t + 1][knot_x] = func(knot_t, knot_x) * t_step + 
                      u[knot_t][knot_x] -
                      t_step / x_step * (u[knot_t][knot_x] - u[knot_t][knot_x - 1]);
    }
}

// Печать количества процессов и времени выполнения
void output_results(int num_knots_t, int num_knots_x, std::vector<std::vector<double>> &u) {
    std::ofstream results("results.csv");
    results << "x,t,u\n";

    for (int knot_t = 0; knot_t < num_knots_t; knot_t++) {
        double t = knot_t * t_step;

        for (int knot_x = 0; knot_x < num_knots_x; knot_x++) {
            double x = knot_x * x_step;
            results << x << "," << t << "," << u[knot_t][knot_x] << "\n";
        }
    }
}

std::vector<std::vector<double>> create_u(int &num_knots_x, int &num_knots_t) {
    std::vector<std::vector<double>> u;
    u.reserve(num_knots_t);

    for (int knot_t = 0; knot_t < num_knots_t; knot_t++) {
            u[knot_t].reserve(num_knots_x);
    }

    return u;
}

void init_u(int &num_knots_x, int &num_knots_t, std::vector<std::vector<double>> &u) {
    for (int m = 0; m < num_knots_x; ++m)
        u[0][m] = fi(m); // инициализация при t = 0

    for (int k = 0; k < num_knots_t; ++k)
        u[k][0] = ksi(k);
}

void common_process_computation(int &rank, int &num_knots_t, std::vector<std::vector<double>> &u, int &x_0, int &x_1, int &knot_t) {
    if (knot_t > 0) {
        // принимаем соседние точки к крайним
        MPI_Recv(&u[knot_t][x_0 - 1], 1, MPI_DOUBLE, rank - 1, knot_t,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&u[knot_t][x_1], 1, MPI_DOUBLE, rank + 1, knot_t, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    fill_layer(knot_t, x_0, u);
    fill_layer(knot_t, x_1 - 1, u);

    if (knot_t < num_knots_t - 2) {
        // отправляем наши вычисления соседям
        MPI_Send(&u[knot_t + 1][x_0], 1, MPI_DOUBLE, rank - 1, knot_t + 1, MPI_COMM_WORLD);
        MPI_Send(&u[knot_t + 1][x_1 - 1], 1, MPI_DOUBLE, rank + 1, knot_t + 1, MPI_COMM_WORLD);
    }
}

void root_process_computation(int &rank, int &num_knots_t, std::vector<std::vector<double>> &u, int &x_1, int &knot_t) {
    if (knot_t > 0)
        MPI_Recv(&u[knot_t][x_1], 1, MPI_DOUBLE, rank + 1, knot_t, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);

    fill_layer(knot_t, x_1 - 1, u);

  if (knot_t < num_knots_t - 2)
    MPI_Send(&u[knot_t + 1][x_1 - 1], 1, MPI_DOUBLE, rank + 1, knot_t + 1, MPI_COMM_WORLD);
}

void last_process_computation(int &rank, int &num_knots_t, std::vector<std::vector<double>> &u, int &x_0, int &knot_t) {
    if (knot_t > 0)
        MPI_Recv(&u[knot_t][x_0 - 1], 1, MPI_DOUBLE, rank - 1, knot_t, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    fill_layer(knot_t, x_0, u);

    if (knot_t < num_knots_t - 2)
        MPI_Send(&u[knot_t + 1][x_0], 1, MPI_DOUBLE, rank - 1, knot_t + 1, MPI_COMM_WORLD);
}

void parallel_solution(int rank, int num_knots_x, int num_knots_t, int size, std::vector<std::vector<double>> &u) {
    // количество вычисляемых точек для всех процессов, кроме главного
    int part = num_knots_x / size;

    // количество сдвига (доп точки для главного процесса)
    int shift = num_knots_x % size;

    // в зависимости от процесса выбор
    // количества вычисляемых точек
    int num_knots = (rank == ROOT) ? (part + shift) : part;

    // первая точка вычислений
    int x_0 = (rank == ROOT) ? (0) : (part * rank + shift);

    // последняя точка вычислений (не включительно)
    int x_1 = x_0 + num_knots;

    for (int knot_t = 0; knot_t < num_knots_t - 1; knot_t++) { // время
        for (int knot_x = x_0 + 1; knot_x < x_1 - 1; knot_x++) { // координата
            // заполняем все внутренние точки, неоходимые
            // данные для вычислений которых процесс уже знает
            fill_layer(knot_t, knot_x, u);
        }

        if (rank != ROOT && rank != size - 1)
            common_process_computation(rank, num_knots_t, u, x_0, x_1, knot_t);

        // главный процесс выполняет вычисления, начиная от первой
        // точки, поэтому ему не нужны сведения о -1 точке, которая
        // не входит в сетку
        else if (rank == ROOT)
            root_process_computation(rank, num_knots_t, u, x_1, knot_t);

        // последнему процессу не нужны сведения о M+1 точке, которая не
        // входит в сетку, так как для точки M вычисленияя выполняются уголком
        else 
            last_process_computation(rank, num_knots_t, u, x_0, knot_t);

        // следующий узел можно начинать считать
        // только при вычислении предыдущего
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    // принимаем окончательные данные главным процессом
    if (rank == ROOT) { 
        for (int r = 1; r < size; r++) {
            for (int k = 1; k < num_knots_t; k++) {
                int first = part * r + shift;
                int last = first + part;
                MPI_Recv(&u[k][first], part, MPI_DOUBLE, r, num_knots_t,
                                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    
    // остальные процессы отправляют вычисленные данные
    else { 
        for (int k = 1; k < num_knots_t; k++) {
            MPI_Send(&u[k][x_0], num_knots, MPI_DOUBLE, ROOT, num_knots_t,
                                MPI_COMM_WORLD);
        }
    }

    if (rank == ROOT)
        output_results(num_knots_t, num_knots_x, u);
}

void serial_solution(int num_knots_t, int num_knots_x, std::vector<std::vector<double>> &u) {
    for (int knot_t = 0; knot_t < num_knots_t - 1; knot_t++)
        for (int knot_x = 1; knot_x < num_knots_x; knot_x++)
            fill_layer(knot_t, knot_x, u); // последовательно заполняем слои

    output_results(num_knots_t, num_knots_x, u);
}

int main(int argc, char **argv) {
    int num_knots_x = x_max / x_step; // количество узлов по x
    int num_knots_t = t_max / t_step; // количество узлов по t

    double start_time,
            stop_time; // время начала и конца выполнения соответственно
    int rank, size; // ранг процесса и количество процессов

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::vector<double>> u = create_u(num_knots_x, num_knots_t);

    if (rank == ROOT) 
        start_time = MPI_Wtime(); // начало отсчета времени
    
    init_u(num_knots_x, num_knots_t, u); // инициализация при x = 0

    // эта часть выполняется в случае параллельной программы
    if (size > 1)
        parallel_solution(rank, num_knots_x, num_knots_t, size, u);

    // выполняется в случае последовательной программы
    else
        serial_solution(num_knots_t, num_knots_x, u);

    MPI_Finalize();

    return 0;
}
