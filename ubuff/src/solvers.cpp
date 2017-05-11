#include "cptstd.h"
#include "matrix.h"
using namespace cpt;

int L = 50;                     // number of interior points in x and y
Matrix<double,2> V(L+2, L+2),   // potential to be found
        rho(L+2, L+2),              // given charge density
        V_new(L+2, L+2);            // new potential after each step

double h;                       // lattice spacing
int steps;                      // number of iteration steps
double accuracy;                // desired accuracy in solution
double omega;                   // overrelaxation parameter

void initialize()
{
    int N = L + 2;
    V = rho = V_new = Matrix<double,2>(N, N);

    h = 1 / double(L + 1);      // assume physical size in x and y = 1
    double q = 10;              // point charge
    int i = N / 2;              // center of lattice
    rho[i][i] = q / (h * h);    // charge density
    steps = 0;
}

void Jacobi() {
    // Jacobi algorithm for a single iterative step
    for (int i = 1; i <= L; i++)
        for (int j = 1; j <= L; j++)
            V_new[i][j] = 0.25 * (V[i - 1][j] + V[i + 1][j] +
                                  V[i][j - 1] + V[i][j + 1] +
                                  h * h * rho[i][j]);
}
double relative_error()
{
    double error = 0;           // average relative error per lattice point
    int n = 0;                  // number of non-zero differences

    for (int i = 1; i <= L; i++)
        for (int j = 1; j <= L; j++) {
            if (V_new[i][j] != 0)
                if (V_new[i][j] != V[i][j]) {
                    error += abs(1 - V[i][j] / V_new[i][j]);
                    ++n;
                }
        }
    if (n != 0)
        error /= n;
    return error;
}

void Gauss_Seidel()
{
    // copy V to V_new
    V_new = V;

    // Gauss-Seidel update in place
    for (int i = 1; i <= L; i++)
        for (int j = 1; j <= L; j++)
            V_new[i][j] = 0.25 * (V_new[i - 1][j] + V_new[i + 1][j] +
                                  V_new[i][j - 1] + V_new[i][j + 1] +
                                  h * h * rho[i][j]);
}

void successive_over_relaxation()   // using red-black checkerboard updating
{
    // update even sites first
    for (int i = 1; i <= L; i++)
        for (int j = 1; j <= L; j++)
            if ((i + j) % 2 == 0)
                V_new[i][j] = (1 - omega) * V[i][j] + omega / 4 *
                                                      (V[i - 1][j] + V[i + 1][j] +
                                                       V[i][j - 1] + V[i][j + 1] +
                                                       h * h * rho[i][j]);

    // update odd sites using updated even sites
    for (int i = 1; i <= L; i++)
        for (int j = 1; j <= L; j++)
            if ((i + j) % 2 != 0)
                V_new[i][j] = (1 - omega) * V[i][j] + omega / 4 *
                                                      (V_new[i - 1][j] + V_new[i + 1][j] +
                                                       V_new[i][j - 1] + V_new[i][j + 1] +
                                                       h * h * rho[i][j]);
}

void iterate(void (*method)())
{
    clock_t t0 = clock();

    while (true) {
        method();
        ++steps;
        double error = relative_error();
        if (error < accuracy)
            break;
        swap(V, V_new);         // use <algorithm> std::swap
    }
    cout << " Number of steps = " << steps << endl;

    clock_t t1 = clock();
    cout << " CPU time = " << double(t1 - t0) / CLOCKS_PER_SEC
         << " sec" << endl;
}
