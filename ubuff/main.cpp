#include "src/cptstd.h"
#include "src/matrix.h"
#include "src/solvers.cpp"

int main() {

    cout << " Iterative solution of Poisson's equation\n"
         << " ----------------------------------------\n";
    cout << " Enter number of interior points in x or y: ";
    cin >> L;

    initialize();

    cout << " Enter desired accuracy in solution: ";
    cin >> accuracy;
    cout << " Enter 1 for Jacobi, 2 for Gauss Seidel, 3 for SOR: ";
    int choice;
    cin >> choice;

    switch (choice) {
    case 1:
        iterate(Jacobi);
        break;
    case 2:
        iterate(Gauss_Seidel);
        break;
    case 3:
        cout << " Enter overrelaxation parameter omega: ";
        cin >> omega;
        iterate(successive_over_relaxation);
        break;
    default:
        cout << " Jacobi: " << endl;
        iterate(Jacobi);
        cout << " Gauss-Seidel: " << endl;
        initialize();
        iterate(Gauss_Seidel);
        omega = 2 / (1 + 4 * atan(1.0) / double(L));
        cout << " Successive Over Relaxation with theoretical optimum omega = "
             << omega << endl;
        initialize();
        iterate(successive_over_relaxation);
        break;
    }

    // write potential to file
    cout << " Potential in file poisson.data" << endl;
    ofstream date_file("poisson.data");
    for (int i = 0; i < L + 2; i++) {
        double x = i * h;
        for (int j = 0; j < L + 2; j++) {
            double y = j * h;
            date_file << x << '\t' << y << '\t' << V[i][j] << '\n';
        }
        date_file << '\n';
    }
    date_file.close();
}