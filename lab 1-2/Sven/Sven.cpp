#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

double f(double x) 
{
    return exp(x - 1.0) - 1.0 - x - x*x/2.0 - x*x*x/6.0;
}

int main()
{
    double x0 = 3.20;
    double delta = 0.5;

    cout << fixed << setprecision(4);
    cout << setw(6)  << "k"
         << setw(12) << "Delta"
         << setw(12) << "xk"
         << setw(14) << "f(xk)"
         << setw(22) << "f(xk)<f(xk-1)?"
         << setw(24) << "[a0, b0]" << "\n";

    double f0 = f(x0);
    double f1 = f(x0 + delta);

    if (f0 <= f1) delta = -delta;

    double xk  = x0;
    double xk1 = x0 + delta;
    double fk  = f0;
    double fk1 = f(xk1);

    int k = 0;
    cout << setw(6)  << k
         << setw(12) << "-"
         << setw(12) << x0
         << setw(14) << f0
         << setw(22) << "-"
         << setw(24) << "-" << "\n";

    while (fk1 < fk) 
    {
        k++;
        cout << setw(6)  << k
             << setw(12) << delta
             << setw(12) << xk1
             << setw(14) << fk1
             << setw(22) << "+"
             << setw(24) << "-" << "\n";

        delta *= 2.0;
        xk  = xk1;
        fk  = fk1;
        xk1 = xk + delta;
        fk1 = f(xk1);
    }

    k++;
    double a = (xk < xk1) ? xk : xk1;
    double b = (xk < xk1) ? xk1 : xk;

    cout << setw(6)  << k
         << setw(12) << delta
         << setw(12) << xk1
         << setw(14) << fk1
         << setw(22) << "-"
         << " [" << a << ", " << b << "]"
         << "\n";

    cout << "\nFinal interval: [" << a << ", " << b << "]\n";
}