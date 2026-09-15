#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double f(double x) 
{
    return exp(x - 1.0) - 1.0 - x - x*x/2.0 - x*x*x/6.0;
}

int main() 
{
    cout << fixed << setprecision(4);
    double a = 3.20;
    double b = 3.45;
    double eps   = 0.001;
    double sigma = eps / 10.0;

    cout << setw(4)  << "k"
         << setw(10) << "x1"
         << setw(10) << "xm"
         << setw(10) << "x2"
         << setw(12) << "f(x1)"
         << setw(12) << "f(xm)"
         << setw(12) << "f(x2)"
         << setw(10) << "Lk" << "\n";

    int k = 0;
    while ((b - a) > eps) 
    {
        double xm = (a + b) / 2.0;      
        double x1 = xm - sigma;
        double x2 = xm + sigma;
        if (x1 < a) x1 = a;
        if (x2 > b) x2 = b;
        double f1 = f(x1);
        double fm = f(xm);
        double f2 = f(x2);
        double Lk = b - a;   

        cout << setw(4)  << k
             << setw(10) << x1
             << setw(10) << xm
             << setw(10) << x2
             << setw(12) << f1
             << setw(12) << fm
             << setw(12) << f2
             << setw(10) << Lk << "\n";

        if (f1 < f2)
            b = x2;
        else
            a = x1;

        k++;
    }

    double xStar = (a + b) / 2.0;
    cout << "\nResult:\n";
    cout << "x* = " << xStar << "\n";
    cout << "f(x*) = " << f(xStar) << "\n";
    cout << "Iterations: " << k << "\n";
    cout << "Interval: [" << a << ", " << b << "]\n";

    return 0;
}