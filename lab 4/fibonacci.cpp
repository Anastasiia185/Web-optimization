#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

double f(double x) {
    return exp(x - 1.0) - 1.0 - x - x*x/2.0 - x*x*x/6.0;
}

struct MethodResult {
    string name;
    int    iterations;
    int    Nf;
    double xStar;
    double fStar;
    double length;
};

int buildFibonacci(double L0, double eps, long long F[]) {
    double R = L0 / eps;
    F[0] = 0;  F[1] = 1;  F[2] = 1;
    int n = 2;
    while ((double)F[n] < R) {
        n++;
        F[n] = F[n-1] + F[n-2];
    }
    return n;
}

MethodResult fibonacci(double a, double b, double eps) {
    double L0 = b - a;
    long long F[50];
    int N = buildFibonacci(L0, eps, F);

    cout << setw(4)  << "k"
         << setw(4)  << "m"
         << setw(12) << "a_k"
         << setw(12) << "b_k"
         << setw(12) << "x1"
         << setw(12) << "x2"
         << setw(12) << "f(x1)"
         << setw(12) << "f(x2)"
         << setw(12) << "L_k"
         << setw(15) << "Decision"
         << "\n";
    cout << string(115, '-') << "\n";

    int Nf = 0;
    int k  = 0;

    int    m  = N;
    double x1 = a + (double)F[m-2] / F[m] * (b - a);
    double x2 = a + (double)F[m-1] / F[m] * (b - a);
    double f1 = f(x1);
    double f2 = f(x2);
    Nf = 2;

    cout << setw(4)  << k
         << setw(4)  << m
         << setw(12) << a
         << setw(12) << b
         << setw(12) << x1
         << setw(12) << x2
         << setw(12) << f1
         << setw(12) << f2
         << setw(12) << (b - a)
         << setw(15) << "-"
         << "\n";

    for (k = 1; k <= N - 2; k++) {
        m = N - k;

        double Lk = b - a;   

        string decision;
        if (f1 <= f2) {
            decision = "f1<=f2 -> b=x2";
            b  = x2;
            x2 = x1;  f2 = f1;         
            // новое x1
            x1 = a + (double)F[m-2] / F[m] * (b - a);
            f1 = f(x1);
            Nf++;
        } else {
            decision = "f1>f2 -> a=x1";
            a  = x1;
            x1 = x2;  f1 = f2;         
            // новое x2
            x2 = a + (double)F[m-1] / F[m] * (b - a);
            f2 = f(x2);
            Nf++;
        }

        cout << setw(4)  << k
             << setw(4)  << m
             << setw(12) << a
             << setw(12) << b
             << setw(12) << x1
             << setw(12) << x2
             << setw(12) << f1
             << setw(12) << f2
             << setw(12) << Lk
             << setw(15) << decision
             << "\n";
    }

    if (f1 <= f2) b = x2;
    else          a = x1;

    double xStar = (a + b) / 2.0;
    double guarantee = L0 / F[N];

    cout << "    x*        = " << xStar << "\n";
    cout << "    f(x*)     = " << f(xStar) << "\n";
    cout << "    Iterations= " << k - 1 << "\n";
    cout << "    Nf        = " << Nf << "\n";
    cout << "    Lk        = " << (b - a) << "\n";
    cout << "    L0/F[N]   = " << guarantee << "\n";
    cout << "    Guarantee Lk <= L0/F[N]? "
         << ((b - a) <= guarantee ? "+" : "-") << "\n";

    return {"Fibonacci", k - 1, Nf, xStar, f(xStar), b - a};
}

int main() 
{
    cout << fixed << setprecision(6);

    double a0  = 2.70;
    double b0  = 3.20;
    double eps = 0.001;

    MethodResult rFib = fibonacci(a0, b0, eps);
    MethodResult rGold = {"Golden Section", 13, 15, 3.199590, f(3.199590), 0.000820};

    MethodResult results[2] = { rFib, rGold };
    cout << "| " << left << setw(18) << "Method"
         << "| " << setw(6)  << "k"
         << "| " << setw(6)  << "Nf"
         << "| " << setw(12) << "x*"
         << "| " << setw(12) << "f(x*)"
         << "| " << setw(12) << "Lk"
         << "| " << setw(12) << "|x*-x*_t|"
         << "|\n";

    cout << "|" << string(19,'-') << "|" << string(7,'-')
         << "|" << string(7,'-')  << "|" << string(13,'-')
         << "|" << string(13,'-') << "|" << string(13,'-')
         << "|" << string(13,'-') << "|\n";

    double xExact = 3.1995;
    for (int i = 0; i < 2; i++) {
        double err = fabs(results[i].xStar - xExact);
        cout << "| " << left << setw(18) << results[i].name
             << "| " << setw(6)  << results[i].iterations
             << "| " << setw(6)  << results[i].Nf
             << "| " << setw(12) << results[i].xStar
             << "| " << setw(12) << results[i].fStar
             << "| " << setw(12) << results[i].length
             << "| " << setw(12) << err
             << "|\n";
    }

    return 0;
}