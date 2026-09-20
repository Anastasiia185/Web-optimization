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
    double eta;
};

MethodResult goldenSection(double a, double b, double eps) {
    const double tau = (sqrt(5.0) - 1.0) / 2.0;

    double x1 = a + (1.0 - tau) * (b - a);
    double x2 = a + tau * (b - a);
    double f1 = f(x1);
    double f2 = f(x2);
    int Nf = 2;
    int k  = 0;


    cout << setw(4)  << "k"
         << setw(11) << "a_k"
         << setw(11) << "b_k"
         << setw(11) << "x1"
         << setw(11) << "x2"
         << setw(11) << "f(x1)"
         << setw(11) << "f(x2)"
         << setw(11) << "L_k"
         << setw(24) << "Decision"
         << "\n";
    cout << string(105, '-') << "\n";

    while (true) {
        double Lk = b - a;
        bool stop = (Lk <= eps);

        string decision;
        if (stop) {
            decision = "STOP";
        } else if (f1 <= f2) {
            decision = "f1 <= f2 -> b = x2";
        } else {
            decision = "f1 >  f2 -> a = x1";
        }

        cout << setw(4)  << k
             << setw(11) << a
             << setw(11) << b
             << setw(11) << x1
             << setw(11) << x2
             << setw(11) << f1
             << setw(11) << f2
             << setw(11) << Lk
             << "  " << left << setw(22) << decision << "\n";

        if (stop) break;   
        if (f1 <= f2) {
            b  = x2;  x2 = x1;  f2 = f1;
            x1 = a + (1.0 - tau) * (b - a);
            f1 = f(x1);  Nf++;
        } else {
            a  = x1;  x1 = x2;  f1 = f2;
            x2 = a + tau * (b - a);
            f2 = f(x2);  Nf++;
        }
        k++;
    }

    double xStar = (a + b) / 2.0;

    cout << "    x*       = " << xStar << "\n";
    cout << "    f(x*)    = " << f(xStar) << "\n";
    cout << "    k        = " << k << "\n";
    cout << "    Nf       = " << Nf << "  (check: ~ k + 1)\n";
    cout << "    [a, b]   = [" << a << ", " << b << "]\n";
    cout << "    Lk       = " << (b - a) << "\n";

    return {"Golden Section", k, Nf, xStar, f(xStar), b - a, 1.0 / Nf};
}

int main() 
{
    cout << fixed << setprecision(6);

    double a0  = 2.70;
    double b0  = 3.20;
    double eps = 0.001;

    MethodResult rGold = goldenSection(a0, b0, eps);

    MethodResult results[3] = {
        {"Dichotomy",      10, 20, 3.1997, f(3.1997), 0.0006, 1.0/20},
        {"Half Division",  10, 20, 3.1997, f(3.1997), 0.0012, 1.0/20},
        rGold
    };

    cout << "| " << left << setw(20) << "Method"
         << "| " << setw(6)  << "k"
         << "| " << setw(6)  << "Nf"
         << "| " << setw(12) << "x*"
         << "| " << setw(12) << "f(x*)"
         << "| " << setw(12) << "Lk"
         << "| " << setw(12) << "eta = 1/Nf"
         << "|\n";

    cout << "|" << string(21,'-')
         << "|" << string(7,'-')
         << "|" << string(7,'-')
         << "|" << string(13,'-')
         << "|" << string(13,'-')
         << "|" << string(13,'-')
         << "|" << string(13,'-')
         << "|\n";

    for (int i = 0; i < 3; i++) {
        cout << "| " << left << setw(20) << results[i].name
             << "| " << setw(6)  << results[i].iterations
             << "| " << setw(6)  << results[i].Nf
             << "| " << setw(12) << results[i].xStar
             << "| " << setw(12) << results[i].fStar
             << "| " << setw(12) << results[i].length
             << "| " << setw(12) << results[i].eta
             << "|\n";
    }

    return 0;
}