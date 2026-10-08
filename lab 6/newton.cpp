#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

double f(double x) {
    return exp(x - 1.0) - 1.0 - x - x*x/2.0 - x*x*x/6.0;
}

double df(double x) {
    return exp(x - 1.0) - 1.0 - x - x*x/2.0;
}

double d2f(double x) {
    return exp(x - 1.0) - 1.0 - x;
}

struct MethodResult {
    string name;
    int    iterations;
    int    Nf;
    double xStar;
    double fStar;
    double error;
};

MethodResult newtonClassic(double x0, double eps) {
    cout << "\nPART A:\n";
    cout << "x0 = " << x0 << ", eps = " << eps << "\n\n";

    cout << setw(4)  << "k"
         << setw(12) << "x_k"
         << setw(14) << "f(x_k)"
         << setw(14) << "f'(x_k)"
         << setw(14) << "|f'(x_k)|"
         << setw(14) << "f''(x_k)"
         << "\n";
    cout << string(72, '-') << "\n";

    double x = x0;
    int Nf = 0;
    int k  = 0;

    cout << setw(4)  << k
         << setw(12) << x
         << setw(14) << f(x)
         << setw(14) << df(x)
         << setw(14) << fabs(df(x))
         << setw(14) << d2f(x)
         << "\n";
    Nf += 3;

    while (fabs(df(x)) > eps) {
        double d1 = df(x);
        double d2 = d2f(x);

        if (d2 < 1e-12) {
            cout << "d2 <= 0 -- stop (would go to maximum)\n";
            break;
        }

        x = x - d1 / d2;
        k++;

        cout << setw(4)  << k
             << setw(12) << x
             << setw(14) << f(x)
             << setw(14) << df(x)
             << setw(14) << fabs(df(x))
             << setw(14) << d2f(x)
             << "\n";
        Nf += 3;
    }

    double xStar = x;
    double fStar = f(x);
    double xExact = 3.1995;
    double err = fabs(xStar - xExact);
    cout << "    x* = " << xStar << "\n";
    cout << "    f(x*)= " << fStar << "\n";
    cout << "    Iterations= " << k << "\n";
    cout << "    Nf= " << Nf << "\n";
    cout << "    |x*-x*_t| = " << err << "\n";

    return {"Newton (classic)", k, Nf, xStar, fStar, err};
}

MethodResult newtonRaphson(double x0, double eps, double alpha) {
    cout << "\nPART B:\n";
    cout << "x0 = " << x0 << ", eps = " << eps
         << ", alpha = " << alpha << "\n\n";

    cout << setw(4)  << "k"
         << setw(12) << "x_k"
         << setw(12) << "f(x_k)"
         << setw(12) << "f'(x_k)"
         << setw(12) << "f''(x_k)"
         << setw(8)  << "t"
         << setw(12) << "x_bar"
         << setw(12) << "f(x_bar)"
         << setw(12) << "f'(x_bar)"
         << setw(12) << "(3.10)"
         << "\n";
    cout << string(108, '-') << "\n";

    double x = x0;
    int Nf = 0;
    int k  = 0;

    cout << setw(4)  << k
         << setw(12) << x
         << setw(12) << f(x)
         << setw(12) << df(x)
         << setw(12) << d2f(x)
         << setw(8)  << "-"
         << setw(12) << "-"
         << setw(12) << "-"
         << setw(12) << "-"
         << setw(12) << "-"
         << "\n";
    Nf += 3;

    bool done = false;
    while (!done && k < 50) {
        double d1 = df(x);
        double d2 = d2f(x);

        if (fabs(d1) <= eps) {
            cout << "|f'| <= eps -- STOP\n";
            break;
        }
        if (d2 <= 1e-12) {
            cout << "d2 <= 0 -- stop\n";
            break;
        }

        double t = 1.0;
        bool accepted = false;
        double xBar = x, fBar = f(x), d1Bar = d1;

        while (true) {
            xBar = x - t * d1 / d2;
            d1Bar = df(xBar);
            fBar = f(xBar);
            Nf += 3;

            double rhs = f(x) - alpha * t * d1 * d1 / d2;

            cout << setw(4)  << k + 1
                 << setw(12) << x
                 << setw(12) << f(x)
                 << setw(12) << d1
                 << setw(12) << d2
                 << setw(8)  << t
                 << setw(12) << xBar
                 << setw(12) << fBar
                 << setw(12) << d1Bar
                 << setw(12) << (fBar <= rhs ? "OK" : "reject")
                 << "\n";

            if (fabs(d1Bar) <= eps) {
                x = xBar;
                d1 = d1Bar;
                accepted = true;
                done = true;
                break;
            }

            if (fBar <= rhs) {
                x = xBar;
                accepted = true;
                break;
            }

            t /= 2.0;
            if (t < 1e-8) {
                cout << "t too small -- stop\n";
                done = true;
                break;
            }
        }

        if (!accepted && !done) break;
        k++;
    }

    double xStar = x;
    double fStar = f(x);
    double xExact = 3.1995;
    double err = fabs(xStar - xExact);
;
    cout << "    x* = " << xStar << "\n";
    cout << "    f(x*)= " << fStar << "\n";
    cout << "    Iterations= " << k << "\n";
    cout << "    Nf= " << Nf << "\n";
    cout << "    |x*-x*_t| = " << err << "\n";

    return {"Newton-Raphson", k, Nf, xStar, fStar, err};
}

int main() 
{
    cout << fixed << setprecision(6);

    double x0    = 2.70;
    double eps   = 0.001;
    double alpha = 0.5;

    cout << "f(x)  = e^(x-1) - 1 - x - x^2/2 - x^3/6\n";
    cout << "f'(x) = e^(x-1) - 1 - x - x^2/2\n";
    cout << "f''(x)= e^(x-1) - 1 - x\n\n";
    cout << "Interval from Swann: [2.70, 3.20]\n";
    cout << "x0 = " << x0 << ", eps = " << eps << "\n";

    MethodResult rNewton = newtonClassic(x0, eps);
    MethodResult rNR = newtonRaphson(x0, eps, alpha);
    MethodResult rDich  = {"Dichotomy",      10, 20, 3.199700, -5.756231, 0.000200};
    MethodResult rHalf  = {"Half Division",  10, 20, 3.199700, -5.756231, 0.000200};
    MethodResult rGold  = {"Golden Section", 13, 15, 3.199520, -5.756178, 0.000020};
    MethodResult rFib   = {"Fibonacci",      13, 15, 3.198770, -5.755954, 0.000730};
    MethodResult rPow   = {"Powell",          4,  7, 3.199700, -5.756167, 0.000200};

    MethodResult results[7] = {
        rNewton, rNR,
        rDich, rHalf, rGold, rFib, rPow
    };

    cout << "| " << left << setw(20) << "Method"
         << "| " << setw(6)  << "k"
         << "| " << setw(6)  << "Nf"
         << "| " << setw(12) << "x*"
         << "| " << setw(12) << "f(x*)"
         << "| " << setw(12) << "|x*-x*_t|"
         << "|\n";

    cout << "|" << string(21,'-') << "|" << string(7,'-')
         << "|" << string(7,'-')  << "|" << string(13,'-')
         << "|" << string(13,'-') << "|" << string(13,'-')
         << "|\n";

    for (int i = 0; i < 7; i++) {
        cout << "| " << left << setw(20) << results[i].name
             << "| " << setw(6)  << results[i].iterations
             << "| " << setw(6)  << results[i].Nf
             << "| " << setw(12) << results[i].xStar
             << "| " << setw(12) << results[i].fStar
             << "| " << setw(12) << results[i].error
             << "|\n";
    }

    return 0;
}