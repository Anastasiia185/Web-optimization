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
    double error;
};

MethodResult powell(double x1, double x2, double x3,
                    double sf, double sx) {
    double f1 = f(x1);
    double f2 = f(x2);
    double f3 = f(x3);
    int Nf = 3; 

    cout << setw(4)  << "k"
         << setw(11) << "x1"
         << setw(11) << "x2"
         << setw(11) << "x3"
         << setw(11) << "f(x1)"
         << setw(11) << "f(x2)"
         << setw(11) << "f(x3)"
         << setw(11) << "a1"
         << setw(11) << "a2"
         << setw(11) << "x_bar"
         << setw(11) << "f(x_bar)"
         << setw(11) << "|fmin-f|"
         << setw(11) << "|xmin-x|"
         << setw(12) << "Decision"
         << "\n";
    cout << string(155, '-') << "\n";

    double xBar = x2;
    double fBar = f2;

    for (int k = 1; k <= 51; k++) {

        int iMin = 0, iMax = 0;
        double xs[3] = { x1, x2, x3 };
        double fs[3] = { f1, f2, f3 };
        for (int i = 1; i < 3; i++) {
            if (fs[i] < fs[iMin]) iMin = i;
            if (fs[i] > fs[iMax]) iMax = i;
        }

        double a1 = (f2 - f1) / (x2 - x1);
        double a2 = ((f3 - f1) / (x3 - x1) - a1) / (x3 - x2);
        if (fabs(x2 - x1) < 1e-12 || fabs(x3 - x2) < 1e-12 || fabs(a2) < 1e-12) {
            cout << "Division by zero — stop. Change dx.\n";
            break;
        }
        if (a2 > 0) {
            xBar = -a1 / (2.0 * a2) + (x1 + x2) / 2.0;
        } else {
            xBar = 2.0 * xs[iMin] - xs[iMax];
        }
        fBar = f(xBar);
        Nf++;

        double dF = fabs(fs[iMin] - fBar);
        double dX = fabs(xs[iMin] - xBar);
        string decision;
        bool converged = false;

        if (dF <= sf && dX <= sx) {
            decision = "both < sigma -> STOP";
            converged = true;
        } else {
            if (fabs(xBar - xs[iMin]) > 1e-12) {
                if (iMax == 0)      { x1 = xBar; f1 = fBar; decision = "x1 <- xbar"; }
                else if (iMax == 1) { x2 = xBar; f2 = fBar; decision = "x2 <- xbar"; }
                else                { x3 = xBar; f3 = fBar; decision = "x3 <- xbar"; }
            } else {
                decision = "xbar==xmin — skip";
            }
        }

        cout << setw(4)  << k
             << setw(11) << x1
             << setw(11) << x2
             << setw(11) << x3
             << setw(11) << f1
             << setw(11) << f2
             << setw(11) << f3
             << setw(11) << a1
             << setw(11) << a2
             << setw(11) << xBar
             << setw(11) << fBar
             << setw(11) << dF
             << setw(11) << dX
             << setw(12) << decision
             << "\n";

        if (converged) break;
    }

    double xStar = xBar;
    double fStar = fBar;

    cout << "\nResult:\n";
    cout << "    x* = " << xStar << "\n";
    cout << "    f(x*)= " << fStar << "\n";
    cout << "    Nf= " << Nf << "\n";

    double xExact = 3.1995;
    double err = fabs(xStar - xExact);
    cout << "   |x*-x*_t| = " << err << "\n";

    return {"Powell", 0, Nf, xStar, fStar, err};
}

int main() 
{
    cout << fixed << setprecision(6);

    double a0  = 2.70;
    double b0  = 3.20;
    double dx  = 0.05;
    double sf  = 0.01;
    double sx  = 0.001;

    double x1 = a0;
    double x2 = a0 + dx;
    double x3 = a0 + 2.0 * dx;

    cout << "f(x) = e^(x-1) - 1 - x - x^2/2 - x^3/6\n";
    cout << "Interval: [" << a0 << ", " << b0 << "]\n";
    cout << "dx = " << dx << "  =>  x1 = " << x1
         << ", x2 = " << x2 << ", x3 = " << x3 << "\n\n";

    MethodResult rPowell = powell(x1, x2, x3, sf, sx);

    MethodResult rDich  = {"Dichotomy",     10, 20, 3.199700, -5.756231, 0.000200};
    MethodResult rHalf  = {"Half Division", 10, 20, 3.199700, -5.756231, 0.000200};
    MethodResult rGold  = {"Golden Section",13, 15, 3.199520, -5.756178, 0.000020};
    MethodResult rFib   = {"Fibonacci",     13, 15, 3.198770, -5.755954, 0.000730};

    MethodResult results[5] = { rPowell, rDich, rHalf, rGold, rFib };
    cout << "| " << left << setw(18) << "Method"
         << "| " << setw(6)  << "k"
         << "| " << setw(6)  << "Nf"
         << "| " << setw(12) << "x*"
         << "| " << setw(12) << "f(x*)"
         << "| " << setw(12) << "|x*-x*_t|"
         << "|\n";

    cout << "|" << string(19,'-') << "|" << string(7,'-')
         << "|" << string(7,'-')  << "|" << string(13,'-')
         << "|" << string(13,'-') << "|" << string(13,'-')
         << "|\n";

    for (int i = 0; i < 5; i++) {
        cout << "| " << left << setw(18) << results[i].name
             << "| " << setw(6)  << results[i].iterations
             << "| " << setw(6)  << results[i].Nf
             << "| " << setw(12) << results[i].xStar
             << "| " << setw(12) << results[i].fStar
             << "| " << setw(12) << results[i].error
             << "|\n";
    }
    return 0;
}