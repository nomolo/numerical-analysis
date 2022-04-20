#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Newton
{
public:
    virtual void solve(const double &x0, const double &e1, const double &e2, int N) = 0;
    virtual double diff(const double &x) = 0;

private:
};

double Legendre(int n, double x)
{
    if (n == 0)
        return 1;
    if (n == 1)
        return x;

    return ((double)(2 * n - 1) / (double)(n)) * x * Legendre(n - 1, x) - ((double)(n - 1) / (double)(n)) * Legendre(n - 2, x);
}

double Laguerre(int n,double x)
{
     if (n == 0)
        return 1;
    if (n == 1)
        return x;

    return ((double)(2*n-1-x))*Laguerre(n-1,x)-((double)(n-1)*(n-1))*Laguerre(n-2,x);
}

double Hermite(int n,double x)
{
    if (n == 0)
        return 1;
    if (n == 1)
        return x;

    return 2*x*Hermite(n-1,x)-2*(n-1)*Hermite(n-2,x);
}

class herNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {

    }
    void solve(const double *x0, const double &e1, const double &e2, int N)
    {
        double F = 0, DF = 0;
        int n = 1;
        for (int i = 0; i < 6; i++)
        {
            double x2 = x0[i], x1 = x0[i];
            for (; n <= N; n++)
            {
                F = Hermite(6, x1);
                DF = diff(x1);
                if (abs(F) < e1)
                {
                    cout << "OK!n = "<<n << endl;
                    cout.precision(7);
                    cout << "x = " << x1 << endl;
                    break;
                }
                if (abs(DF) < e2)
                {
                    cout << "Error!The diff of f(x0) is too small." << endl;
                    break;
                }
                x2 = x1 - F / DF;
                if (abs(x2 - x1) < e1)
                {
                    cout << "OK!n = "<<n << endl;
                    cout.precision(7);
                    cout << "x = " << x2 << endl;
                    break;
                }
                x1 = x2;
            }
            if (n >= N)
                cout << "Failed to solve." << endl;
        }
    }

    double diff(const double &x)
    {
        return (Hermite(6, x + 0.00001) - Hermite(6, x)) / 0.00001;
    }
};

class lagNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {

    }
    void solve(const double *x0, const double &e1, const double &e2, int N)
    {
        double F = 0, DF = 0;
        int n = 1;
        for (int i = 0; i < 5; i++)
        {
            double x2 = x0[i], x1 = x0[i];
            for (; n <= N; n++)
            {
                F = Laguerre(5, x1);
                DF = diff(x1);
                if (abs(F) < e1)
                {
                    cout << "OK!n = "<<n << endl;
                    cout.precision(7);
                    cout << "x = " << x1 << endl;
                    break;
                }
                if (abs(DF) < e2)
                {
                    cout << "Error!The diff of f(x0) is too small." << endl;
                    break;
                }
                x2 = x1 - F / DF;
                if (abs(x2 - x1) < e1)
                {
                    cout << "OK!n = "<<n << endl;
                    cout.precision(7);
                    cout << "x = " << x2 << endl;
                    break;
                }
                x1 = x2;
            }
            if (n >= N)
                cout << "Failed to solve." << endl;
        }
    }

    double diff(const double &x)
    {
        return (Laguerre(5, x + 0.00001) - Laguerre(5, x)) / 0.00001;
    }
};

class legNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {

    }
    void solve(const double *x0, const double &e1, const double &e2, int N)
    {
        double F = 0, DF = 0;
        int n = 1;
        for (int i = 0; i < 6; i++)
        {
            double x2 = x0[i], x1 = x0[i];
            for (; n <= N; n++)
            {
                F = Legendre(6, x1);
                DF = diff(x1);
                if (abs(F) < e1)
                {
                    cout << "OK!n = "<<n << endl;
                    cout.precision(7);
                    cout << "x = " << x1 << endl;
                    break;
                }
                if (abs(DF) < e2)
                {
                    cout << "Error!The diff of f(x0) is too small." << endl;
                    break;
                }
                x2 = x1 - F / DF;
                if (abs(x2 - x1) < e1)
                {
                    cout << "OK!n = "<<n<< endl;
                    cout.precision(7);
                    cout << "x = " << x2 << endl;
                    break;
                }
                x1 = x2;
            }
            if (n >= N)
                cout << "Failed to solve." << endl;
        }
    }

    double diff(const double &x)
    {
        return (Legendre(6, x + 0.00001) - Legendre(6, x)) / 0.00001;
    }
};

class cosNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {
        double F, DF;
        double x2 = x0, x1 = x0;
        int n = 1;
        for (; n <= N; n++)
        {
            F = cos(x1) - x1;
            DF = diff(x1);
            if (abs(F) < e1)
            {
                cout << "OK!" << endl;
                cout << "x = " << x1 << endl;
                break;
            }
            if (abs(DF) < e2)
            {
                cout << "Error!The diff of f(x0) is too small." << endl;
                break;
            }
            x2 = x1 - F / DF;
            if (abs(x2 - x1) < e1)
            {
                cout << "OK!" << endl;
                cout.precision(7);
                cout << "x = " << x2 << endl;
                break;
            }
            x1 = x2;
        }
        if (n >= N)
            cout << "Failed to solve." << endl;
    }

    double diff(const double &x)
    {
        return (-sin(x) - 1);
    }
};

class exNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {
        double F = 0, DF = 0;
        double x2 = x0, x1 = x0;
        int n = 1;
        for (; n <= N; n++)
        {
            F = x1 - exp(-x1);
            DF = diff(x1);
            if (abs(F) < e1)
            {
                cout << "OK!" << endl;
                cout << "x = " << x1 << endl;
                break;
            }
            if (abs(DF) < e2)
            {
                cout << "Error!The diff of f(x0) is too small." << endl;
                break;
            }
            x2 = x1 - F / DF;
            if (abs(x2 - x1) < e1)
            {
                cout << "OK!" << endl;
                cout.precision(7);
                cout << "x = " << x2 << endl;
                break;
            }
            x1 = x2;
        }
        if (n >= N)
            cout << "Failed to solve." << endl;
    }

    double diff(const double &x)
    {
        return (1 + exp(-x));
    }
};

class exsinNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {
        double F = 0, DF = 0;
        double x2 = x0, x1 = x0;
        int n = 1;
        for (; n <= N; n++)
        {
            F = exp(-x1) - sin(x1);
            DF = diff(x1);
            if (abs(F) < e1)
            {
                cout << "OK!" << endl;
                cout << "x = " << x1 << endl;
                break;
            }
            if (abs(DF) < e2)
            {
                cout << "Error!The diff of f(x0) is too small." << endl;
                break;
            }
            x2 = x1 - F / DF;
            if (abs(x2 - x1) < e1)
            {
                cout << "OK!" << endl;
                cout.precision(7);
                cout << "x = " << x2 << endl;
                break;
            }
            x1 = x2;
        }
        if (n >= N)
            cout << "Failed to solve." << endl;
    }

    double diff(const double &x)
    {
        return (-exp(-x) - cos(x));
    }
};

class exxNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {
        double F = 0, DF = 0;
        double x2 = x0, x1 = x0;
        int n = 1;
        for (; n <= N; n++)
        {
            F = x1 * x1 - 2 * x1 * exp(-x1) + exp(-2 * x1);
            DF = diff(x1);
            if (abs(F) < e1)
            {
                cout << "OK!" << endl;
                cout << "x = " << x1 << endl;
                break;
            }
            if (abs(DF) < e2)
            {
                cout << "Error!The diff of f(x0) is too small." << endl;
                break;
            }
            x2 = x1 - F / DF;
            if (abs(x2 - x1) < e1)
            {
                cout << "OK!" << endl;
                cout.precision(7);
                cout << "x = " << x2 << endl;
                break;
            }
            x1 = x2;
        }
        if (n >= N)
            cout << "Failed to solve." << endl;
    }

    double diff(const double &x)
    {
        return (2 * x - 2 * exp(-x) + 2 * x * exp(-x) - 2 * exp(-2 * x));
    }
};

int main()
{
    double e1 = 1e-6, e2 = 1e-4;

    // cosNewton mycosslove;
    // exNewton myexsolve;
    // exsinNewton myexsinsolve;
    // exxNewton myexxsolve;
    //legNewton mylegsolve;
    lagNewton mylagsolve;
    herNewton myhersolve;

    double legx0[6]={-0.9,-0.2,0.2,0.9,0.6,-0.6};
    double lagx0[5]={0.4,1,3.0,7.0,12.0};
    double herx0[6]={-4.0,-1.6,-0.6,1.3,2.0,3.0};
    // myexxsolve.solve(0.5, e1, e2, 20);
    // myexsinsolve.solve(0.6, e1, e2, 10);
    // mycosslove.solve(M_PI / 4, e1, e2, 10);
    // myexsolve.solve(0.5, e1, e2, 10);
    //mylegsolve.solve(legx0, e1, e1, 20);
    //mylagsolve.solve(lagx0,e1,e1,500);
    myhersolve.solve(herx0,e1,e1,200);


    return 0;
}
