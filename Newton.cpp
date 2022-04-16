#include <iostream>
#include <vector>
#include <cmath>
#include <boost/ptr_container/indirect_fun.hpp>

using namespace std;

class Newton
{
public:
    virtual void solve(const double &x0, const double &e1, const double &e2, int N) = 0;
    virtual double diff(const double &x) = 0;

private:
};

class cosNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {
        double F = cos(x0) - x0, DF = diff(x0);
        double x2 =x0, x1=x0;
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
                cout.precision(10);
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
        double F=0, DF=0 ;
        double x2 =x0, x1=x0;
        int n = 1;
        for (; n <= N; n++)
        {
            F = x1-exp(-x1);
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
                cout.precision(10);
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
        return (1+exp(-x));
    }
};

class exsinNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {
        double F=0, DF=0 ;
        double x2 =x0, x1=x0;
        int n = 1;
        for (; n <= N; n++)
        {
            F = exp(-x1)-sin(x1);
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
                cout.precision(10);
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
        return (-exp(-x)-cos(x));
    }
};

class exNewton : public Newton
{
public:
    void solve(const double &x0, const double &e1, const double &e2, int N)
    {
        double F=0, DF=0 ;
        double x2 =x0, x1=x0;
        int n = 1;
        for (; n <= N; n++)
        {
            F = x1*x1-2*x1*exp(-x1)+exp(-2*x1);
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
                cout.precision(10);
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
        return (2*x-2*exp(-x)+2*x*exp(-x)-2*exp(-2*x));
    }
};

int main()
{
    double x0 = M_PI/4;
    double e1 = 1e-6, e2 = 1e-4;
    int N = 10;

    cosNewton mycosslove;
    exNewton myexsolve;

    mycosslove.solve(x0, e1, e2, N);
    myexsolve.solve(0.5,e1,e2,10);

    return 0;
}
