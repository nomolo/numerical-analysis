#include <iostream>
#include <cmath>

using namespace std;

//继承时需要自己写构造函数,输入区间[a,b],初值alpha,和取点数N

class RungeKutta
{
public:
    void printResult(double *x, double *y)
    {
        cout << "_________________________________________________" << endl;
        cout << "N = "<<mN<<"    计算结果：" << endl
             << endl;
        for (int i = 0; i < mN; i++)
        {
            cout << "x= " << x[i] << "\t | "
                 << "y = " << y[i] << endl;
        }
    }

    void solve(double *x, double *y)
    {
        double x0 = ma, y0 = malpha, x1, y1;
        double h = ((mb - ma) / (double)mN);
        double K1, K2, K3, K4;

        for (int i = 0; i < this->mN; i++)
        {
            K1 = func(x0, y0);
            K2 = func(x0 + h / 2.0, y0 + h * K1 / 2.0);
            K3 = func(x0 + h / 2.0, y0 + h * K2 / 2.0);
            K4 = func(x0 + h, y0 + h * K3);
            x1 = x0 + h;
            y1 = y0 + h * (K1 + 2 * K2 + 2 * K3 + K4) / 6.0;
            x[i] = x1;
            y[i] = y1;
            x0 = x1;
            y0 = y1;
        }
    }

private:
    virtual double func(const double &x, const double &y) = 0;

protected:
    double ma, mb, malpha;

public:
    int mN;
};

class rungekutta1_1 : public RungeKutta
{
public:
    rungekutta1_1(double a, double b, double alpha, int N)
    {
        this->ma = a;
        this->mb = b;
        this->mN = N;
        this->malpha = alpha;
    }

private:
    virtual double func(const double &x, const double &y)
    {
        return x + y;
    }
};

class rungekutta1_2 : public RungeKutta
{
public:
    rungekutta1_2(double a, double b, double alpha, int N)
    {
        this->ma = a;
        this->mb = b;
        this->mN = N;
        this->malpha = alpha;
    }

private:
    virtual double func(const double &x, const double &y)
    {
        return -y * y;
    }
};

class rungekutta2_1 : public RungeKutta
{
public:
    rungekutta2_1(double a, double b, double alpha, int N)
    {
        this->ma = a;
        this->mb = b;
        this->mN = N;
        this->malpha = alpha;
    }

private:
    virtual double func(const double &x, const double &y)
    {
        return (2 * y / x + x * x * exp(x));
    }
};

class rungekutta2_2 : public RungeKutta
{
public:
    rungekutta2_2(double a, double b, double alpha, int N)
    {
        this->ma = a;
        this->mb = b;
        this->mN = N;
        this->malpha = alpha;
    }

private:
    virtual double func(const double &x, const double &y)
    {
        return (y * y + y) / x;
    }
};

class rungekutta3_1 : public RungeKutta
{
public:
    rungekutta3_1(double a, double b, double alpha, int N)
    {
        this->ma = a;
        this->mb = b;
        this->mN = N;
        this->malpha = alpha;
    }

private:
    virtual double func(const double &x, const double &y)
    {
        return -20 * (y - x * x) + 2 * x;
    }
};

class rungekutta3_2 : public RungeKutta
{
public:
    rungekutta3_2(double a, double b, double alpha, int N)
    {
        this->ma = a;
        this->mb = b;
        this->mN = N;
        this->malpha = alpha;
    }

private:
    virtual double func(const double &x, const double &y)
    {
        return -20 * y + 20 * sin(x) + cos(x);
    }
};

class rungekutta3_3 : public RungeKutta
{
public:
    rungekutta3_3(double a, double b, double alpha, int N)
    {
        this->ma = a;
        this->mb = b;
        this->mN = N;
        this->malpha = alpha;
    }

private:
    virtual double func(const double &x, const double &y)
    {
        return -20 * (y - exp(x) * sin(x)) + exp(x) * (sin(x) + cos(x));
    }
};

int main()
{
    double x_rk[25], y_rk[25]; //用来储存结果数组
    rungekutta1_1 prob1_1(0.0, 1.0, -1.0, 5);
    rungekutta1_2 prob1_2(0.0, 1.0, 1.0, 5);
    rungekutta2_1 prob2_1(1.0, 3.0, 0.0, 5);
    rungekutta2_2 prob2_2(1.0, 3.0, -2.0, 5);
    rungekutta3_1 prob3_1(0.0, 1.0, 1.0 / 3.0, 5);
    rungekutta3_2 prob3_2(0.0, 1.0, 1.0, 5);
    rungekutta3_3 prob3_3(0.0, 1.0, 0.0, 5);

    cout << "+++++++++问题1:+++++++++" << endl;
    prob1_1.solve(x_rk, y_rk);
    prob1_1.printResult(x_rk, y_rk);
    prob1_1.mN = 10;
    prob1_1.solve(x_rk, y_rk);
    prob1_1.printResult(x_rk, y_rk);
    prob1_1.mN = 20;
    prob1_1.solve(x_rk, y_rk);
    prob1_1.printResult(x_rk, y_rk);

    prob1_2.solve(x_rk, y_rk);
    prob1_2.printResult(x_rk, y_rk);
    prob1_2.mN = 10;
    prob1_2.solve(x_rk, y_rk);
    prob1_2.printResult(x_rk, y_rk);
    prob1_2.mN = 20;
    prob1_2.solve(x_rk, y_rk);
    prob1_2.printResult(x_rk, y_rk);

    cout << "+++++++++问题2:+++++++++" << endl;

    prob2_1.solve(x_rk, y_rk);
    prob2_1.printResult(x_rk, y_rk);
    prob2_1.mN = 10;
    prob2_1.solve(x_rk, y_rk);
    prob2_1.printResult(x_rk, y_rk);
    prob2_1.mN = 20;
    prob2_1.solve(x_rk, y_rk);
    prob2_1.printResult(x_rk, y_rk);

    prob2_2.solve(x_rk, y_rk);
    prob2_2.printResult(x_rk, y_rk);
    prob2_2.mN = 10;
    prob2_2.solve(x_rk, y_rk);
    prob2_2.printResult(x_rk, y_rk);
    prob2_2.mN = 20;
    prob2_2.solve(x_rk, y_rk);
    prob2_2.printResult(x_rk, y_rk);

    cout << "+++++++++问题3:+++++++++" << endl;

    prob3_1.solve(x_rk, y_rk);
    prob3_1.printResult(x_rk, y_rk);
    prob3_1.mN = 10;
    prob3_1.solve(x_rk, y_rk);
    prob3_1.printResult(x_rk, y_rk);
    prob3_1.mN = 20;
    prob3_1.solve(x_rk, y_rk);
    prob3_1.printResult(x_rk, y_rk);

    prob3_2.solve(x_rk, y_rk);
    prob3_2.printResult(x_rk, y_rk);
    prob3_2.mN = 10;
    prob3_2.solve(x_rk, y_rk);
    prob3_2.printResult(x_rk, y_rk);
    prob3_2.mN = 20;
    prob3_2.solve(x_rk, y_rk);
    prob3_2.printResult(x_rk, y_rk);

    prob3_3.solve(x_rk, y_rk);
    prob3_3.printResult(x_rk, y_rk);
    prob3_3.mN = 10;
    prob3_3.solve(x_rk, y_rk);
    prob3_3.printResult(x_rk, y_rk);
    prob3_3.mN = 20;
    prob3_3.solve(x_rk, y_rk);
    prob3_3.printResult(x_rk, y_rk);

    return 0;
}
