#include <iostream>
#include <vector>
#include <cmath>
#include <string.h>
using namespace std;

template <class T>
class Lagrange
{
public:
    void lagrange(const int &n, const int &p, T *x, T *y, vector<T> &xi, vector<T> &fx) //输入一组已知点，计算未知点,p为未知点个数，n+1为已知点个数
    {
        for (int m = 0; m < p; m++)
        {
            for (int k = 0; k <= n; k++)
            {
                for (int j = 0; j <= n; j++)
                {
                    if (j != k)
                        l = l * (x[m] - xi.at(j)) / (xi.at(k) - xi.at(j));
                }
                y[m] = y[m] + l * fx.at(k);
                l = 1.0;
            }
        }
    }

private:
    double l = 1.0;
};

int main()
{
    int n = 5;

    vector<double> xi[4], fx[3];
    vector<double> ex[3];
    Lagrange<double> lag;

    for (int q = 0; q < 3; q++)
    {
        if (q == 0)
            n = 5;
        if (q == 1)
            n = 10;
        if (q == 2)
            n = 20;
        for (int i = 0; i <= n; i++)
        {
            xi[q].push_back((-5.0 + i * 10.0 / n));
            fx[q].push_back(1.0 / (1.0 + xi[q].at(i) * xi[q].at(i)));
        }
    }

    double y[4][5] = {}, x[5] = {0.75, 1.75, 2.75, 3.75, 4.75};

    cout << endl
         << endl
         << "问题1:" << endl;

    for (int p = 0; p < 3; p++)
    {
        if (p == 0)
            n = 5;
        if (p == 1)
            n = 10;
        if (p == 2)
            n = 20;
        lag.lagrange(n, 5, x, y[p], xi[p], fx[p]);
    }
    cout.precision(10);
    cout << "1/(1+x^2),插值区间[-5,5]" << endl;
    for (int t = 0; t < 3; t++)
    {
        cout << "____________________________________________" << endl;
        if (t == 0)
        {
            cout << "n = 5:" << endl;
        }
        if (t == 1)
        {
            cout << "n = 10:" << endl;
        }
        if (t == 2)
        {
            cout << "n = 20:" << endl;
        }
        for (int i = 0; i < 5; i++)
            std::cout << "x = " << x[i] << ", "
                      << "y = " << y[t][i] << std::endl;
    }
    cout << "===========================================";
    cout << endl
         << "真实值为" << endl;
    for (int i = 0; i < 5; i++)
    {
        cout << "x = " << x[i] << "y = " << (1.0 / (1.0 + x[i] * x[i])) << endl;
    }

    for (int i = 0; i < 3; i++)
    {
        xi[i].clear();
    }

    for (int i = 0; i < 5; i++)
    {
        x[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int v = 0; v < 5; v++)
            y[i][v] = 0;
    }

    for (int q = 0; q < 3; q++)
    {
        if (q == 0)
            n = 5;
        if (q == 1)
            n = 10;
        if (q == 2)
            n = 20;
        for (int i = 0; i <= n; i++)
        {
            xi[q].push_back((-1.0 + i * 2.0 / n));
            ex[q].push_back(exp(xi[q].at(i)));
        }
    }

    x[0] = -0.95;
    x[1] = -0.05;
    x[2] = 0.05;
    x[3] = 0.95;

    for (int p = 0; p < 3; p++)
    {
        if (p == 0)
            n = 5;
        if (p == 1)
            n = 10;
        if (p == 2)
            n = 20;
        lag.lagrange(n, 4, x, y[p], xi[p], ex[p]);
    }
    cout.precision(10);
    cout << endl
         << endl
         << "e^x,插值区间[-1,1]" << endl;
    for (int t = 0; t < 3; t++)
    {
        cout << "____________________________________________" << endl;
        if (t == 0)
        {
            cout << "n = 5:" << endl;
        }
        if (t == 1)
        {
            cout << "n = 10:" << endl;
        }
        if (t == 2)
        {
            cout << "n = 20:" << endl;
        }
        for (int i = 0; i < 4; i++)
            std::cout << "x = " << x[i] << ", "
                      << "y = " << y[t][i] << std::endl;
    }

    cout << "===========================================";
    cout << endl
         << "真实值为" << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << "x = " << x[i] << "y = " << (exp(x[i])) << endl;
        ;
    }

    cout << endl
         << endl
         << "问题2:" << endl;

    for (int i = 0; i < 3; i++)
    {
        xi[i].clear();
        fx[i].clear();
    }

    for (int i = 0; i < 5; i++)
    {
        x[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int v = 0; v < 5; v++)
            y[i][v] = 0;
    }

    x[0] = -0.95;
    x[1] = -0.05;
    x[2] = 0.05;
    x[3] = 0.95;

    for (int q = 0; q < 3; q++)
    {
        if (q == 0)
            n = 5;
        if (q == 1)
            n = 10;
        if (q == 2)
            n = 20;
        for (int i = 0; i <= n; i++)
        {
            xi[q].push_back((-1.0 + i * 2.0 / n));
            fx[q].push_back(1.0 / (1.0 + xi[q].at(i) * xi[q].at(i)));
        }
    }

    for (int p = 0; p < 3; p++)
    {
        if (p == 0)
            n = 5;
        if (p == 1)
            n = 10;
        if (p == 2)
            n = 20;
        lag.lagrange(n, 4, x, y[p], xi[p], fx[p]);
    }
    cout.precision(10);
    cout << endl
         << endl
         << "1/(1+x^2),插值区间[-1,1]" << endl;
    for (int t = 0; t < 3; t++)
    {
        cout << "____________________________________________" << endl;
        if (t == 0)
        {
            cout << "n = 5:" << endl;
        }
        if (t == 1)
        {
            cout << "n = 10:" << endl;
        }
        if (t == 2)
        {
            cout << "n = 20:" << endl;
        }
        for (int i = 0; i < 4; i++)
            std::cout << "x = " << x[i] << ", "
                      << "y = " << y[t][i] << std::endl;
    }

    cout << "===========================================";
    cout << endl
         << "真实值为" << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << "x = " << x[i] << "y = " << (1.0 / (1.0 + x[i] * x[i])) << endl;
    }

    for (int i = 0; i < 3; i++)
    {
        xi[i].clear();
        ex[i].clear();
    }
    for (int i = 0; i < 5; i++)
    {
        x[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int v = 0; v < 5; v++)
            y[i][v] = 0;
    }

    for (int q = 0; q < 3; q++)
    {
        if (q == 0)
            n = 5;
        if (q == 1)
            n = 10;
        if (q == 2)
            n = 20;
        for (int i = 0; i <= n; i++)
        {
            xi[q].push_back((-5.0 + i * 10.0 / n));
            ex[q].push_back(exp(xi[q].at(i)));
        }
    }

    x[0] = -4.75;
    x[1] = -0.25;
    x[2] = 0.25;
    x[3] = 4.75;

    for (int p = 0; p < 3; p++)
    {
        if (p == 0)
            n = 5;
        if (p == 1)
            n = 10;
        if (p == 2)
            n = 20;
        lag.lagrange(n, 4, x, y[p], xi[p], ex[p]);
    }
    cout.precision(10);
    cout << endl
         << endl
         << "e^x,插值区间[-5,5]" << endl;
    for (int t = 0; t < 3; t++)
    {
        cout << "____________________________________________" << endl;
        if (t == 0)
        {
            cout << "n = 5:" << endl;
        }
        if (t == 1)
        {
            cout << "n = 10:" << endl;
        }
        if (t == 2)
        {
            cout << "n = 20:" << endl;
        }
        for (int i = 0; i < 4; i++)
            std::cout << "x = " << x[i] << ", "
                      << "y = " << y[t][i] << std::endl;
    }

    cout << "===========================================";
    cout << endl
         << "真实值为" << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << "x = " << x[i] << "y = " << (exp(x[i])) << endl;
        ;
    }

    cout << endl
         << endl
         << "问题3:" << endl;

    for (int i = 0; i < 3; i++)
    {
        xi[i].clear();
        fx[i].clear();
    }
    for (int i = 0; i < 5; i++)
    {
        x[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int v = 0; v < 5; v++)
            y[i][v] = 0;
    }

    x[0] = -0.95;
    x[1] = -0.05;
    x[2] = 0.05;
    x[3] = 0.95;

    for (int q = 0; q < 3; q++)
    {
        if (q == 0)
            n = 5;
        if (q == 1)
            n = 10;
        if (q == 2)
            n = 20;
        for (int i = 0; i <= n; i++)
        {
            xi[q].push_back(cos(((2 * i + 1) * M_PI) / (2 * (n + 1))));
            fx[q].push_back(1.0 / (1.0 + xi[q].at(i) * xi[q].at(i)));
        }
    }

    for (int p = 0; p < 3; p++)
    {
        if (p == 0)
            n = 5;
        if (p == 1)
            n = 10;
        if (p == 2)
            n = 20;
        lag.lagrange(n, 4, x, y[p], xi[p], fx[p]);
    }
    cout.precision(10);
    cout << endl
         << endl
         << "1/(1+x^2),插值区间[-1,1]" << endl;
    for (int t = 0; t < 3; t++)
    {
        cout << "____________________________________________" << endl;
        if (t == 0)
        {
            cout << "n = 5:" << endl;
        }
        if (t == 1)
        {
            cout << "n = 10:" << endl;
        }
        if (t == 2)
        {
            cout << "n = 20:" << endl;
        }
        for (int i = 0; i < 4; i++)
            std::cout << "x = " << x[i] << ", "
                      << "y = " << y[t][i] << std::endl;
    }

    cout << "===========================================";
    cout << endl
         << "真实值为" << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << "x = " << x[i] << "y = " << (1.0 / (1.0 + x[i] * x[i])) << endl;
    }

    for (int i = 0; i < 3; i++)
    {
        xi[i].clear();
        ex[i].clear();
    }

    for (int i = 0; i < 5; i++)
    {
        x[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int v = 0; v < 5; v++)
            y[i][v] = 0;
    }

    for (int q = 0; q < 3; q++)
    {
        if (q == 0)
            n = 5;
        if (q == 1)
            n = 10;
        if (q == 2)
            n = 20;
        for (int i = 0; i <= n; i++)
        {
            xi[q].push_back(cos(((2 * i + 1) * M_PI) / (2 * (n + 1))));
            ex[q].push_back(exp(xi[q].at(i)));
        }
    }

    x[0] = -0.95;
    x[1] = -0.05;
    x[2] = 0.05;
    x[3] = 0.95;

    for (int p = 0; p < 3; p++)
    {
        if (p == 0)
            n = 5;
        if (p == 1)
            n = 10;
        if (p == 2)
            n = 20;
        lag.lagrange(n, 4, x, y[p], xi[p], ex[p]);
    }
    cout.precision(10);
    cout << endl
         << endl
         << "e^x,插值区间[-1,1]" << endl;
    for (int t = 0; t < 3; t++)
    {
        cout << "____________________________________________" << endl;
        if (t == 0)
        {
            cout << "n = 5:" << endl;
        }
        if (t == 1)
        {
            cout << "n = 10:" << endl;
        }
        if (t == 2)
        {
            cout << "n = 20:" << endl;
        }
        for (int i = 0; i < 4; i++)
            std::cout << "x = " << x[i] << ", "
                      << "y = " << y[t][i] << std::endl;
    }

    cout << "===========================================";
    cout << endl
         << "真实值为" << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << "x = " << x[i] << "y = " << (exp(x[i])) << endl;
        ;
    }

    memset(x, 0, sizeof(x));
    memset(y, 0, sizeof(y));
    for (int i = 0; i < 4; i++)
    {
        xi[i].clear();
    }

    cout << endl
         << endl
         << "问题4:" << endl;

    vector<double> sx[4];
    {
        xi[0].push_back(1.0);
        xi[0].push_back(4.0);
        xi[0].push_back(9.0);
    }
    {
        xi[1].push_back(36);
        xi[1].push_back(49);
        xi[1].push_back(64);
    }
    {
        xi[2].push_back(100);
        xi[2].push_back(121);
        xi[2].push_back(144);
    }
    {
        xi[3].push_back(169);
        xi[3].push_back(196);
        xi[3].push_back(225);
    }

    for (int q = 0; q < 4; q++)
    {
        for (int i = 0; i <= 2; i++)
        {
            sx[q].push_back(sqrt(xi[q].at(i)));
        }
    }

    x[0]=5;
    x[1]=50;
    x[2]=115;
    x[3]=185;

    for (int p = 0; p < 4; p++)
    {
        lag.lagrange(2, 4, x, y[p], xi[p], sx[p]);
    }

    cout.precision(10);
    cout << endl
         << endl
         << "sqrt(x)" << endl;
    for (int t = 0; t < 4; t++)
    {
        cout << "____________________________________________" << endl;
        for (int i = 0; i < 4; i++)
            std::cout << "x = " << x[i] << ", "
                      << "y = " << y[t][i] << std::endl;
    }

    cout << "===========================================";
    cout << endl
         << "真实值为" << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << "x = " << x[i] << "y = " << (sqrt(x[i])) << endl;
        ;
    }

    

    return 0;
}
