#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <ctime>

#define _USE_MATH_DEFINES

using namespace std;

const int dim = 30, iterations = 100;
const double leftMargin = -5.12, rightMargin = 5.12;

int precision = 5;
const int bitLen = ceil(log2((rightMargin - leftMargin) * pow(10, precision)));



double DeJong (double *X) {
    double sum = 0;
    for (int i = 0; i < dim; ++i) sum += X[i] * X[i];
    return sum;
}

double Schwefel (double *X) {
    double sum = 0;
    for (int i = 0; i < dim; ++i) sum += (-X[i])*sin(sqrt(abs(X[i])));
    return sum;
}

double Rastrigin (double *X) {
    double sum = 10 * dim;
    for (int i = 0; i < dim; ++i) sum += X[i]*X[i] - 10*cos(2 * M_PI * X[i]);
    return sum;
}

double SixHump(double *X) {
    return (4 - 2.1 * X[0] * X[0] + X[0] * X[0] * X[0] * X[0] / 3) * X[0] * X[0] + X[0] * X[1] + (-4 + 4 * X[1] * X[1]) * X[1] * X[1];
}



void randomExploration (bool *X) {
    for (int j = 0; j < dim * bitLen; ++j) X[j] = rand() & 1;
}

double *decode (bool *X) {
    double *realValue = new double[dim]();

    for (int i = 0; i < dim; ++i) {
        int64_t scaleFactor = 0;

        for (int j = 0; j < bitLen; ++j) {
            scaleFactor *= 2;
            scaleFactor += *(X + i*bitLen + j);
        }

        realValue[i] = leftMargin + (double)scaleFactor / ((1LL << bitLen) - 1) * (rightMargin - leftMargin);
    }

    return realValue;
}



void BIHC_Exploitation (bool *X, double (*f)(double *), bool *bestX, double &bestY) {
    for (int i = 0; i < dim*bitLen; ++i)
    {
        X[i] = !X[i];

        double *decX = decode(X);
        double currY = f(decX);
        delete []decX;

        //retin in bestY "best improvement neighbour"
        if (currY < bestY) bestY = currY, memcpy (bestX, X, dim*bitLen);

        X[i] = !X[i];
    }
}

double BIHC (double (*f)(double *)) {
    double best = 1e9, Y, bestY;
    bool X[dim * bitLen], bestX[dim * bitLen];

    for (int i = 1; i <= iterations; ++i) {
        //generez aleatoriu punctul de start
        randomExploration(X);

        //retin valoarea din acest punct
        double *decX = decode(X);
        Y = bestY = f(decX);
        delete []decX;

        bool local = false;
        //ma deplasez prin spatiul vecinilor la distanta Hamming 1 cat timp acestia imi ofera o sol mai buna
        do {
            BIHC_Exploitation (X, f, bestX, bestY);

            if (bestY < Y) Y = bestY, memcpy (X, bestX, dim * bitLen);
            else local = true;
        } while (!local);

        //verific daca optimul iteratiei curente imbunatateste optimul global
        best = min(best, Y);
    }

    return best;
}



bool foundBetter;
void FIHC_Exploitation (bool *X, double (*f)(double *), bool *bestX, double &bestY) {
    for (int i = 0; i < dim*bitLen && !foundBetter; ++i)
    {
        X[i] = !X[i];

        double *decX = decode(X);
        double currY = f(decX);
        delete []decX;

        //retin in bestY "first improvement neighbour" si setez conditia de oprire a iteratiei prin vecini
        if (currY < bestY) bestY = currY, memcpy (bestX, X, dim*bitLen), foundBetter = true;

        X[i] = !X[i];
    }
}

double FIHC (double (*f)(double *)) {
    double best = 1e9, Y, bestY;
    bool X[dim * bitLen], bestX[dim * bitLen];

    for (int i = 1; i <= iterations; ++i) {
        //generez aleatoriu punctul de start
        randomExploration(X);

        //retin valoarea din acest punct
        double *decX = decode(X);
        Y = bestY = f(decX);
        delete []decX;

        //ma deplasez prin spatiul vecinilor la distanta Hamming 1 cat timp acestia imi ofera o sol mai buna
        do {
            foundBetter = false;
            FIHC_Exploitation (X, f, bestX, bestY);

            if (bestY < Y) Y = bestY, memcpy (X, bestX, dim * bitLen);
        } while (foundBetter);

        //verific daca optimul iteratiei curente imbunatateste optimul global
        best = min(best, Y);
    }

    return best;
}



double T;
void SA_Exploitation (bool *X, double (*f)(double *), bool *bestX, double &bestY) {
    for (int i = 0; i < dim*bitLen && !foundBetter; ++i)
    {
        X[i] = !X[i];

        double *decX = decode(X);
        double currY = f(decX);
        delete []decX;

        //retin in bestY "first improvement neighbour" si setez conditia de oprire a iteratiei prin vecini
        if (currY < bestY) bestY = currY, memcpy (bestX, X, dim*bitLen), foundBetter = true;
        /*retin in bestY o solutie candidat mai slaba (cu speranta ca ma va conduce spre un optim local neexplorat,
        poate chiar cel global) daca trece testul urmator si setez conditia de oprire a iteratiei prin vecini*/
        else if ((rand() / (RAND_MAX + 1e-6) * 1.0) < exp(-abs(currY - bestY) / T))
            bestY = currY, memcpy (bestX, X, dim*bitLen), foundBetter = true;

        X[i] = !X[i];
    }
}

double SA (double (*f)(double *)) {
    double best = 1e9, Y, bestY;
    bool X[dim * bitLen], bestX[dim * bitLen];

    //initializez temperatura
    double T = 5;

    for (int i = 1; i <= iterations; ++i) {
        //generez aleatoriu punctul de start
        randomExploration(X);

        //retin valoarea din acest punct
        double *decX = decode(X);
        Y = bestY = f(decX);
        delete []decX;

        //ma deplasez prin spatiul vecinilor la distanta Hamming 1 cat timp acestia imi ofera o sol mai "buna"
        do {
            foundBetter = false;
            SA_Exploitation (X, f, bestX, bestY);

            if (bestY < Y) Y = bestY, memcpy (X, bestX, dim * bitLen);
        } while (foundBetter);

        //verific daca optimul iteratiei curente imbunatateste optimul global
        best = min(best, Y);
        //asigur o scadere treptata a temperaturii dupa fiecare iteratie prin inmultirea cu o val subunitara
        T *= 0.9;
    }

    return best;
}

int main()
{
    srand(time(0));
    int runs = 30;

    double curr, mean = 0, best = 1e9, worst = -1e9;
    double optim[runs+1], sd = 0;

    for (int i = 1; i <= runs; ++i) {
        curr = SA(Rastrigin);

        mean += curr; optim[i] = curr;
        best = min(best, curr);
        worst = max(worst, curr);
    }

    mean /= runs;
    for (int i = 1; i <= runs; ++i) sd += (optim[i] - mean) * (optim[i] - mean);
    sd = sqrt(sd / (runs-1));

    cout << "best: " << fixed << best << "\n";
    cout << "worst: " << fixed << worst << "\n";
    cout << "mean: " << fixed << mean << "\n";
    cout << "standard deviation: " << fixed << sd << "\n";

    return 0;
}

/*void BIHC_Exploitation (int currDim, bool *X, double (*f)(double *), bool *bestX, double &bestY) {
    if (currDim == dim) {
        double *decX = decode(X);
        double currY = f(decX);
        delete []decX;

        if (currY < bestY) bestY = currY, memcpy (bestX, X, dim*bitLen);
        return;
    }

    BIHC_Exploitation (currDim+1, X, f, bestX, bestY);
    for (int i = 0; i < bitLen; ++i) {
        X[currDim*bitLen + i] = !X[currDim*bitLen + i];
        BIHC_Exploitation (currDim+1, X, f, bestX, bestY);
        X[currDim*bitLen + i] = !X[currDim*bitLen + i];
    }
}*/

/*void FIHC_Exploitation (int currDim, bool *X, double (*f)(double *), bool *bestX, double &bestY) {
    if (foundBetter) return;

    if (currDim == dim) {
        double *decX = decode(X);
        double currY = f(decX);
        delete []decX;

        if (currY < bestY) bestY = currY, memcpy (bestX, X, dim*bitLen), foundBetter = true;
        return;
    }

    FIHC_Exploitation (currDim+1, X, f, bestX, bestY);
    for (int i = 0; i < bitLen; ++i) {
        X[currDim*bitLen + i] = !X[currDim*bitLen + i];
        FIHC_Exploitation (currDim+1, X, f, bestX, bestY);
        X[currDim*bitLen + i] = !X[currDim*bitLen + i];
    }
}*/
