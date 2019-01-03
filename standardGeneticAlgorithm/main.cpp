#include <bits/stdc++.h>
#define MAXCANDIDATESIZE 900

using namespace std;

const double M_PI = acos(-1);

const int dim = 2;
const double leftMargin = -2, rightMargin = 2;

const int precision = 5;
const int bitLen = ceil(log2((rightMargin - leftMargin) * pow(10, precision)));

const int popSize = 100; //"standard" parameters
const double mutationProbability = 0.01, crossOverProbability = 0.25;

struct population {
    bool candidateSol[popSize][MAXCANDIDATESIZE];
    double candidateFitness[popSize], popFitness, wheelStart[popSize];
};
population *currPop = new population();

double Rastrigin(double *X);
double RastriginFitness(double *X);
double Griewangk(double *X);
double GriewangkFitness(double *X);
double Rosenbrock(double *X);
double RosenbrockFitness(double *X);
double SixHump(double *X);
double SixHumpFitness(double *X);

void initializePopulation();
double *decode(bool *X);

void evaluatePopulation(double (*f)(double *));
int binarySearch(int left, int right, double val);
void selectPopulation();

void mutatePopulation();
void crossOverPopulation();
void uniformCrossOverPopulation();

double bestFitness;
int noChange, bestCandidate;
void geneticAlgorithm (double (*f)(double *)) { //f - fitness function
    initializePopulation();
    evaluatePopulation(f);

    noChange = 0;
    bestFitness = currPop->candidateFitness[bestCandidate];

    /*termination condition
    reaching a maximum of iMax iterations/ no improvement in the last noChangeMax iterations*/
    for (int i = 1; i <= 1000 && noChange < 100; ++i) {
        selectPopulation(); //select P(i) from P(i-1)

        //alter P(i)
        mutatePopulation();
        crossOverPopulation();
        //uniformCrossOverPopulation();

        evaluatePopulation(f);
        if (currPop->candidateFitness[bestCandidate] <= bestFitness) ++noChange;
        else bestFitness = currPop->candidateFitness[bestCandidate], noChange = 0;
    }
}

int main()
{
    srand(time(0));
    int runs = 30;

    double curr, mean = 0, best = 1e9, worst = -1e9;
    double optim[runs], sd = 0;

    for (int i = 0; i < runs; ++i) {
            geneticAlgorithm(SixHumpFitness);

            double *decodeCandidate = decode(currPop->candidateSol[bestCandidate]);
            curr = SixHump(decodeCandidate);
            delete []decodeCandidate;

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



double Rastrigin(double *X) {
    double sum = 10 * dim;
    for (int i = 0; i < dim; ++i) sum += X[i]*X[i] - 10*cos(2 * M_PI * X[i]);
    return sum;
}

double RastriginFitness(double *X) {
    double sum = 10 * dim;
    for (int i = 0; i < dim; ++i) sum += X[i]*X[i] - 10*cos(2 * M_PI * X[i]);
    return 1 / (sum + 0.1);
}

double Griewangk(double *X) {
    double sum = 0, prod = 1;
    for (int i = 0; i < dim; ++i) sum += X[i]*X[i] / 4000;
    for (int i = 0; i < dim; ++i) prod *= cos(X[i] / sqrt(i+1));
    return sum - prod + 1;
}

double GriewangkFitness(double *X) {
    double sum = 0, prod = 1;
    for (int i = 0; i < dim; ++i) sum += X[i]*X[i] / 4000;
    for (int i = 0; i < dim; ++i) prod *= cos(X[i] / sqrt(i+1));
    return 1 / (sum - prod + 1 + 0.1);
}

double Rosenbrock(double *X) {
    double sum = 0;
    for (int i = 0; i < dim; ++i) sum += 100 * (X[i+1] - X[i]*X[i]) * (X[i+1] - X[i]*X[i]) + (1 - X[i]) * (1 - X[i]);
    return sum;
}

double RosenbrockFitness(double *X) {
    double sum = 0;
    for (int i = 0; i < dim; ++i) sum += 100 * (X[i+1] - X[i]*X[i]) * (X[i+1] - X[i]*X[i]) + (1 - X[i]) * (1 - X[i]);
    return 1 / (sum + 0.1);
}

double SixHump(double *X) {
    return (4 - 2.1 * X[0]*X[0] + X[0]*X[0]*X[0]*X[0] / 3) * X[0]*X[0] + X[0]*X[1] + (-4 + 4 * X[1]*X[1]) * X[1]*X[1];
}

double SixHumpFitness(double *X) {
    return 1 / (((4 - 2.1 * X[0]*X[0] + X[0]*X[0]*X[0]*X[0] / 3) * X[0]*X[0] + X[0]*X[1] + (-4 + 4 * X[1]*X[1]) * X[1]*X[1]) + 1.1);
}



void initializePopulation() {
    for (int i = 0; i < popSize; ++i)
        for (int j = 0; j < dim*bitLen; ++j) currPop->candidateSol[i][j] = rand() & 1;
}

double *decode(bool *X) {
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



void evaluatePopulation(double (*f)(double *)) {
    double *decodeCandidate = decode(currPop->candidateSol[0]);
    currPop->candidateFitness[0] = f(decodeCandidate);
    delete []decodeCandidate;

    currPop->popFitness += currPop->candidateFitness[0];
    bestCandidate = 0;

    //prepare the "fortune wheel"
    for (int i = 1; i < popSize; ++i) {
        //compute the fitness of the individual i
        double *decodeCandidate = decode(currPop->candidateSol[i]);
        currPop->candidateFitness[i] = f(decodeCandidate);
        delete []decodeCandidate;

        //compute the total fitness of the population
        currPop->popFitness += currPop->candidateFitness[i];

        if (currPop->candidateFitness[i] > currPop->candidateFitness[bestCandidate]) bestCandidate = i;
    }

    //set the starting point on the "wheel" for every individual
    for (int i = 1; i < popSize; ++i)
        currPop->wheelStart[i] = currPop->wheelStart[i-1] + //probability of selecting the individual i
        (currPop->candidateFitness[i-1] / currPop->popFitness);
}

int binarySearch(int left, int right, double val) {
    if (left > right) return -1;

    int mid = (left + right) / 2;
    if ((val >= currPop->wheelStart[mid] && val < currPop->wheelStart[mid+1]) || mid == popSize-1) return mid;

    if (val < currPop->wheelStart[mid]) return binarySearch (left, mid-1, val);
    else return binarySearch (mid+1, right, val);
}

void selectPopulation() {
    population *newPop = new population();
    for (int i = 0; i < popSize; ++i) {
        double fortune = rand() / (RAND_MAX + 1.);
        int luckyCandidate = binarySearch(0, popSize-1, fortune);
        memcpy(newPop->candidateSol[i], currPop->candidateSol[luckyCandidate], dim*bitLen * sizeof(bool));
    }

    delete currPop;
    currPop = newPop;
}



void mutatePopulation() {
    for (int i = 0; i < popSize; ++i)
        for (int j = 0; j < dim*bitLen; ++j)
            if (rand() / (RAND_MAX + 1.) < mutationProbability) currPop->candidateSol[i][j] = !currPop->candidateSol[i][j];
}

void crossOverPopulation() {
    int pos[2], curr = 1;
    for (int i = 0; i < popSize; ++i)
        if (rand() / (RAND_MAX + 1.) < crossOverProbability) {
            curr = (curr+1) % 2;
            pos[curr] = i;

            if (curr) {
                int cut = rand() % (dim*bitLen - 1) + 1;
                for (int j = cut; j < dim*bitLen; ++j) swap(currPop->candidateSol[pos[0]][j], currPop->candidateSol[pos[1]][j]);
            }
        }
}

void uniformCrossOverPopulation() {
    int pos[2], curr = 1;
    for (int i = 0; i < popSize; ++i)
        if (rand() / (RAND_MAX + 1.) < crossOverProbability) {
            curr = (curr+1) % 2;
            pos[curr] = i;

            if (curr) {
                for (int j = 0; j < dim*bitLen; ++j)
                    if (rand() / (RAND_MAX + 1.) < 0.5) swap(currPop->candidateSol[pos[0]][j], currPop->candidateSol[pos[1]][j]);
            }
        }
}

