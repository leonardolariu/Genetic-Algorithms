#include <bits/stdc++.h>

using namespace std;

FILE *fin = fopen("WesternSahara.txt", "r"); const int actualBest = 27603, citiesNum = 29;
//FILE *fin = fopen("Djibouti.txt", "r"); const int actualBest = 6656, citiesNum = 38;
//FILE *fin = fopen("Qatar.txt", "r"); const int actualBest = 9352, citiesNum = 194;

const int popSize = 1000;
const double mutationProbability = 0.1, crossOverProbability = 0.15;

struct point {
    double x, y;
} city[citiesNum];

double dist[citiesNum][citiesNum];

struct population {
    int candidateSol[popSize][citiesNum];
    double candidateFitness[popSize], popFitness, wheelStart[popSize];
}; population *currPop = new population();



struct cityRank {
    int city;
    double randomNum;
} cityList[citiesNum];

int compare (cityRank a, cityRank b)
{
    if (a.randomNum < b.randomNum) return 1;
    return 0;
}

void initializePopulation();
void evaluatePopulation();

int binarySearch(int, int, double);
void selectPopulation();

void crossOverPopulation();
void mutatePopulation();



double lastBestFitness;
int noChange, bestCandidate;
void geneticAlgorithm () {
    initializePopulation();
    evaluatePopulation();

    noChange = 0;
    lastBestFitness = currPop->candidateFitness[bestCandidate];

    /*termination condition
    no improvement in the last noChangeMax iterations*/
    while (noChange < 50) {
        selectPopulation(); //select P(i) from P(i-1)

        //alter P(i)
        crossOverPopulation();
        mutatePopulation();

        evaluatePopulation();
        if (currPop->candidateFitness[bestCandidate] <= lastBestFitness) ++noChange;
        else noChange = 0;
        lastBestFitness = currPop->candidateFitness[bestCandidate];
    }
}

int main()
{
    int useless;
    for (int i = 0; i < citiesNum; ++i)
        fscanf(fin, "%d %lf %lf", &useless, &city[i].x, &city[i].y);
    fclose(fin);

    for (int i = 0; i < citiesNum-1; ++i)
        for (int j = i+1; j < citiesNum; ++j)
            dist[i][j] = dist[j][i] = sqrt((city[i].x - city[j].x) * (city[i].x - city[j].x) +
                                           (city[i].y - city[j].y) * (city[i].y - city[j].y));

    srand(time(0));
    int runs = 1;

    double curr, mean = 0, best = 1e9, worst = -1e9;
    double optim[runs], sd = 0;

    for (int i = 0; i < runs; ++i) {
        geneticAlgorithm();

        curr = 1 / currPop->candidateFitness[bestCandidate];

        mean += curr; optim[i] = curr;
        best = min(best, curr);
        worst = max(worst, curr);
    }

    mean /= runs;
    for (int i = 0; i < runs; ++i) sd += (optim[i] - mean) * (optim[i] - mean);
    sd = sqrt(sd / (runs-1));

    printf("best: %f\n", best);
    printf("worst: %f\n", worst);
    printf("mean: %f\n", mean);
    printf("standard deviation: %f\n", sd);
    printf("actual best: %d\n", actualBest);

    return 0;
}

void initializePopulation() {
    for (int i = 0; i < popSize; ++i) {
        for (int j = 0; j < citiesNum; ++ j) {
            cityList[j].city = j;
            cityList[j].randomNum = rand() / (RAND_MAX + 1.);
        }

        sort(cityList, cityList + citiesNum, compare);

        for (int j = 0; j < citiesNum; ++j) currPop->candidateSol[i][j] = cityList[j].city;
    }
}

void evaluatePopulation() {
    bestCandidate = 0;
    //prepare the "fortune wheel"
    for (int i = 0; i < popSize; ++i) {
        //compute the fitness of the individual i
        for (int j = 0; j < citiesNum-1; ++j)
            currPop->candidateFitness[i] += dist[currPop->candidateSol[i][j]][currPop->candidateSol[i][j+1]];
        currPop->candidateFitness[i] += dist[currPop->candidateSol[i][citiesNum-1]][currPop->candidateSol[i][0]];

        //assure that shorter length determines greater fitness
        currPop->candidateFitness[i] = 1 / currPop->candidateFitness[i];

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
        memcpy(newPop->candidateSol[i], currPop->candidateSol[luckyCandidate], citiesNum * sizeof(int));
    }

    delete currPop;
    currPop = newPop;
}



void crossOverPopulation() {
    int parent[2], parentCount = 1;
    for (int i = 0; i < popSize; ++i)
        if (rand() / (RAND_MAX + 1.) < crossOverProbability) {
            parentCount = (parentCount+1) % 2;
            parent[parentCount] = i;

            if (parentCount) {
                int cut1 = rand() % (citiesNum), cut2;
                do cut2 = rand() % (citiesNum); while (cut1 == cut2);

                if (cut1 > cut2) swap(cut1, cut2);

                int parentCopy[citiesNum];
                memcpy(parentCopy, currPop->candidateSol[parent[0]], citiesNum * sizeof(int));

                bool used[citiesNum] = {};
                for (int j = cut1; j <= cut2; ++j) used[currPop->candidateSol[parent[0]][j]] = 1;

                int currPos = 0;
                for (int j = 0; j < citiesNum; ++j) {
                    if (currPos == cut1) currPos = cut2+1;
                    if (!used[currPop->candidateSol[parent[1]][j]]) {
                        currPop->candidateSol[parent[0]][currPos++] = currPop->candidateSol[parent[1]][j];
                        used[currPop->candidateSol[parent[1]][j]] = 1;
                    }
                }

                used[citiesNum] = {};
                for (int j = cut1; j <= cut2; ++j) used[currPop->candidateSol[parent[1]][j]] = 1;

                currPos = 0;
                for (int j = 0; j < citiesNum; ++j) {
                    if (currPos == cut1) currPos = cut2+1;
                    if (!used[parentCopy[j]]) {
                        currPop->candidateSol[parent[1]][currPos++] = parentCopy[j];
                        used[parentCopy[j]] = 1;
                    }
                }
            }
        }
}

void mutatePopulation() {
    for (int i = 0; i < popSize; ++i)
        if (rand() / (RAND_MAX + 1.) < mutationProbability) {
            int cut1 = rand() % (citiesNum), cut2;
            do cut2 = rand() % (citiesNum); while (cut1 == cut2);

            swap(currPop->candidateSol[i][cut1], currPop->candidateSol[i][cut2]);
        }
}
