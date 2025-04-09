#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;
int nrCromozomi; // numarul de cromozomi
float capatSt; // capatul stang al intervalului
float capatDr;  // capatul drept al intervalului
int a, b, c, d; // coeficienti pentru functia de fitness
int precizie; // precizia pentru functia de fitness
int probCrossover; // probabilitatea de crossover
int probMutatie; // probabilitatea de mutatie
int nrGeneratii; // numarul de generatii
int lungimeCromozom;

ifstream in("input");
ofstream out("output");

void readFile()
{
    in >> nrCromozomi;
    in >> capatSt >> capatDr;
    in >> d >> a >> b >> c;
    in >> precizie;
    in >> probCrossover;
    in >> probMutatie;
    in >> nrGeneratii;
    
    lungimeCromozom = ceil(log2((capatDr-capatSt) * pow(10, precizie)));
}

typedef vector<int> Individ;
Individ genereazaIndivid(){
    Individ individ(lungimeCromozom);
    for (int i = 0; i < lungimeCromozom; ++i) {
        individ[i] = rand() % 2;
    }
    return individ;
}
double decodare(Individ individ) {
    float x = 0;
    for (int i = 0; i < lungimeCromozom; ++i) {
        x += individ[i] * pow(2, lungimeCromozom - i - 1);
    }
    return capatSt + (capatDr - capatSt) * x / (pow(2, lungimeCromozom) - 1);
}
// Overload <<
ostream& operator<<(ostream& os, const Individ& individ) {
    for (int i = 0; i < individ.size(); ++i) {
        os << individ[i];
    }
    return os;
}

void genereazaPopulatie(vector<Individ>& populatie){
    for (int i = 0; i < nrCromozomi; ++i) {
        populatie.push_back(genereazaIndivid());
    }
}

double fitnessFunction(double x)
{
    return d * pow(x, 4) + a * x * x + b * x + c;
}

void calculeazaProbabilitatiSelectie(const vector<Individ>& populatie, vector<double>& probabilitatiSelectie, double& performantaTotala)
{
    performantaTotala = 0.0;
    for (int i = 0; i < nrCromozomi; ++i) {
        double x = decodare(populatie[i]);
        probabilitatiSelectie[i] = fitnessFunction(x);
        performantaTotala += probabilitatiSelectie[i];
    }
    for (int i = 0; i < nrCromozomi; ++i) {
        probabilitatiSelectie[i] /= performantaTotala;
    }
}

void genereazaIntervaleSelectie(const vector<double>& probabilitatiSelectie, vector<double>& intervaleSelectie)
{
    intervaleSelectie.push_back(0);
    for (int i = 0; i < nrCromozomi; ++i) {
        intervaleSelectie.push_back(probabilitatiSelectie[i] + intervaleSelectie[i]);
    }
}

void printPopulatieInitiala(const vector<Individ>& populatie)
{
    out << "    Populatia initiala: \n";
    for (int i = 0; i < nrCromozomi; ++i, out << endl) {
        double x = decodare(populatie[i]);
        out << i+1 << ": " << populatie[i] << " ";
        out << "x= " << x << " ";
        out << "f(x)= " << fitnessFunction(x) << " ";
    }
    out << endl;
}

void printProbabilitatiSelectie(vector<double>& probabilitatiSelectie){
    out<<"\n\n    Probabilitati selectie: \n";
    for (int i = 0; i < nrCromozomi; ++i, out << endl) {
        out << "cromozom    " << i+1 << " probabilitate " << probabilitatiSelectie[i] << " ";
    }
}

void printIntervaleSelectie(vector<double>& intervaleSelectie){
    out<<"\n\n    Intervale selectie: \n";
    for(int i = 0; i < intervaleSelectie.size(); ++i){
        out << intervaleSelectie[i] << " ";
    }
    out << endl;
}

int binarySearch(const vector<double>& intervaleSelectie, double u) {
    int left = 0;
    int right = intervaleSelectie.size() - 1;
    while (left < right) {
        int mid = left + (right - left) / 2;
        if (u < intervaleSelectie[mid]) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }
    return left;
}

void procesSelectie(const vector<Individ>& populatie, const vector<double>& intervaleSelectie, vector<Individ>& indiviziSelectati)
{
    // selectie proportionala
    for(int i=0; i < nrCromozomi; ++i){
        double u = (double)rand() / RAND_MAX; // numar aleator in [0, 1]
        out << "u= " << u;
        int index = binarySearch(intervaleSelectie, u); // cautam intervalul in care se afla u
        if(index == nrCromozomi)
            index--; // daca u este 1, luam ultimul interval
        out << "    selectam cromozomul " << index + 1 << endl;
        indiviziSelectati.push_back(populatie[index]); // selectam cromozomul
        }

    out << "\n\n    Dupa selectie: \n";
    for (int i = 0; i < nrCromozomi; ++i, out << endl) {
        double x = decodare(indiviziSelectati[i]);
        Individ individSelectat = indiviziSelectati[i];
        out << i+1 << ": " << individSelectat << " ";
        out << "x= " << x << " ";
        out << "f(x)= " << fitnessFunction(x) << " ";
    }
}

void probabilitateIncrucisare(vector<Individ>& indiviziSelectati, vector<int>& indiviziIncrucisati)
{
    out << "\n\n    Probabilitatea de incrucisare: " << probCrossover / 100.0 << endl;
    for (int i = 0; i < nrCromozomi; ++i, out << endl) {
        double u = (double)rand() / RAND_MAX; // numar aleator in [0, 1]
        double x = decodare(indiviziSelectati[i]);
        out << i+1 << ": " << indiviziSelectati[i]<< " ";
        out << "u= " << u << " ";
        if (u < probCrossover / 100.0) {
            out<< "<" << probCrossover / 100.0 << " participa";
            indiviziIncrucisati.push_back(i); // adaugam individul la lista de indivizi incrucisati
        }
    }
    
    if(indiviziIncrucisati.size() % 2 != 0) {
        indiviziIncrucisati.pop_back(); // eliminam ultimul individ daca numarul de indivizi incrucisati este impar
    }


    for(int i = 0; i < indiviziIncrucisati.size(); i+=2) {
        out << "Recombinare dintre cromozomul " << indiviziIncrucisati[i] + 1 << " si cromozomul " << indiviziIncrucisati[i+1] + 1 << ": \n";
    
        int punctIncrucisare = rand() % lungimeCromozom; // punct de incrucisare
        Individ individ1 = indiviziSelectati[indiviziIncrucisati[i]];
        Individ individ2 = indiviziSelectati[indiviziIncrucisati[i+1]];
        out << individ1 << " " << individ2 << " punctIncrucisare= " << punctIncrucisare << endl;
    
        for(int j = 0; j < punctIncrucisare; j++) {
            swap(individ1[j], individ2[j]); // facem swap la cromozom
        }

        indiviziSelectati[indiviziIncrucisati[i]] = individ1; // actualizam individul 1
        indiviziSelectati[indiviziIncrucisati[i+1]] = individ2; // actualizam individul 2
        out << "Rezultat    "<< individ1 << " " << individ2 << endl;
        }
}

void printRecombinare(const vector<Individ>& indiviziSelectati)
{
    out << "\n\n    Dupa recombinare: \n";
    for (int i = 0; i < nrCromozomi; ++i, out << endl) {
        double x = decodare(indiviziSelectati[i]);
        Individ individSelectat = indiviziSelectati[i];
        out << i+1 << ": " << individSelectat << " ";
        out << "x= " << x << " ";
        out << "f(x)= " << fitnessFunction(x) << " ";
    }
}

void printMutatie(vector<Individ>& indiviziSelectati)
{
    out << "\n\n    Probabilitate de mutatie pentru fiecare gena: " << probMutatie/100.0 << " \n";
    out << "\n\n    Au fost modificati cromozomii: \n";
    // folosim mutatia rara
    for(int i = 0; i < nrCromozomi; ++i) {
        double u = (double)rand() / RAND_MAX; // numar aleator in [0, 1]
        if( u < probMutatie / 100.0) {
            int pozitie = rand() % lungimeCromozom; // pozitia la care se face mutatia
            indiviziSelectati[i][pozitie] = 1 - indiviziSelectati[i][pozitie]; // facem mutatia prin inversarea bitului
            out << i+1 << endl;
        }
    }
    out << "\n\n    Dupa mutatie: \n";
    for (int i = 0; i < nrCromozomi; ++i, out << endl) {
        double x = decodare(indiviziSelectati[i]);
        Individ individSelectat = indiviziSelectati[i];
        out << i+1 << ": " << individSelectat << " ";
        out << "x= " << x << " ";
        out << "f(x)= " << fitnessFunction(x) << " ";
    }
}

void printUrmatoareleGeneratii()
{
    out << "\n\n   Evolutia maximului";
}

int main()
{
    srand(time(0)); // Initializam generatorul cu timpul curent
    readFile();
    vector<Individ> populatie;

    genereazaPopulatie(populatie);
    printPopulatieInitiala(populatie);

    vector<double> probabilitatiSelectie(nrCromozomi);
    vector<double> fitness(nrCromozomi);
    double performantaTotala;
    calculeazaProbabilitatiSelectie(populatie, probabilitatiSelectie, performantaTotala); // e buna
    
    printProbabilitatiSelectie(probabilitatiSelectie); // e buna
    
    vector<double> intervaleSelectie;
    genereazaIntervaleSelectie(probabilitatiSelectie, intervaleSelectie); // e buna
    printIntervaleSelectie(intervaleSelectie);

    vector<Individ> indiviziSelectati;
    procesSelectie(populatie, intervaleSelectie, indiviziSelectati);

    vector<int> indiviziIncrucisati;
    probabilitateIncrucisare(indiviziSelectati, indiviziIncrucisati);

    printRecombinare(indiviziSelectati);

    printMutatie(indiviziSelectati);
    

    populatie = indiviziSelectati; // populatia devine populatia selectata
        
    out << "\n\n";
    // Generatii
    for(int i=1; i<nrGeneratii;++i)
    {
        vector<double> probabilitatiSelectie(nrCromozomi);
        vector<double> fitness(nrCromozomi);
        double performantaTotala;
        calculeazaProbabilitatiSelectie(populatie, probabilitatiSelectie, performantaTotala);

        vector<double> intervaleSelectie;
        genereazaIntervaleSelectie(probabilitatiSelectie, intervaleSelectie);

        vector<Individ> indiviziSelectati;
        for(int i=0; i < nrCromozomi; ++i){
            double u = (double)rand() / RAND_MAX; // numar aleator in [0, 1]
            int index = binarySearch(intervaleSelectie, u); // cautam intervalul in care se afla u
            if(index == nrCromozomi)
                index--; // daca u este 1, luam ultimul interval
            indiviziSelectati.push_back(populatie[index]); // selectam cromozomul
            }

    
        vector<int> indiviziIncrucisati;
        for (int i = 0; i < nrCromozomi; ++i) {
            double u = (double)rand() / RAND_MAX; // numar aleator in [0, 1]
            double x = decodare(indiviziSelectati[i]);
            if (u < probCrossover / 100.0) {
                indiviziIncrucisati.push_back(i); // adaugam individul la lista de indivizi incrucisati
            }
        }
        
        if(indiviziIncrucisati.size() % 2 != 0) {
            indiviziIncrucisati.pop_back(); // eliminam ultimul individ daca numarul de indivizi incrucisati este impar
        }


        for(int i = 0; i < indiviziIncrucisati.size(); i+=2) {

            int punctIncrucisare = rand() % lungimeCromozom; // punct de incrucisare
            Individ individ1 = indiviziSelectati[indiviziIncrucisati[i]];
            Individ individ2 = indiviziSelectati[indiviziIncrucisati[i+1]];
        
            for(int j = 0; j < punctIncrucisare; j++) {
                swap(individ1[j], individ2[j]); // facem swap la cromozom
            }

            indiviziSelectati[indiviziIncrucisati[i]] = individ1; // actualizam individul 1
            indiviziSelectati[indiviziIncrucisati[i+1]] = individ2; // actualizam individul 2
            }

        // folosim mutatia rara
        for(int i = 0; i < nrCromozomi; ++i) {
            double u = (double)rand() / RAND_MAX; // numar aleator in [0, 1]
            if( u < probMutatie / 100.0) {
                int pozitie = rand() % lungimeCromozom; // pozitia la care se face mutatia
                indiviziSelectati[i][pozitie] = 1 - indiviziSelectati[i][pozitie]; // facem mutatia prin inversarea bitului
            }
        }

        double maxFitness = fitnessFunction(decodare(indiviziSelectati[0]));
        double meanFitness = 0.0;
        for(int j = 1; j < nrCromozomi; ++j) {
            double x = decodare(indiviziSelectati[j]);
            double fitnessVal = fitnessFunction(x);
            meanFitness += fitnessVal;
            if(fitnessVal > maxFitness) {
                maxFitness = fitnessVal;
            }
        }
        meanFitness /= nrCromozomi;
        out << "Generatia " << i+1 << " max fitness: " << maxFitness << endl;
        
        populatie = indiviziSelectati; // populatia devine populatia selectata
        
    }
    return 0;
}


// -x^4 + 4x^2 + 2x + 4