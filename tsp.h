#ifndef H_TSP
#define H_TSP

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<random>
#include<cmath>
#include<algorithm>


using namespace std;

random_device rd;
mt19937 rng(rd());

namespace TSP {
    
    class Graph {
        
        public:
        string name; // name of problem
        int nDim; // num of nodes
        int gDim; // num of 
        vector<vector<double>> weightMatrix;
        vector<vector<int>> nodeGroup;
    
        Graph() = default;
        Graph(int nDim, int gDim, vector<vector<double>> weightMatrix, vector<vector<int>> nodeGroup);
    
    };
    
    class Parser {

        public:
        vector<pair<double,double>> coordMatrix;
        vector<vector<double>> weightMatrix;
        vector<vector<int>> nodeGroup;
        ifstream parsee;
        string name;
        string type = "GTSP";
        int nDim = 1; // num of nodes
        int gDim = 1; // num of groups
        int capacity;
        string edgeWeightType;
        string edgeWeightFormat;
        string edgeDataFormat;
        string nodeCoordType = "NO_COORDS";
        string displayDataType = "NO_DISPLAY";
        
        void parseNodeCoord();
        void parseEdgeWeight();
        void parseNodeGroup();
       
        Parser(string fileName);
        string parseInfo();
        Graph run();

    };

    class Genetic {
        
        public:
        using Chromosome = vector<pair<int,int>>;
        Graph graph;
        int size;
        vector<Chromosome> population;
        double selectionRate;
        double crossoverRate;
        double mutationRate;
        int repeatingStopLimit; // times for when top cost appear to stop
        int longStopLimit; // num of max generations for simulation
       
        double calculateEdgeCost(pair<int,int> geneA, pair<int,int> geneB);
        double calculateTotalCost(Chromosome chromo);
        void sort();
        void select(); // remove worst-fits
        void crossover(vector<Chromosome>& newPopulation, Chromosome chromoA, Chromosome chromoB); // greedy select crossover GSX
        void twoOptMutation(vector<Chromosome>& newPopulation, Chromosome chromo); // 2opt only for consec pairs
        void nnsMutation(vector<Chromosome>& newPopulation, Chromosome chromo); // nearest neighbor search on a 
        void create(vector<Chromosome>& newPopulation); // randomize new chromosomes
        void newPopulation(); // new pop for new generation

        Genetic() = default;
        Genetic(Graph& graph, int size, double selectionRate, double crossoverRate, double mutationRate, int repeatingStopLimit, int longStopLimit);
        void run();

    };
}

#endif