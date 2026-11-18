#include "tsp.h"

TSP::Parser::Parser(string fileName) {
    parsee = ifstream(fileName);
}

string TSP::Parser::parseInfo() {
    string parsit = "";
    size_t it;

    while (parsit == "") {
        parsee >> parsit;
        while ((it = parsit.find(":")) != string::npos) {
            parsit.erase(it);
        };
        if (parsit != "") {
            return parsit;
        };
    }
}

TSP::Graph TSP::Parser::run() {
    string parsit;
    
    while (parsee.tellg() != parsee.end) {
        parsit = parseInfo();
        if (parsit == "NAME") { // info section
            name = parseInfo();
        } else if (parsit == "TYPE") {
            type = parseInfo();
        } else if (parsit == "DIMENSION") {
            nDim = stoi(parseInfo());
        } else if (parsit == "GTSP_SETS") {
            gDim = stoi(parseInfo());
        } else if (parsit == "CAPACITY") {
            capacity = stoi(parseInfo());
        } else if (parsit == "EDGE_WEIGHT_TYPE") {
            edgeWeightType = parseInfo();
        } else if (parsit == "EDGE_WEIGHT_FORMAT") {
            edgeWeightFormat = parseInfo();
        } else if (parsit == "EDGE_DATA_FORMAT") {
            edgeDataFormat = parseInfo();
        } else if (parsit == "NODE_COORD_TYPE") {
            nodeCoordType = parseInfo();
        } else if (parsit == "EDGE_WEIGHT_FORMAT") {
            edgeWeightFormat = parseInfo();
        } else if (parsit == "DISPLAY_DATA_TYPE") {
            displayDataType = parseInfo();
        } else if (parsit == "NODE_COORD_SECTION") { // data section
            TSP::Parser::parseNodeCoord();
        } else if (parsit == "EDGE_WEIGHT_SECTION") {
            TSP::Parser::parseEdgeWeight();
        } else if (parsit == "GTSP_SET_SECTION") {
            TSP::Parser::parseNodeGroup();
        } else if (parsit == "EOF") {
            break;
        }
    }
    
    return TSP::Graph::Graph(nDim, gDim, weightMatrix, nodeGroup);
};

void TSP::Parser::parseNodeCoord() {
    string parsit;

    if (edgeWeightType == "EUC_2D") {
        for (int i = 0; i < nDim; i++) {
            parsee >> parsit;
            coordMatrix.push_back(pair<double,double>());
            parsee >> parsit;
            coordMatrix[i].first = stod(parsit);
            parsee >> parsit;
            coordMatrix[i].second = stod(parsit);
        }
        weightMatrix = vector<vector<double>>(nDim, vector<double>(nDim));
        double cost;
        for (int i = 0; i < nDim ; i++) {
            for (int j = 0; j <= i; j++) {
                cost = (sqrt(pow(coordMatrix[i].first - coordMatrix[j].first, 2) + pow(coordMatrix[i].second - coordMatrix[j].second, 2)));
                weightMatrix[i][j] = cost;
                weightMatrix[j][i] = cost;
            }
        }
    }
    
}

void TSP::Parser::parseEdgeWeight() {
    string parsit;
    
    if (edgeWeightFormat == "FULL_MATRIX") {
        for (int i = 0; i < nDim; i++) {
            weightMatrix.push_back(vector<double>());
            for (int j = 0; j < nDim; j++) {
                parsee >> parsit;
                weightMatrix[i].push_back(stod(parsit));
            }
        }
    } else return;
}

void TSP::Parser::parseNodeGroup() {
    string parsit;

    for (int i = 0; i < gDim; i++) {
        parsee >> parsit;
        nodeGroup.push_back(vector<int>());
        while (parsit != "-1") {
            parsee >> parsit;
            if (parsit != "-1") {
                nodeGroup[i].push_back(stoi(parsit) - 1);
            }
        };
    }
}

TSP::Graph::Graph(int nDim, int gDim, vector<vector<double>> weightMatrix, vector<vector<int>> nodeGroup) {
    this->nDim = nDim;
    this->gDim = gDim;
    this->weightMatrix = weightMatrix;
    this->nodeGroup = nodeGroup;
}

TSP::Genetic::Genetic(TSP::Graph& graph, int size = 100, double selectionRate = 0.95, double crossoverRate = 0.50, double mutationRate = 0.10, int repeatingStopLimit = 250, int longStopLimit = 2500) {
    this->graph = graph;
    this->size = size;
    this->selectionRate = selectionRate;
    this->crossoverRate = crossoverRate;
    this->mutationRate = mutationRate;
    this->repeatingStopLimit = repeatingStopLimit;
    this->longStopLimit = longStopLimit;
}

double TSP::Genetic::calculateEdgeCost(pair<int,int> geneA, pair<int,int> geneB) {
    return graph.weightMatrix[graph.nodeGroup[geneA.first][geneA.second]][graph.nodeGroup[geneB.first][geneB.second]];
}

double TSP::Genetic::calculateTotalCost(Chromosome chromo) {
    double cost;
    for (int i = 0; i < graph.gDim; i++) { // loop around the chromosome
        cost += calculateEdgeCost(chromo[i], chromo[(i + 1) % graph.gDim]);
    }
    return cost;
}

void TSP::Genetic::sort() {
    std::sort(population.begin(), population.end(), [this](const Chromosome& chromoA, const Chromosome& chromoB) {
        return this->calculateTotalCost(chromoA) < this->calculateTotalCost(chromoB);
    });
}

void TSP::Genetic::select() {
    population.erase(population.begin() + round(size * selectionRate), population.end());
}

void TSP::Genetic::crossover(vector<Chromosome>& newPopulation, Chromosome chromoA, Chromosome chromoB) {
    vector<Chromosome> parents = {chromoA, chromoB};
    Chromosome genePool, newChromo;

    vector<int> tempGroup;
    
    uniform_int_distribution<int> distr(0, graph.gDim - 1);
    
    int greedy = distr(rng);
    
    for (auto chromo : parents) {
        for (int i = 0; i < graph.gDim; i++) {
            if (chromo[i].first == greedy) {
                newChromo.push_back(chromo[i]);
                tempGroup.push_back(greedy);
                break;
            }
        }
        while (newChromo.size() < graph.gDim) {
            for (auto parent : parents) {
                for (int i = 0; i < graph.gDim; i++) {
                    if (parent[i].first == greedy) {
                        genePool.push_back(parent[(i+graph.gDim-1)%graph.gDim]);
                        genePool.push_back(parent[(i+1)%graph.gDim]);
                        break;
                    }
                }
            }
            std::sort(genePool.begin(), genePool.end(), [this, newChromo](const pair<int,int> geneA, const pair<int,int> geneB) {
                return this->calculateEdgeCost(newChromo[newChromo.size() - 1], geneA) < this->calculateEdgeCost(newChromo[newChromo.size() - 1], geneB);
            });
            while (find(tempGroup.begin(), tempGroup.end(), genePool[0].first) != tempGroup.end()) {
                genePool.erase(genePool.begin());
                if (genePool.empty() != false) break;
            }
            if (genePool.empty() != true) {
                newChromo.push_back(genePool[0]);
                greedy = genePool[0].first;
                tempGroup.push_back(greedy);
            } else {
                do {
                    greedy = distr(rng);
                } while (find(tempGroup.begin(), tempGroup.end(), greedy) != tempGroup.end());
                for (auto parent : parents) {
                    for (int i = 0; i < graph.gDim; i++) {
                        if (parent[i].first == greedy) {
                            genePool.push_back(parent[i]);
                            break;
                        }
                    }
                }
                std::sort(genePool.begin(), genePool.end(), [this, newChromo](const pair<int,int> geneA, const pair<int,int> geneB) {
                    return this->calculateEdgeCost(newChromo[newChromo.size() - 1], geneA) < this->calculateEdgeCost(newChromo[newChromo.size() - 1], geneB);
                });
                newChromo.push_back(genePool[0]);
                tempGroup.push_back(greedy);
            }
            genePool.clear();
        }
        newPopulation.push_back(newChromo);
        newChromo.clear();
        tempGroup.clear();
    }
}

void TSP::Genetic::twoOptMutation(vector<Chromosome>& newPopulation, Chromosome chromo) {
    for (int i = 0; i < graph.gDim; i++) {
        double delta =  - calculateEdgeCost(chromo[i], chromo[(i + 1) % graph.gDim])
                        - calculateEdgeCost(chromo[(i + 2) % graph.gDim], chromo[(i + 3) % graph.gDim])
                        + calculateEdgeCost(chromo[i], chromo[(i + 2) % graph.gDim])
                        + calculateEdgeCost(chromo[(i + 1) % graph.gDim], chromo[(i + 3) % graph.gDim]);

        if (delta < 0) {
            swap(chromo[(i+1)%(graph.gDim)], chromo[(i+2)%(graph.gDim)]);
        }
    }
    newPopulation.push_back(chromo);
}

void TSP::Genetic::nnsMutation(vector<Chromosome>& newPopulation, Chromosome chromo) {
    uniform_int_distribution<int> distr(0, graph.gDim - 1);
    int a = distr(rng), b = distr(rng);

    for (int i = 0; i < graph.gDim; i++) {
        if (i < min(a, b) || i > max(a, b)) continue;
        else {
            double nearestCost = INFINITY;
            int neighbor;
            for (int j = i + 1; j <= max(a, b); j++) {
                if (calculateEdgeCost(chromo[i], chromo[j]) < nearestCost) {
                    nearestCost = calculateEdgeCost(chromo[i], chromo[j]);
                    neighbor = j;
                }
            }
            swap(chromo.begin() + i + 1, chromo.begin() + neighbor);
        }
    }
    newPopulation.push_back(chromo);
}

void TSP::Genetic::create(vector<Chromosome>& newPopulation) {
    Chromosome newChromo;
    vector<int> tempGroup(graph.gDim);
    int tempNode;
    for (int i = 0; i < graph.gDim; i++) tempGroup[i] = i;
    while (newPopulation.size() < size) {
        shuffle(tempGroup.begin(), tempGroup.end(), rng);
        for (auto j : tempGroup) {
            uniform_int_distribution<int> distr(0, (int) graph.nodeGroup[j].size() - 1);
            tempNode = distr(rng);
            newChromo.push_back({j, tempNode});
        }
        newPopulation.push_back(newChromo);
        newChromo.clear();
    }
}

void TSP::Genetic::newPopulation() {
    vector<Chromosome> newPopulation, crossoverPool;
    discrete_distribution<int> distr({crossoverRate, mutationRate / 2, mutationRate / 2, 1 - crossoverRate - mutationRate});
    int chromoFate;

    select();
    newPopulation.push_back(population[0]);
    for (auto chromo : population) {
        chromoFate = distr(rng);
        if (chromoFate == 0) {
            crossoverPool.push_back(chromo);
        } else if (chromoFate == 1) {
            twoOptMutation(newPopulation, chromo);
        } else if (chromoFate == 2) {
            nnsMutation(newPopulation, chromo);
        } else if (chromoFate == 3) {
            newPopulation.push_back(chromo);
        }
    }
    while (crossoverPool.size() >= 2) {
        crossover(newPopulation, crossoverPool[0], crossoverPool[1]);
        crossoverPool.erase(crossoverPool.begin(), crossoverPool.begin() + 2);
    }
    if (crossoverPool.size() == 1) newPopulation.push_back(crossoverPool[0]);
    create(newPopulation);
    population = newPopulation;
}

void TSP::Genetic::run() {
    int repeatingCount = 0, longCount = 0;
    double bestCost = INFINITY;
    
    create(population);

    while (repeatingCount < repeatingStopLimit && longCount < longStopLimit) {
        longCount++;
        sort();
        if (calculateTotalCost(population[0]) < bestCost) {
            bestCost = calculateTotalCost(population[0]);
            repeatingCount = 1;
        } else repeatingCount++;
        newPopulation();
    }

    cout    << "\n"
            << "Number of Generations: " << longCount << "\n"
            << "Best obtained cost: " << bestCost << "\n";
}