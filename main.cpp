#include "tsp.cpp"
#include<chrono>
#include<iomanip>

int main() {

    TSP::Parser parsee("resources/80rd400.gtsp");

    TSP::Graph graph = parsee.run();

    cout << "\n" << "Simulation " << parsee.name << "\n"; 
    
    TSP::Genetic simulation(graph);

    const auto time1 = chrono::steady_clock::now();
    
    simulation.run();
    
    const auto time2 = chrono::steady_clock::now();

    cout << "\n" << "Time: " << chrono::duration_cast<chrono::microseconds>(time2 - time1).count() << "\n";
    
    
    return 0;
}