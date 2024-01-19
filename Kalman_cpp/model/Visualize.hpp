#pragma once
#include "vector.hpp"

class Visualize {
private:
    std::vector <std::vector<my_Vector>> traces;
    std::string id;
public:
    Visualize(std::string id) {
        this->id = id;
    }
    Visualize init() {
        std::cout<<id<<std::endl;
        return *this;
    }
    void addTrace(std::vector<my_Vector> config) {
        this->traces.push_back(config);
    }

    void extendsTraceByVec(my_Vector vec) {
        this->traces[this->traces.size() - 1].push_back(vec);
    }
    void printTraces() {

    }
};