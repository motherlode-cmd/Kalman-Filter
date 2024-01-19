#pragma once
#include "vector.hpp"
#include "vector"
/**
 * @description Receiver class is used to simulate the receiver.
 * It has a position in 3D space and can accept signals from the aircraft.
 * It controlled by Navigation class, like the aircraft, like a receiver in real life.
 */
class Receiver {
public:
    std::vector <double> signals;
    std::string id;
    my_Vector pos;
    Receiver(my_Vector pos) {
        this->pos = pos;
    }
    my_Vector get_pos() {
        return this->pos;
    }

    /** @description accept signal form the aircraft (navigation post) you have to divide difference by TIMEOUT_FACTOR before using */
    void acceptSignal(double comingTime) {
        this->signals.clear();
        this->signals.push_back(comingTime);
        //std::cout<<signals.size()<<" "<<comingTime<<"sig catch"<<std::endl;
    }

    void clear() {
        this->signals.clear();
    }
    std::vector <double> get_signals() {
        return this->signals;
    }
};
