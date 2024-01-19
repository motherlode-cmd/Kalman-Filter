#pragma once
#include "../const/time.hpp"
#include "Reciever.hpp"
#include "Aircraft.hpp"
#include "CoordsCalc.hpp"
#include "Filter.hpp"
#include  "../const/kalmanFilter.hpp"
#include  "matrix.hpp"
#include "chrono"
#include "../types/navigation.hpp"
#include <fstream>
//time
double getNow() {
    std::chrono::milliseconds currentTime = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );
    int64_t currentTimeInMilliseconds = currentTime.count();
    std::ifstream ifs;
    ifs.open("time.txt");
    double time_start;
    ifs>>time_start;
    //std::cout<<time_start;
   // std::cout<<(currentTimeInMilliseconds / 1000 - time_start)<<std::endl;
    return (currentTimeInMilliseconds / 1000 - time_start);
}

//navigation

class Navigation {
private:
    std::vector <Receiver> _receivers;
    Aircraft _aircraft;
    std::vector <my_Vector> pathHistory;
    KalmanFilter Filter;
    double lastCheckTime = 0;
    double deltaTime = 0;
    double startTime = 0;
public:
    Navigation() {
        TInputData settings;
        settings.F = getFMatrix(0, Matrix(9, 9));
        settings.H = getHMatrix();
        settings.Q = getQMatrix(0, Sigma, Matrix(9, 9));
        settings.R = getRMatrix(Sigma);
        settings.P = getPMatrix(Sigma);
        settings.x = my_Vector(9);
        settings.x.fill(0, 9);
        KalmanFilter kal(settings);
        this->Filter = kal;
    }

    my_Vector getLastFilterResult() {
        return this->Filter.get_state();
    }

    Receiver createReceiver(my_Vector pos) {
        _receivers.push_back(Receiver(pos));
        //std::cout<<"Created rec size = "<<_receivers.size()<<std::endl;
        pos.print();
        return this->_receivers[this->_receivers.size() - 1];
    }

    Aircraft createAircraft(my_Vector pos) {
        this->_aircraft = Aircraft(pos);
        //std::cout<<"Created air  "<<std::endl;
        pos.print();
        return this->_aircraft;
    }

    void makeCheck() {
        double dateNow = getNow();
        this->deltaTime = dateNow - this->startTime - this->lastCheckTime;
        this->lastCheckTime = dateNow - this->startTime;
        //std::cout<<"DELTA TIME "<< this->deltaTime<<std::endl;
        this->_aircraft.moveByDt(this->deltaTime);
        //? Model limitations. We do not simulate the flight of light and get just delay
        std::vector <double> delays;
        for (int i = 0; i < this->_receivers.size(); i++) {
            Receiver bcn = this->_receivers[i];
            my_Vector pos = bcn.pos;
            //check aircraft
            double delay = this->_aircraft.getLightDelay(my_Vector(pos));
            delays.push_back(delay);
        }

        //? Just as if the package from the receivers was sent to the aircraft, and all the clocks on the receivers are synchronized
        // ? ms/1000 = s
        for (int i = 0; i < _receivers.size(); i++) {
            double del = delays[i];
            double signalTime = dateNow + del;
            //std::cout<<"Signal "<<del<<std::endl;
            _receivers[i].acceptSignal(signalTime);
        }
    }

    // Generate Array<TDOA> from Array<TOA>
    std::vector <TTDOAMeasurement> getTDOA (std::vector <TMeasurement> & measurements) {
        std::vector <TTDOAMeasurement> result;
        for (int i = 0; i < measurements.size() - 1; i++) {
            for (int j = i + 1; j < measurements.size(); j++) {
                //? E(noise) == 0
                //? D(noie) == 75
                double noise = (random() % 10 - 0.5) / 80000000;
                TTDOAMeasurement tmp;
                tmp.TDOA = measurements[i].TOA -  measurements[j].TOA + noise;
                tmp.recievers.push_back(measurements[i].receiver);
                tmp.recievers.push_back(measurements[j].receiver);
                result.push_back(tmp);
            }
        }
        return result;
    }

    my_Vector findCoords() {
        CoordsCalc coords;
        my_Vector firstResult;
        //std::cout<<"Here find coords::Navigation"<<std::endl;

        if (this->pathHistory.size() == 0) {
            //std::cout<< this->_receivers.size()<<std::endl;
            firstResult = coords.filterResult(coords.getClosedSolution(this->_receivers));
        }
        my_Vector approx = !(this->pathHistory.empty()) ? this->pathHistory[this->pathHistory.size() - 1] : firstResult;
        my_Vector MLSResult;
        //std::cout<<"Approx "<<std::endl;
        //approx.print();
        for (int i = 0; i < TIMES_TO_MLS; i++) {
            std::vector <TMeasurement> measurements;
            for(auto el : this->_receivers) {
                TMeasurement meas;
                meas.TOA = el.signals.empty() ? 0 : el.signals[0];
                meas.receiver = (el.pos.cut(3)).vec;
                measurements.push_back(meas);
            }
            MLSResult = coords.getIterativeSolutionByMLS(
                    this->getTDOA(measurements),
                    approx.cut(3)
            );
            approx = MLSResult;
        }

        this->pushHistory(MLSResult);
        //std::cout<<"Got History"<<std::endl;
        //MLSResult.print();
        return MLSResult;
    }

    my_Vector filter(my_Vector z) {
        //std::cout<<"Delta time is "<<this->deltaTime<<std::endl;
        this->Filter.set_FMatrix(getFMatrix(this->deltaTime, this->Filter.get_FMatrix()));
        this->Filter.set_QMatrix(getQMatrix(this->deltaTime, Sigma, this->Filter.get_QMatrix()));

        my_Vector result = this->Filter.update(z);
        return result;
    }

    void startAircraft() {
        this->startTime = getNow();
    }

    void print_Air(){
        this->_aircraft.print_pos();
    }

    my_Vector get_air_pos(){
        return _aircraft.get_pos();
    }

    void stopAircraft() {
        this->startTime = -1;
        this->lastCheckTime = 0;
        this->deltaTime = 0;
    }

    void pushHistory(my_Vector point) {
        my_Vector result = point.cut(3);
        this->pathHistory.push_back(result);
        this->filter(result);
    }

    double getPositionMod(my_Vector vec) {
        return this->_aircraft.getPositionMod(vec);
    }

    std::vector <Receiver> get_receivers() {
        return this->_receivers;
    }
};
