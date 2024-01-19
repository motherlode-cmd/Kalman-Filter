#pragma once
#include  "Navigation.hpp"
#include "vector.hpp"
#include "future"
#include <chrono>
#include <thread>
#include <fstream>
int counter = 0;
void detection(bool stop, Navigation & CommandCenter, std::ofstream & real,
               std::ofstream & mls_result, std::ofstream & kalman_result, std::ofstream & error) {
    if (stop) CommandCenter.stopAircraft();
    CommandCenter.makeCheck();
    my_Vector MLS = CommandCenter.findCoords();
    std::cout<<"Real Pos :"<<std::endl;
    CommandCenter.print_Air();
    CommandCenter.get_air_pos().print(real);
    std::cout<<"MLS res: "<<std::endl;
    MLS.print();
    MLS.print(mls_result);
    error<<counter<<" "<<CommandCenter.getPositionMod(MLS) / CommandCenter.get_air_pos().len()<<" ";
    auto filterRes = CommandCenter.getLastFilterResult();
    //std::cout<<"HERE "<<filterRes.vec.size()<<std::endl;
    //filterRes.print();
    std::cout<<"Filter res : " <<std::endl;
    filterRes.cut(3).print();
    filterRes.cut(3).print(kalman_result);
    error<<CommandCenter.getPositionMod(filterRes) / CommandCenter.get_air_pos().len() << std::endl;
    counter++;
    CommandCenter.startAircraft();
}
