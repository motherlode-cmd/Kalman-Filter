#pragma once
#include "../model/Reciever.hpp"

struct TMeasurement  {//вышка
    double TOA;//Time Of Arrival
    std::vector <double> receiver;
};

struct TTDOAMeasurement {//Разница между двумя вышками
    double TDOA;//Time d of arivel - разниа времен между 1 и 2
    std::vector <std::vector <double> > recievers;//double [][]
    //[x1][y1][z1]
    //[x2][y2][z2]
};
