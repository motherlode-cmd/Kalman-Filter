#pragma once
#include  "../const/physic.hpp"
#include  "../types/navigation.hpp"
#include "matrix.hpp"

template<typename T>
std::vector<T> slice(std::vector<T> v, int m, int n){
    std::vector <T> vec;
    for(int i = m; i < n; i++)
        vec.push_back(v[i]);
    return vec;
}

double TtoR( double t) {
    return t * LIGHT_SPEED;
}

class CoordsCalc {
public:
    CoordsCalc() = default;

    std::vector <std::vector<my_Vector>> getClosedSolution(std::vector <Receiver> mes) {
        if (mes.size() != 4) {
            std::vector <std::vector<my_Vector>> results;
            for (int i = 0; i < (mes.size() / 4) ; i++) {
                std::cout<<i<<std::endl;
                std::vector <std::vector<my_Vector>> solve = getClosedSolution(slice(mes, i * 4, (i + 1) * 4));
                for(auto i : solve)
                    results.push_back(i);
            }
            return results;
        } else {
            std::cout<<"calc in CoordCalc::getClosedSolution"<<std::endl;
            // calculations
            Matrix A = Matrix(4, 4);
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 3; j++) {
                    A.set(i, j, mes[i].pos.vec[j]);
                }
                A.set(i, 3, -TtoR(mes[i].signals.empty() ? 0 : mes[i].signals[0]));
            }//Матрица Адля решения
            A.print();
            std::vector <my_Vector> s_i(4);//Вектора для вышек
            for(int i = 0; i < 4; i++) {
                s_i[i].fill(1, 4);
                for(int j = 0; j < 4; j++) {
                    s_i[i].set(j, A.matrix[i][j]);
                }
            }
            //const b = new Vec(...s_i.map((s) => s.lorentzianProduct(s)));
            my_Vector b(4);
            for (int i = 0; i < 4; i++) {
                b.set(i, s_i[i].lorentz_product(s_i[i]));
            }
            Matrix invA = A.invert();
            my_Vector tmp;//[1 1 1 1]
            tmp.fill(1, 4);
            my_Vector d = 0.5 * invA * tmp;
            my_Vector e = 0.5 * invA * b;
            double alpha = d.lorentz_product(d);
            double betta = 2 * d.lorentz_product( e) - 1;
            double gamma = e.lorentz_product(e);//коэф квадрат уравнения
            //рещение кв уравнения
            double lambda1 = (-betta + sqrt(betta * betta - 4 * alpha * gamma)) / (2 * alpha);
            double lambda2 = (-betta - sqrt(betta * betta - 4 * alpha * gamma)) / (2 * alpha);
            std::vector <std::vector <my_Vector> > results;
            std::vector <my_Vector> sub_result;
            sub_result.push_back(lambda1 * d + e);
            sub_result.push_back(lambda2 * d + e);
            results.push_back(sub_result);
            std::cout<<"Closed solution "<<std::endl;
            sub_result[0].print();
            sub_result[1].print();
            return results;
        }
    }

    double rCalc(std::vector <double> coords, std::vector <double> reciever) {
        double d0 = coords[0] - reciever[0];
        double d1 = coords[1] - reciever[1];
        double d2 = coords[2] - reciever[2];
        return sqrt(d0 * d0 + d1 * d1 + d2 * d2);
    }

    double dmCalc(std::vector <double> coords, std::vector <std::vector <double>> receivers) {
        return rCalc(coords, receivers[0]) - rCalc(coords, receivers[1]);
    }

    my_Vector getIterativeSolutionByMLS(std::vector <TTDOAMeasurement> measurements,my_Vector approx) {
        Matrix A(measurements.size(), 3);

        std::vector <double> coords = approx.vec;

        for(int i = 0; i < A.matrix.size(); i++) {
            for(int j = 0; j < A.matrix[i].size(); j++) {
                std::vector <std::vector<double>> receivers = measurements[i].recievers;
                double r_i = rCalc(coords, receivers[0]);
                double r_j = rCalc(coords, receivers[1]);
                A.set(i, j, (receivers[1][j] - coords[j]) / r_j - (receivers[0][j] - coords[j]) / r_i);
            }
        }
        // ? subtract the real measurements by the measurements derived from the approximate solution m(x_a)
        my_Vector dm;
        for (int i = 0; i < measurements.size(); i++) {
            double el = TtoR(measurements[i].TDOA) - dmCalc(coords, measurements[i].recievers);
            dm.vec.push_back(el);
        }
        //std::cout<<"dm size "<<dm.vec.size()<<std::endl;
        //dm.print();
        Matrix A_T = A.transpose();
        my_Vector BVec = A_T * dm; // b
        Matrix MatrixDecomposition = (A_T * A); // L^TL = A
        my_Vector solution = MatrixDecomposition.solve_system(BVec);
        my_Vector AdjustmentSolution = approx + solution;
        return AdjustmentSolution;
    }

    bool is_one(std::vector <my_Vector> & vecArr, my_Vector & vec1) {
        for(auto vec : vecArr)
            if((vec1 - vec).len()<= E)
                return true;
        return false;
    }
    bool for_all(std::vector <std::vector <my_Vector>> & from, my_Vector & pattern_vec) {
        for(auto vecArr : from) {
            if(!is_one(vecArr, pattern_vec)) return false;
        }
        return true;
    }
    my_Vector filterResult(std::vector <std::vector <my_Vector>> from) {
        std::vector <my_Vector> pattern = from[0];
        for(auto el : pattern) {
            if(for_all(from, el)) return el;
        }
        return pattern[0];
    }
};
