#pragma once
#include "math.h"
#include <iostream>
#include <vector>
struct my_Vector {
    std::vector <double> vec;
    my_Vector() = default;
    ~my_Vector() {vec.clear();}
    my_Vector(int size) {vec = std::vector <double>(size);}
    my_Vector(const my_Vector & copyable) {
        for(auto i : copyable.vec) {
            vec.push_back(i);
        }
    }
    my_Vector(my_Vector && copyable) {
        for(auto i : copyable.vec) {
            vec.push_back(i);
        }
        copyable.vec.clear();
    }
    my_Vector & operator = (const my_Vector & copyable) {
        vec.clear();
        for(auto i : copyable.vec) {
            vec.push_back(i);
        }
        return *this;
    }
    my_Vector & operator = (my_Vector && copyable)  noexcept {
        vec.clear();
        for(auto i : copyable.vec) {
            vec.push_back(i);
        }
        copyable.vec.clear();
        return *this;
    }

    my_Vector reverse() {
        my_Vector new_vec(*this);
        for(int i = 0; i < vec.size(); i++) {
            new_vec.vec[i] = vec[vec.size() - 1 - i];
        }
        return new_vec;
    }

    friend my_Vector operator + (const my_Vector & v1,const my_Vector & v2) {
        my_Vector add;
        for(int i= 0; i < v1.vec.size(); i++) {
            add.vec.push_back(v1.vec[i] + v2.vec[i]);
        }
        return add;
    }

    friend my_Vector operator - (const my_Vector & v1,const my_Vector & v2) {
        my_Vector add;
        for(int i= 0; i < v1.vec.size(); i++) {
            add.vec.push_back(v1.vec[i] - v2.vec[i]);
        }
        return add;
    }

    void fill(double a, int count) {
        vec.clear();
        for(int i = 0; i < count; i++) {
            vec.push_back(a);
        }
    }

    friend my_Vector operator * (const my_Vector & v1, double k) {
        my_Vector mul;
        for(auto i : v1.vec) {
            mul.vec.push_back(i * k);
        }
        return mul;
    }
    friend my_Vector operator * (double k, const my_Vector & v1) {
        my_Vector mul;
        for(auto i : v1.vec) {
            mul.vec.push_back(i * k);
        }
        return mul;
    }
    double lorentz_product(my_Vector v2) const {
        if(vec.size() != v2.vec.size()) return 0;
        double result = 0;
        for(int i = 0; i < vec.size()- 1; i++) {
            result += vec[i] * v2.vec[i];
        }
        result -= vec[vec.size() - 1] * v2.vec[v2.vec.size() - 1];
        return result;
    }

    double len() {
        double result = 0;
        for(auto i : vec) {
            result += i * i;
        }
        return sqrt(result);
    }
    my_Vector cut(int count) {
        std::vector <double> copy(vec);
        vec.clear();
        for(int i = 0; i < count; i++) {
            vec.push_back(copy[i]);
        }
        copy.clear();
        return *this;
    }
    std::vector <double> get_vec() {
        return vec;
    }
    void set(int i, double val) {
        vec[i] = val;
    }
    void add_coord(int i, double val) {
        vec[i] += val;
    }
    void print() {
        for(auto i : vec) {
            std::cout<<" "<<i;
        }
        std::cout<<std::endl;
    }
    void print(std::ostream & os) {
        for(auto i : vec) {
            os<<" "<<i;
        }
        os<<std::endl;
    }
};
