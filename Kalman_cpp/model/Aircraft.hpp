#pragma once
#include "../const/physic.hpp"
#include  "vector.hpp"

class Aircraft {
    std::string id;
    my_Vector pos;
    my_Vector velocity;
public:
    Aircraft(){
        pos = my_Vector(3);
        pos.fill(0,3);
        velocity = my_Vector(3);
    }
    Aircraft(my_Vector pos, my_Vector velocity) {
        this->pos = pos;
        this->velocity = velocity;
    }
    Aircraft(my_Vector pos) {
        this->pos = pos;
        velocity = my_Vector(3);
        velocity.fill(10, 3);
    }

    void moveByDt(double dt) {
        this->pos = this->pos + (this->velocity * dt);
        my_Vector dv(3);
        dv.set(0, 0.5); dv.set(1, 0.1); dv.set(2, 1);
        this->velocity = this->velocity + (dv * dt);
    }
    /** @description get the time in seconds it takes light to reach the aircraft */
    double getLightDelay(my_Vector fromPos) {
        return (this->pos - fromPos.cut(3)).len() / LIGHT_SPEED;
    }

    double getPositionMod(my_Vector fromPos) {
        my_Vector slice_from(3);
        for(int i = 0; i < 3; i++)
            slice_from.set(i, pos.vec[i] - fromPos.vec[i]);
        return slice_from.len();
    }
    my_Vector get_pos(){
        return pos;
    }
    void print_pos(){
        pos.print();
    }
};
