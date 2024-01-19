//
// Created by motherlode on 19.01.23.
//
#ifndef UNTITLED4_SYSTEM_H
#define UNTITLED4_SYSTEM_H
#include <iostream>
using namespace std;
typedef struct Vector {
    double x = 0;
    double y = 0;
    double z = 0;
    Vector() = default;
    Vector(double tx, double ty, double tz):x(tx), y(ty), z(tz){}
    Vector(Vector const & newX):x(newX.x), y(newX.y), z(newX.z){}
    Vector & operator = (const Vector &newX){
        x = newX.x;
        y = newX.y;
        z = newX.z;
        return *this;
    }
    friend ostream & operator << (ostream& os, const Vector& vec){
        os<< vec.x<<" "<<vec.y<<" "<<vec.z<<' ';
        return os;
    }
    friend Vector operator + (const Vector& x1, const Vector& x2){
        return Vector(x1.x + x2.x, x1.y + x2.y, x1.z + x2.z );
    }
    friend Vector operator * (const Vector& x1, double k){
        return Vector(x1.x * k, x1.y * k, x1.z * k );
    }
    friend double operator * (const Vector& x1, const Vector& x2){
        return (x1.x * x2.x + x1.y * x2.y + x1.z * x2.z);
    }
    friend Vector operator / (const Vector& x1, double k) {
        if(k != 0)
            return Vector(x1.x / k, x1.y / k, x1.z / k );
        else
            return Vector(x1);
    }
    double length() {
        return sqrt(x*x + y*y + z*z);
    }
} Vector;

Vector vec_mult(const Vector & v1,const Vector & v2) {
    /*
     * i j k
     * 1 2 3
       4 5 6 */
    Vector v(v1.y * v2.z - v1.z * v2.y, -(v1.x * v2.z - v1.z * v2.x), v1.x * v2.y - v1.y * v2.x);
    return v;
}

typedef struct Matrix {
    double  matrix[3][3]= {{1,1,1},
                           {0,1,0},
                           {1,0,1}};
    Matrix() = default;
    Matrix(Matrix const & newX){
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                matrix[i][j] = newX.matrix[i][j];
            }
        }
    }
    friend ostream & operator << (ostream& os, const Matrix& x){
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                os<<x.matrix[i][j]<<' ';
            }
        }
        return os;
    }
    Matrix & operator = (const Matrix &newX){
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                matrix[i][j] = newX.matrix[i][j];
            }
        }
        return *this;
    }
    friend Matrix operator + (const Matrix& x1, const Matrix& x2){
        Matrix matr;
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                matr.matrix[i][j] = x1.matrix[i][j] + x2.matrix[i][j];
            }
        }
        return matr;
    }
    friend Matrix operator * (const Matrix& x1, double k){
        Matrix matr;
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                matr.matrix[i][j] = x1.matrix[i][j] * k;
            }
        }
        return matr;
    }
    friend Vector operator * (const Matrix& x1, const Vector& v1){
        Vector vec;
        vec.x = x1.matrix[0][0] * v1.x + x1.matrix[0][1] * v1.y + x1.matrix[0][2] * v1.z;
        vec.y = x1.matrix[1][0] * v1.x + x1.matrix[1][1] * v1.y + x1.matrix[1][2] * v1.z;
        vec.z = x1.matrix[2][0] * v1.x + x1.matrix[2][1] * v1.y + x1.matrix[2][2] * v1.z;
        return vec;
    }
    friend Vector operator * (const Vector & v1,const Matrix& x1){
        Vector vec;
        vec.x = x1.matrix[0][0] * v1.x + x1.matrix[1][0] * v1.y + x1.matrix[2][0] * v1.z;
        vec.y = x1.matrix[0][1] * v1.x + x1.matrix[1][1] * v1.y + x1.matrix[2][1] * v1.z;
        vec.z = x1.matrix[0][2] * v1.x + x1.matrix[1][2] * v1.y + x1.matrix[2][2] * v1.z;
        return vec;
    }
    friend Matrix operator * (const Matrix& x1, const Matrix &x2){
        Matrix matr;
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                matr.matrix[i][j] = 0;
                for(int k = 0; k < 3; k++) {
                    matr.matrix[i][j] += x1.matrix[i][k] * x2.matrix[k][j];
                }
            }
        }
        return matr;
    }
    double det() const{
        /*
         * 1 2 3
         * 4 5 6
         * 7 8 9
         * */
        return matrix[0][0] * matrix[1][1] * matrix[2][2] +
                   matrix[1][0] * matrix[2][1] * matrix[0][2] +
                   matrix[0][1] * matrix[1][2] * matrix[2][0] -
                   matrix[0][2] * matrix[1][1] * matrix[2][0] -
                   matrix[0][1] * matrix[1][0] * matrix[2][2] -
                   matrix[2][1] * matrix[1][2] * matrix[0][0];
    }
    double minor(int n, int m) const{
        double a[2][2];
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                if(i != n && j != m) {
                    int ci = i > n ? i - 1 : i;
                    int cj = j > m ? j - 1 : j;
                    a[ci][cj] = matrix[i][j];
                }
            }
        }
        return a[0][0] * a[1][1] - a[0][1] * a[1][0];
    }
} Matrix;

Matrix transpose(const Matrix& m) {
    Matrix mt;
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            mt.matrix[i][j] = m.matrix[j][i];
        }
    }
    return mt;
}

Matrix invert(const Matrix& m) {
    Matrix alg_add;
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            alg_add.matrix[i][j] = m.minor(i,j) * ((i + j) % 2 == 0 ? 1 : -1);
        }
    }
    if(m.det() != 0)
        return transpose(alg_add) * (1 / m.det());

    return m;
}
Vector proj(const Vector& b,const Vector& a) {
    double pr = (b*a)/(b*b);
    Vector c = b * pr;
    return c;
}
void process_GS(Matrix & matr, Vector l) {
    if(l.x != 0 || l.y != 0 || l.z != 0) {
        Vector a1 = Vector(matr.matrix[0][0], matr.matrix[1][0], matr.matrix[2][0]);
        Vector a2 = Vector(matr.matrix[0][1], matr.matrix[1][1], matr.matrix[2][1]);
        Vector a3 = Vector(matr.matrix[0][2], matr.matrix[1][2], matr.matrix[2][2]);
        Vector b1 = a1;
        Vector b2 = a2 + (proj(b1, a2) * (-1));
        Vector b3 = a3 + (proj(b2, a3) * (-1)) + (proj(b2, a3) * (-1));
        matr.matrix[0][0] = b1.x / b1.length();
        matr.matrix[1][0] = b1.y / b1.length();
        matr.matrix[2][0] = b1.y / b1.length();
        matr.matrix[0][1] = b2.x / b2.length();
        matr.matrix[1][1] = b2.y / b2.length();
        matr.matrix[2][1] = b2.z / b2.length();
        matr.matrix[0][2] = b3.x / b3.length();
        matr.matrix[1][2] = b3.y / b3.length();
        matr.matrix[2][2] = b3.z / b3.length();
    }

}

typedef struct consts {
    double mass = 1; /* mass M */
    Matrix Ibody; /* Ibody */
    Matrix Ibodyinv; /* I−1 body (inverse of Ibody) */
    consts() {
        Ibody.matrix[0][0] = 1; Ibody.matrix[0][1] = 0; Ibody.matrix[0][2] = 0;
        Ibody.matrix[1][0] = 0; Ibody.matrix[1][1] = 1; Ibody.matrix[1][2] = 0;
        Ibody.matrix[2][0] = 0; Ibody.matrix[2][1] = 0; Ibody.matrix[2][2] = 1;
        Ibody = Ibody * (mass/12);
        Ibodyinv = invert(Ibody);
    }
} Constants;

typedef struct X {
    Vector C = Vector(100,100,200);//C' = v(t)
    Vector p = Vector(0,0,0);//p = F(t)
    Matrix R;//R' = ω(t)∗R(t)
    Vector L = Vector(0,0,0);//L' = τ(t)
    /* Derived quantities (auxiliary variables) */
    Matrix Iinv; /* I−1(t) */
    Vector v; /* v(t) */
    Vector omega; /* ω(t) */
    /* Computed quantities */
    Vector force; /* F(t) */
    Vector torque; /* τ(t) */

    X() = default;
    void calculate_derived(const Constants * data) {
        v = p / data->mass;
        Iinv = R * data->Ibodyinv * transpose(R);
        omega = Iinv * L;
    }
    X(Vector tC, Vector tp, Matrix tR, Vector tL) {
        C = tC;
        p = tp;
        R = tR;
        L = tL;
    }
    X(X const & new_X){
        C = new_X.C;
        p = new_X.p;
        R = new_X.R;
        L = new_X.L;
    }
    X & operator = (const X &new_X){
        C = new_X.C;
        p = new_X.p;
        R = new_X.R;
        L = new_X.L;
        return *this;
    }
    friend X operator + (const X& x1, const X& x2){
        return X(x1.C+ x2.C, x1.p + x2.p, x1.R+ x2.R, x1.L + x2.L);
    }
    friend X operator * (const X& x1, double & k){
        return X(x1.C * k, x1.p * k, x1.R * k, x1.L * k);
    }
} SystemState;


#endif //UNTITLED4_SYSTEM_H
