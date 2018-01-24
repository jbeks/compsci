#ifndef MYVECTOR_H
#define MYVECTOR_

#include <vector>
#include <stdexcept>
#include <math.h>
#include <iostream>

using namespace std;

class myVector {
    public:
        myVector(unsigned int);
        myVector(vector<double>);
        vector<double> get() const;
        unsigned int size() const;
        double operator[](int) const;
        myVector operator+(const double) const;
        myVector operator-(const double) const;
        myVector operator*(const double) const;
        myVector operator/(const double) const;
        myVector& operator=(const vector<double>&);
        myVector operator+(const myVector&) const;
        myVector operator-(const myVector&) const;
        myVector operator*(const myVector&) const;
        myVector operator/(const myVector&) const;
        myVector& operator=(const myVector&);
        double dot(myVector);
        double norm();
    private:
        vector<double> vec;
};

myVector::myVector(unsigned int k) {
    vector<double> v;
    for (unsigned int i = 0; i < k; i++) {
        v.push_back(0);
    }
    vec = v;
}

myVector::myVector(vector<double> v) {
    vec = v;
}

vector<double> myVector::get() const {
    return vec;
}

unsigned int myVector::size() const {
    return vec.size();
}

double myVector::operator[](int i) const {
    return vec[i];
}

myVector myVector::operator+(const double k) const {
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] + k);
    }
    return myVector(tmp);
}

myVector myVector::operator-(const double k) const {
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] - k);
    }
    return myVector(tmp);
}

myVector myVector::operator*(const double k) const {
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] * k);
    }
    return myVector(tmp);
}

myVector myVector::operator/(const double k) const {
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] / k);
    }
    return myVector(tmp);
}

myVector& myVector::operator=(const vector<double>& v) {
    vec = v;
    return *this;
}

myVector myVector::operator+(const myVector& v) const {
    if (vec.size() != v.size()) {
        throw logic_error("Can't add vectors of different sizes.");
    }
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] + v[i]);
    }
    return myVector(tmp);
}

myVector myVector::operator-(const myVector& v) const {
    if (vec.size() != v.size()) {
        throw logic_error("Can't subtract vectors of different sizes.");
    }
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] - v[i]);
    }
    return myVector(tmp);
}

myVector myVector::operator*(const myVector& v) const {
    if (vec.size() != v.size()) {
        throw logic_error("Can't multiply vectors of different sizes.");
    }
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] * v[i]);
    }
    return myVector(tmp);
}

myVector myVector::operator/(const myVector& v) const {
    if (vec.size() != v.size()) {
        throw logic_error("Can't divide vectors of different sizes.");
    }
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] / v[i]);
    }
    return myVector(tmp);
}

myVector& myVector::operator=(const myVector& v) {
    if (this == &v) {
        return *this;
    }
    vec = v.get();
    return *this;
}

double myVector::dot(myVector v) {
    if (vec.size() != v.size()) {
        throw logic_error(
            "Can't take dot product of vectors of different sizes."
        );
    }
    double tmp = 0;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp = tmp + vec[i] * v[i];
    }
    return tmp;
}

double myVector::norm() {
    return sqrt(this->dot(*this));
}

myVector operator+(double k, const myVector& v) {
    return v + k;
}

myVector operator*(double k, const myVector& v) {
    return v * k;
}

ostream& operator<<(ostream& os, const myVector& v) {
    os << "[";
    for (unsigned int i = 0; i < v.size(); i++) {
        if (v[i] >= 0) {
            os << " ";
        }
        os << v[i] << " ";
    }
    return os << "]";
}

#endif /* MYVECTOR_H */

