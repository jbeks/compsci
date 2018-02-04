#ifndef MYVECTOR_H
#define MYVECTOR_

#include <vector>
#include <stdexcept>
#include <math.h>
#include <iostream>

using namespace std;

// vector for easy vector calculations
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

// create empty vector of size k
myVector::myVector(unsigned int k) {
    vector<double> v;
    for (unsigned int i = 0; i < k; i++) {
        v.push_back(0);
    }
    vec = v;
}

// create vector from c-vector
myVector::myVector(vector<double> v) {
    vec = v;
}

// get vectors c-vector
vector<double> myVector::get() const {
    return vec;
}

// get vectors size
unsigned int myVector::size() const {
    return vec.size();
}

// get vecotr-value at index i
double myVector::operator[](int i) const {
    return vec[i];
}

// add double to vector
myVector myVector::operator+(const double k) const {
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] + k);
    }
    return myVector(tmp);
}

// subtract double from vector
myVector myVector::operator-(const double k) const {
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] - k);
    }
    return myVector(tmp);
}

// multiply vector by double
myVector myVector::operator*(const double k) const {
    vector<double> tmp;
    for (unsigned int i = 0; i < vec.size(); i++) {
        tmp.push_back(vec[i] * k);
    }
    return myVector(tmp);
}

// divide vector by diybke
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

// add vector to vector
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

// subtract vector from vector
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

// multiply vector by vector
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

// divide vector by vector
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

// assign vector to vector
myVector& myVector::operator=(const myVector& v) {
    if (this == &v) {
        return *this;
    }
    vec = v.get();
    return *this;
}

// calculate dotproduct of vector with vector
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

// calculate norm of vector
double myVector::norm() {
    return sqrt(this->dot(*this));
}

// add vector to double
myVector operator+(double k, const myVector& v) {
    return v + k;
}

// multiply double by vector
myVector operator*(double k, const myVector& v) {
    return v * k;
}

// insert vector into osstream
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

