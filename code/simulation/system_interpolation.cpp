#include <Python.h>
#include <vector>
#include <myVector.h>
#include <math.h>

using namespace std;

static PyObject* interpolate(PyObject *self, PyObject *args);

// set methods included in module
static PyMethodDef system_interpolation_methods[] = {
    {"interpolate", interpolate, METH_VARARGS,
    "Interpolate positions in system"},
    {NULL, NULL, 0, NULL}
};

// set method parameters
static struct PyModuleDef system_interpolation_definition = {
    PyModuleDef_HEAD_INIT,
    "system_interpolation",
    "Python module with interpolation functions for nbody simulation",
    -1,
    system_interpolation_methods
};

// initialize module
PyMODINIT_FUNC PyInit_system_interpolation(void) {
    Py_Initialize();
    return PyModule_Create(&system_interpolation_definition);
}


void euler(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
);
void verlet(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
);
void rk2(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
);
void rk4(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
);
void hermite(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
);

// interpolation method for module
static PyObject* interpolate(PyObject *self, PyObject *args) {
    char *itype;
    double dt;
    double G;
    PyObject *obj;
    // read interpolation method, timestep,
    // gravity constant and simulation data object
    if (!PyArg_ParseTuple(args, "sddO", &itype, &dt, &G, &obj)) {
        PyErr_SetString(PyExc_Exception, "Invalid input.");
        return NULL;
    }
    // check if simulation data object is valid
    if (!PyList_Check(obj)) {
        PyErr_SetString(PyExc_Exception, "Invalid input.");
        return NULL;
    }
    PyObject *iter = PyObject_GetIter(obj);
    if (!iter) {
        PyErr_SetString(PyExc_Exception, "Could not create iterator.");
        return NULL;
    }
    // read all simulation data
    vector<double> mass;
    vector<myVector> pos;
    vector<myVector> vel;
    while (true) {
        PyObject *next = PyIter_Next(iter);
        // stop if all data has been read
        if (!next) { break; }
        // check for valid input
        if (!PyList_Check(next)) {
            PyErr_SetString(PyExc_Exception, "Invalid input.");
            return NULL;
        }
        if (PyList_Size(next) != 3) {
            PyErr_SetString(PyExc_Exception, "Invalid input.");
            return NULL;
        }
        PyObject *subiter = PyObject_GetIter(next);
        if (!subiter) {
            PyErr_SetString(PyExc_Exception, "Could not create iterator.");
            return NULL;
        }
        // read body data
        for (int i = 0; i < 3; i++) {
            PyObject *subnext = PyIter_Next(subiter);
            if (!subnext) {
                PyErr_SetString(PyExc_Exception, "Invalid input.");
                return NULL;
            }
            // first data is mass
            if (i == 0) {
                if (!PyFloat_Check(subnext)) {
                    PyErr_SetString(PyExc_Exception, "Invalid input.");
                    return NULL;
                }
                mass.push_back(PyFloat_AsDouble(subnext));
            }
            // other data are vectors
            else {
                if (!PyList_Check(subnext)) {
                    PyErr_SetString(PyExc_Exception, "Invalid input.");
                    return NULL;
                }
                PyObject *subsubiter = PyObject_GetIter(subnext);
                if (!subsubiter) {
                    PyErr_SetString(
                        PyExc_Exception, "Could not create iterator."
                    );
                    return NULL;
                }
                // collect vector data
                vector<double> v;
                while (true) {
                    PyObject *subsubnext = PyIter_Next(subsubiter);
                    if (!subsubnext) { break; }
                    if (!PyFloat_Check(subsubnext)) {
                        PyErr_SetString(PyExc_Exception, "Invalid input.");
                        return NULL;
                    }
                    v.push_back(PyFloat_AsDouble(subsubnext));
                }
                myVector myV = myVector(v);
                // fist vector (second data) is position
                if (i == 1) {
                    pos.push_back(myV);
                }
                // second vector (third data) is velocity
                else {
                    vel.push_back(myV);
                }
            }
        }
    }
    // apply given interpolation method
    if (!strcmp(itype, "euler")) {
        euler(dt, G, mass, pos, vel);
    }
    else if (!strcmp(itype, "verlet") || !strcmp(itype, "leapfrog")) {
        verlet(dt, G, mass, pos, vel);
    }
    else if (!strcmp(itype, "rk2")) {
        rk2(dt, G, mass, pos, vel);
    }
    else if (!strcmp(itype, "rk4")) {
        rk4(dt, G, mass, pos, vel);
    }
    else if (!strcmp(itype, "hermite")) {
        hermite(dt, G, mass, pos, vel);
    }
    else {
        PyErr_SetString(PyExc_Exception, "Unsupported interpolation type.");
        return NULL;
    }
    // build python object form position data
    PyObject *pos_lst = PyList_New(0);
    for (unsigned int i = 0; i < pos.size(); i++) {
        vector<double> vec = pos[i].get();
        PyObject *tmp = PyList_New(0);
        for (unsigned int j = 0; j < vec.size(); j++) {
            PyList_Append(tmp, PyFloat_FromDouble(vec[j]));
        }
        PyList_Append(pos_lst, tmp);
    }
    // build python object from velocity data
    PyObject *vel_lst = PyList_New(0);
    for (unsigned int i = 0; i < vel.size(); i++) {
        vector<double> vec = vel[i].get();
        PyObject *tmp = PyList_New(0);
        for (unsigned int j = 0; j < vec.size(); j++) {
            PyList_Append(tmp, PyFloat_FromDouble(vec[j]));
        }
        PyList_Append(vel_lst, tmp);
    }
    // create and return python object with positions and velocity
    PyObject *ret = PyTuple_New(2);
    PyTuple_SetItem(ret, 0, pos_lst);
    PyTuple_SetItem(ret, 1, vel_lst);
    return ret;
}

// calculate acceleration using Newtons laws of motion
myVector acceleration(
    unsigned int k, double G,
    vector<double> &m, vector<myVector> &p
) {
    myVector a = myVector(p[0].size());
    for (unsigned int i = 0; i < m.size(); i++) {
        if (i != k) {
            myVector vec = p[i] - p[k];
            a = a + vec * m[i] / pow(vec.norm(), 3);
        }
    }
    return a * G;
}

// calculate jerk (derivative of acceleration)
myVector jerk(
    unsigned int k, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
) {
    myVector j = myVector(p[0].size());
    for (unsigned int i = 0; i < m.size(); i++) {
        if (i != k) {
            myVector p_vec = p[i] - p[k];
            myVector v_vec = v[i] - v[k];
            double p_norm = p_vec.norm();
            j = j + m[i] * (v_vec / pow(p_norm, 3)
              - 3 * p_vec.dot(v_vec) * p_vec / pow(p_norm, 5));
        }
    }
    return j * G;
}

// implementations of different interpolation methods

void euler(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
) {
    vector<myVector> a;
    for (unsigned int i = 0; i < m.size(); i++) {
        a.push_back(acceleration(i, G, m, p));
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        p[i] = p[i] + v[i] * dt;
        v[i] = v[i] + a[i] * dt;
    }
}

void verlet(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
) {
    vector<myVector> a;
    for (unsigned int i = 0; i < m.size(); i++) {
        a.push_back(acceleration(i, G, m, p));
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        v[i] = v[i] + .5 * a[i] * dt;
        p[i] = p[i] + v[i] * dt;
    }
    a.clear();
    for (unsigned int i = 0; i < m.size(); i++) {
        a.push_back(acceleration(i, G, m, p));
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        v[i] = v[i] + .5 * a[i] * dt;
    }
}

void rk2(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
) {
    vector<myVector> p0 = p;
    vector<myVector> v0 = v;
    vector<myVector> a;
    for (unsigned int i = 0; i < m.size(); i++) {
        a.push_back(acceleration(i, G, m, p));
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        p[i] = p[i] + .5 * v[i] * dt;
        v[i] = v[i] + .5 * a[i] * dt;
    }
    a.clear();
    for (unsigned int i = 0; i < m.size(); i++) {
        a.push_back(acceleration(i, G, m, p));
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        p[i] = p0[i] + v[i] * dt;
        v[i] = v0[i] + a[i] * dt;
    }
}

void rk4(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
) {
    vector<myVector> p0 = p;
    vector<myVector> k1;
    vector<myVector> k2;
    vector<myVector> k3;
    vector<myVector> a;
    for (unsigned int i = 0; i < m.size(); i++) {
        a.push_back(acceleration(i, G, m, p));
        k1.push_back(a[i] * dt);
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        p[i] = p0[i] + .5 * (v[i] + .25 * k1[i]) * dt;
    }
    a.clear();
    for (unsigned int i = 0; i < m.size(); i++) {
        a.push_back(acceleration(i, G, m, p));
        k2.push_back(a[i] * dt);
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        p[i] = p0[i] + (v[i] + .5 * k2[i]) * dt;
    }
    a.clear();
    for (unsigned int i = 0; i < m.size(); i++) {
        a.push_back(acceleration(i, G, m, p));
        k3.push_back(a[i] * dt);
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        p[i] = p0[i] + (v[i] + .5 * (k1[i] + 2 * k2[i]) / 3) * dt;
        v[i] = v[i] + (k1[i] + 4 * k2[i] + k3[i]) / 6;
    }
}

void hermite(
    double dt, double G,
    vector<double> &m, vector<myVector> &p, vector<myVector> &v
) {
    vector<myVector> p0 = p;
    vector<myVector> v0 = v;
    vector<myVector> a0;
    vector<myVector> a1;
    vector<myVector> j0;
    vector<myVector> j1;
    for (unsigned int i = 0; i < m.size(); i++) {
        a0.push_back(acceleration(i, G, m, p));
        j0.push_back(jerk(i, G, m, p, v));
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        p[i] = p0[i] + v0[i] * dt + .5 * a0[i] * pow(dt, 2)
             + j0[i] * pow(dt, 3) / 6;
        v[i] = v0[i] + a0[i] * dt + .5 * j0[i] * pow(dt, 2);
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        a1.push_back(acceleration(i, G, m, p));
        j1.push_back(jerk(i, G, m, p, v));
    }
    for (unsigned int i = 0; i < m.size(); i++) {
        v[i] = v0[i] + .5 * (a0[i] + a1[i]) * dt
             + (j0[i] - j1[i]) * pow(dt, 2) / 12;
        p[i] = p0[i] + .5 * (v0[i] + v[i]) * dt
             + (a0[i] - a1[i]) * pow(dt, 2) / 12;
    }
}

