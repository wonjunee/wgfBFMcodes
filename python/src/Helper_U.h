#ifndef HELPER_U
#define HELPER_U

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <vector>

namespace py = pybind11;
using namespace std;

void calculate_DUstar(  py::array_t<double, py::array::c_style | py::array::forcecast> DUstar_np,
                        py::array_t<double, py::array::c_style | py::array::forcecast> V_np,
                        py::array_t<double, py::array::c_style | py::array::forcecast> phi_np,
                        const int n1, const int n2, const double tau){
    py::buffer_info DUstar_buf= DUstar_np.request();
    py::buffer_info V_buf     = V_np.request();
    py::buffer_info phi_buf   = phi_np.request();
    double *DUstar            = static_cast<double *>(DUstar_buf.ptr);
    double *V                 = static_cast<double *>(V_buf.ptr);
    double *phi               = static_cast<double *>(phi_buf.ptr);
    
    double sum=0;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            // DUstar[i*n1+j]=exp(-(phi[i*n1+j])/tau - V[i*n1+j]);
            DUstar[i*n1+j]=exp(-phi[i*n1+j] - V[i*n1+j]);
            sum+=DUstar[i*n1+j];
        }
    }
    sum/=n1*n2;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            DUstar[i*n1+j]*=1.0/sum;
        }
    }
}

#endif