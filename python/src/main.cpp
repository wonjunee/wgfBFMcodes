#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <vector>
#include "Pushforward.h"
#include "Accessory.h"
#include "FLT.h"
#include "Helper_U.h"
#include "BFM.h"

namespace py = pybind11;
using namespace std;

void calculate_gradient(py::array_t<double, py::array::c_style | py::array::forcecast> gradx_np,
                        py::array_t<double, py::array::c_style | py::array::forcecast> grady_np,
                        py::array_t<double, py::array::c_style | py::array::forcecast> u_np,
                        const int n1, const int n2){
    py::buffer_info gradx_buf = gradx_np.request();
    py::buffer_info grady_buf = grady_np.request();
    py::buffer_info u_buf     = u_np.request();
    double *gradx             = static_cast<double *>(gradx_buf.ptr);
    double *grady             = static_cast<double *>(grady_buf.ptr);
    double *u                 = static_cast<double *>(u_buf.ptr);

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1-1;++j){
            gradx[i*n1+j]=1.0*n1*(u[i*n1+j+1]-u[i*n1+j]);
        }
    }
    for(int j=0;j<n1;++j){
        for(int i=0;i<n2-1;++i){
            grady[i*n1+j]=1.0*n2*(u[(i+1)*n1+j]-u[i*n1+j]);
        }
    }
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


PYBIND11_MODULE(bfmgf, m) {
    // optional module docstring
    m.doc() = "C++ wrapper for Fast diffusion (Wasserstein gradient flows) code";

    py::class_<FLT2D>(m, "FLT2D")
        .def(py::init<int, int> () )
        .def("find_c_concave", &FLT2D::find_c_concave);

    py::class_<BFM>(m, "BFM")
        .def(py::init<int, int, double> () )
        .def("compute_push_forth", &BFM::compute_push_forth)
        .def("compute_pull_back", &BFM::compute_pull_back)
        .def("update_sigma", &BFM::update_sigma);

    m.def("calculate_DUstar", &calculate_DUstar);
    m.def("calculate_gradient", &calculate_gradient);
}
