/**
 * 
 * c-transform using Fast Legendre Transform algorithm.
 *
 */

#ifndef FLT_BFM_H
#define FLT_BFM_H

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <iostream>

namespace py = pybind11;

class FLT2D{
public:
    int *indices;
    int hullCount;
    int n1;
    int n2;

    int* argmin;
    double* temp;

    FLT2D(int n1, int n2){
        this->n1 = n1;
        this->n2 = n2;

        int n=fmax(n1,n2);
        indices   = new int[n];
        hullCount = 0;

        argmin  = new int[n];    
        temp = new double[n1*n2];
    }

    ~FLT2D(){
        delete[] indices;
        delete[] argmin;
        delete[] temp;
    }  



    // void find_c_concave(const double* phi, double* phi_c, const double tau){
    //     for(int i=0;i<n2;++i){
    //         for(int j=0;j<n1;++j){
    //             double x=(j+0.5)/(1.0*n1);
    //             double y=(i+0.5)/(1.0*n2);
    //             temp[i*n1+j]=0.5*(x*x+y*y)-tau*phi[i*n1+j];
    //         }
    //     }

    //     compute_2d_dual(phi_c,temp);

    //     for(int i=0;i<n2;++i){
    //         for(int j=0;j<n1;++j){
    //             double x=(j+0.5)/(1.0*n1);
    //             double y=(i+0.5)/(1.0*n2);
    //             phi_c[i*n1+j]=1.0/tau*(0.5*(x*x+y*y)-phi_c[i*n1+j]);
    //         }
    //     }
    // }


    void find_c_concave(py::array_t<double, py::array::c_style | py::array::forcecast> phi_c_np, py::array_t<double, py::array::c_style | py::array::forcecast> phi_np, const double tau){
        // const double tau = 1;
        py::buffer_info phi_c_buf = phi_c_np.request();
        py::buffer_info phi_buf   = phi_np.request();
        double *phi_c             = static_cast<double *>(phi_c_buf.ptr);
        double *phi               = static_cast<double *>(phi_buf.ptr);
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                double x=(j+0.5)/(1.0*n1);
                double y=(i+0.5)/(1.0*n2);
                temp[i*n1+j]=0.5*(x*x+y*y)-tau*phi[i*n1+j];
            }
        }

        compute_2d_dual(phi_c,temp);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                double x=(j+0.5)/(1.0*n1);
                double y=(i+0.5)/(1.0*n2);
                phi_c[i*n1+j]=1.0/tau*(0.5*(x*x+y*y)-phi_c[i*n1+j]);
            }
        }
    }

    void convexify(double *phi, double *dual){
        compute_2d_dual(dual, phi);
        compute_2d_dual(phi, dual);
    }


    void transpose_doubles(double *transpose, double *data, int n1, int n2){    
        for(int i=0;i<n2;i++){
            for(int j=0;j<n1;j++){
             
                transpose[j*n2+i]=data[i*n1+j];
            }
        }
    }

    void compute_2d_dual(double *dual, const double *u){
        
        int pcount=n1*n2;
        
        int n=fmax(n1,n2);
        
        memcpy(temp, u, pcount*sizeof(double));
        
        for(int i=0;i<n2;i++){
            compute_dual(&dual[i*n1], &temp[i*n1], argmin, n1);
            
        }
        transpose_doubles(temp, dual, n1, n2);
        for(int i=0;i<n1*n2;i++){
            dual[i]=-temp[i];
        }
        for(int j=0;j<n1;j++){
            compute_dual(&temp[j*n2], &dual[j*n2], argmin, n2);
            
        }
        transpose_doubles(dual, temp, n2, n1);
    }



    void compute_dual(double *dual, double *u, int *dualIndicies, int n){
        
        get_convex_hull(u, n);
    
        compute_dual_indicies(dualIndicies, u, n);
        
        for(int i=0;i<n;i++){
            double s=(i+.5)/(n*1.0);
            int index=dualIndicies[i];
            double x=(index+.5)/(n*1.0);
            double v1=s*x-u[dualIndicies[i]];
            double v2=s*(n-.5)/(n*1.0)-u[n-1];
            if(v1>v2){
                dual[i]=v1;
            }else{
                dualIndicies[i]=n-1;
                dual[i]=v2;
            }
            
        }   
    }

    void get_convex_hull(double *u, int n){
        
        indices[0]=0;
        indices[1]=1;
        hullCount=2;
        
        for(int i=2;i<n;i++){    
            add_point(u, i);
        }
    }


    void add_point(double *u, int i){
        
        if(hullCount<2){
            indices[1]=i;
            hullCount++;
        }else{
            int hc=hullCount;
            int ic1=indices[hc-1];
            int ic2=indices[hc-2];
            
            double oldSlope=(u[ic1]-u[ic2])/(ic1-ic2);
            double slope=(u[i]-u[ic1])/(i-ic1);
            
            if(slope>=oldSlope){
                int hc=hullCount;
                indices[hc]=i;
                hullCount++;
            }else{
                hullCount--;
                add_point(u, i);
            }
        }
    }



    //Flavien: Q indicies or indices?
    void compute_dual_indicies(int *dualIndicies, double *u, int n){
        
        int counter=1;
        int hc=hullCount;
        
        for(int i=0;i<n;i++){
           
            double s=(i+.5)/(n*1.0);
            int ic1=indices[counter];
            int ic2=indices[counter-1];
            
            double slope=n*(u[ic1]-u[ic2])/(ic1-ic2);
            //printf("%f %f %d\n",s, slope, hc);
            while(s>slope&&counter<hc-1){
                counter++;
                ic1=indices[counter];
                ic2=indices[counter-1];
                slope=n*(u[ic1]-u[ic2])/(ic1-ic2);
            }
            dualIndicies[i]=indices[counter-1];
            
        }
    }



};

#endif