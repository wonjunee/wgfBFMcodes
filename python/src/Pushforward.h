#ifndef PUSHFORWARD_H
#define PUSHFORWARD_H

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>

namespace py = pybind11;
using namespace std;

class vector2d{
public:
    double x;
    double y;
    vector2d():x(0),y(0){};
    vector2d(double x_, double y_):x(x_),y(y_){};
}; // vector 2d

class Pushforward_mapping{
public:
    double* indices_x;
    double* indices_y;
    int n1;
    int n2;
    double dx;
    double dy;

    Pushforward_mapping(){
        indices_x=nullptr;
        indices_y=nullptr;
    }
    Pushforward_mapping(int n1, int n2){
        // initialize(n1,n2);
        this->n1=n1;
        this->n2=n2;
        dx=1.0/n1;
        dy=1.0/n2;
        indices_x=new double[n1*n2];
        indices_y=new double[n1*n2];
    }

    void initialize(int n1,int n2){
        this->n1=n1;
        this->n2=n2;
        dx=1.0/n1;
        dy=1.0/n2;
        indices_x=new double[n1*n2];
        indices_y=new double[n1*n2];
    }

    void destroy_all(){
        delete[] indices_x;
        delete[] indices_y;
    }

    // virtual ~Pushforward_mapping(){
    //     delete[] indices_x;
    //     delete[] indices_y;
    // }

    double cross(vector2d& o, vector2d& a, vector2d& b){
        return 1.0*(a.x-o.x)*(b.y-o.y) - 1.0*(a.y-o.y)*(b.x-o.x);
    }

    void color_the_area(double* final_rho, const double* initial_rho, vector<vector2d>& hull_vertices, double original_area){

        int left=fmax(0,fmin(hull_vertices[0].x,fmin(hull_vertices[1].x,fmin(hull_vertices[2].x,hull_vertices[3].x))));
        int right=fmin(n1-1,fmax(hull_vertices[0].x,fmax(hull_vertices[1].x,fmax(hull_vertices[2].x,hull_vertices[3].x))));
        int bottom=fmax(0,fmin(hull_vertices[0].y,fmin(hull_vertices[1].y,fmin(hull_vertices[2].y,hull_vertices[3].y))));
        int top=fmin(n2-1,fmax(hull_vertices[0].y,fmax(hull_vertices[1].y,fmax(hull_vertices[2].y,hull_vertices[3].y))));

        int ind=0;

        for(int i=bottom;i<=top;++i){
            for(int j=left;j<=right;++j){
                bool inside=true;
                vector2d point(j,i);
                for(int k=0;k<4;++k){
                    int k_next=k%4;
                    int k_current=(k-1+4)%4;
                    double res=cross(hull_vertices[k_next],hull_vertices[k_current],point);
                    double det_value=fmin(abs(hull_vertices[k_next].x-hull_vertices[k_current].x),abs(hull_vertices[k_next].y-hull_vertices[k_current].y));
                    // double det_value=0;
                    if(res>det_value+0.1){
                        inside=false;
                        break;
                    }
                }

                if(inside){
                    indices_x[ind]=point.x;
                    indices_y[ind]=point.y;
                    ind++;    
                    
                }
            }
        }

        for(int i=0;i<ind;++i){
            int ii=indices_y[i];
            int jj=indices_x[i];
            final_rho[ii*n1+jj]+=1.0*original_area/ind;        
        }
    }

    void run_pushforward_on_density(const double* initial_rho,double* final_rho,const double* pushforward_x,const double* pushforward_y){

        for(int i=0;i<n1*n2;++i){
            final_rho[i]=0;
        }

        int i,j;
        double x,y;

        vector<vector2d> hull_vertices;
        hull_vertices.push_back(vector2d(0,0));
        hull_vertices.push_back(vector2d(0,0));
        hull_vertices.push_back(vector2d(0,0));
        hull_vertices.push_back(vector2d(0,0));

        double x00, x01, x10, x11, x20, x21;
        double y00, y01, y02, y10, y11, y12;

        for(i=0;i<n2;++i){
            for(j=0;j<n1;++j){

                double original_area=initial_rho[i*n1+j];

                double gradx=0;
                double grady=0;

                if(original_area>0){
                    // create hull vertices

                    if(j==0){
                        x00=0;
                        x10=0;
                        x20=0;
                        x01=pushforward_x[int(fmax(0,(i-1))*n1+j)];
                        x11=pushforward_x[(i)*n1+j];
                        x21=pushforward_x[int(fmin(n2-1,(i+1))*n1+j)];
                    }else if(j==n1-1){
                        x00=pushforward_x[int(fmax(0,(i-1))*n1+(j-1))];
                        x10=pushforward_x[(i)*n1+(j-1)];
                        x20=pushforward_x[int(fmin(n2-1,(i+1))*n1+(j-1))];
                        x01=0;
                        x11=0;
                        x21=0;

                    }else{
                        x00=pushforward_x[int(fmax(0,(i-1))*n1+(j-1))];
                        x10=pushforward_x[(i)*n1+(j-1)];
                        x20=pushforward_x[int(fmin(n2-1,(i+1))*n1+(j-1))];
                        x01=pushforward_x[int(fmax(0,(i-1))*n1+j)];
                        x11=pushforward_x[(i)*n1+j];
                        x21=pushforward_x[int(fmin(n2-1,(i+1))*n1+j)];
                    }

                    if(i==0){
                        y00=0;
                        y01=0;
                        y02=0;
                        y10=pushforward_y[int(i*n1+fmax(0,(j-1)))];
                        y11=pushforward_y[i*n1+(j)];
                        y12=pushforward_y[int(i*n1+fmin(n1-1,(j+1)))];
                    }else if(i==n2-1){
                        y00=pushforward_y[int((i-1)*n1+fmax(0,(j-1)))];
                        y01=pushforward_y[(i-1)*n1+(j)];
                        y02=pushforward_y[int((i-1)*n1+fmin(n1-1,(j+1)))];
                        y10=0;
                        y11=0;
                        y12=0;
                    }else{
                        y00=pushforward_y[int((i-1)*n1+fmax(0,(j-1)))];
                        y01=pushforward_y[(i-1)*n1+(j)];
                        y02=pushforward_y[int((i-1)*n1+fmin(n1-1,(j+1)))];
                        y10=pushforward_y[int(i*n1+fmax(0,(j-1)))];
                        y11=pushforward_y[i*n1+(j)];
                        y12=pushforward_y[int(i*n1+fmin(n1-1,(j+1)))];
                    }


                    gradx=0.5*(x00+x10)/dx;
                    grady=0.5*(y00+y01)/dy;

                    hull_vertices[0].x=(j-0.5-gradx);
                    hull_vertices[0].y=(i-0.5-grady);    

                    gradx=0.5*(x01+x11)/dx;
                    grady=0.5*(y01+y02)/dy;

                    hull_vertices[1].x=(j+0.5-gradx);
                    hull_vertices[1].y=(i-0.5-grady);    

                    gradx=0.5*(x11+x21)/dx;
                    grady=0.5*(y11+y12)/dy;

                    hull_vertices[2].x=(j+0.5-gradx);
                    hull_vertices[2].y=(i+0.5-grady);    

                    gradx=0.5*(x10+x20)/dx;
                    grady=0.5*(y10+y11)/dy;

                    hull_vertices[3].x=(j-0.5-gradx);
                    hull_vertices[3].y=(i+0.5-grady); 

                    color_the_area(final_rho, initial_rho, hull_vertices, original_area);    
                }
            }
        }                 
    }
};



class poisson_solver{
public:
    fftw_plan planIn;
    fftw_plan planOut;
    double *workspace; // This will contain the coefficients for complex numbers
    double *kernel;

    int n1;
    int n2;
    double dx;
    double dy;

    double coeff;

    poisson_solver(){
        workspace=nullptr;
        kernel=nullptr;
    }

    poisson_solver(int n1, int n2) {
        // initialize(n1,n2);
        this->n1=n1;
        this->n2=n2;

        workspace =(double*) fftw_malloc(n1*n2*sizeof(double));

        planIn = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
        planOut = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);
        
        kernel=new double[n1*n2];
        create_negative_laplacian_kernel_2d();
    }

    void initialize(int n1, int n2) {
        this->n1=n1;
        this->n2=n2;

        workspace =(double*) fftw_malloc(n1*n2*sizeof(double));

        planIn = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
        planOut = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);
        
        kernel=new double[n1*n2];
        create_negative_laplacian_kernel_2d();
    }

    void create_negative_laplacian_kernel_2d(){
        int pcount=n1*n2;
        
        for(int i=0;i<n2;i++){
            for(int j=0;j<n1;j++){
                double negativeLaplacian=2.0*n1*n1*(1-cos(M_PI*(j)/n1)) + 2*n2*n2*(1-cos(M_PI*(i)/n2));
                kernel[i*n1+j]=negativeLaplacian;
            }
        }
    }

    void perform_inverse_laplacian(const double* push_mu, const double* DFstar){

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                workspace[i*n1+j]=-push_mu[i*n1+j]+DFstar[i*n1+j];
            }
        }

        fftw_execute(planIn);
        // workspace[0]=0;

        for(int i=0;i<n1*n2;++i){
            workspace[i]/=4*(n1)*(n2)*(1+coeff*kernel[i]);
        }

        fftw_execute(planOut);
    }


    void destroy_all(){
        delete[] workspace;
        delete[] kernel;
        fftw_destroy_plan(planIn);
        fftw_destroy_plan(planOut);
    }

    // virtual ~poisson_solver(){
    //     destroy_all();
    // }


}; // Poisson Solver

#endif