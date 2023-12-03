#ifndef ACCESSORY_H
#define ACCESSORY_H

#include <iostream>
#include <fstream>

using namespace std;

void create_csv_parameters(int n1,int n2,int nt,double tau,double gamma,double kappa){
    ofstream outfile;
    outfile.open("parameters.csv");
    outfile<<n1<<"\n"<<n2<<"\n"<<nt<<"\n"<<tau<<"\n"<<gamma<<"\n"<<kappa;
    outfile.close();
}

void create_csv_file(const double* A, int size, string filename){
    ofstream outfile;
    outfile.open(filename);
    for(int i=0;i<size;++i){
        outfile<<A[i];
        if(i<size-1){
            outfile<<"\n";
        }
    }
    outfile.close();
}

void create_bin_file(const double* A, int size, string filename){
    ofstream out(filename, ios::out | ios::binary);
    if(!out) {
        cout << "Cannot open file.";
        return;
    }

    out.write((char *) A, size*sizeof(double));
    out.close();
}


void prepare_pixels(double* rho, unsigned char* pixels, int pcount, double max){

    max=0;
    for(int i=0;i<pcount;++i){
        max=fmax(max,rho[i]);
    }
    for(int i=0;i<pcount;i++){
        double val=255*rho[i]/max;
        pixels[i]=fmin(255,val);
    }
}


#endif