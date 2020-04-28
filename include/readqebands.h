#ifndef READQEBANDS_H
#define READQEBANDS_H
#include <string>
#include <iostream>
#include <complex>
void readbands(double** bands,int kpoints,int bandnumber,std::string nscf);
void readoccupation(double** occupationnumber,int kpoints,int bandnumber,std::string nscf);
void readvmatrix(std::complex<double>*** kpoint_product,int kpoints,int bandnumber,std::string pmat);
void readkpoints(double** kpoints,double* kweight,int kpoints_count,double& volume,std::string outnscf);
void readdimension(int& kpoints,int& bandnumber,std::string pmat);
#endif
