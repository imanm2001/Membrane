#ifndef CTS_H
#define CTS_H
#include <qobject.h>
#include <Phys/vecd3d.h>
#include "utils/classes.h"
#include <gsl/gsl_linalg.h>
#include <QRandomGenerator>
#endif // CTS_H
QT_BEGIN_NAMESPACE

namespace Physics {
class CTS;
}
QT_END_NAMESPACE


using Physics::VecD3d;
class Physics::CTS:public QObject
{
private :

    gsl_matrix *_res,*_coeffs,*_B,*_gslMat;
    gsl_vector *_temp1,*_temp2,*_strainV,*_d;
    int _dim;
    VecD3d **_vs1,**_vs2;
    const int INDEX[6]={2,1,3,1,3,2};
    gsl_permutation * p;
    QRandomGenerator* _rand;
    int _Brows,_Bcolumns;
public:
    void printM(gsl_matrix *m){
        int index=0;
        for(int i=0;i<m->size1;i++){
            for(int j=0;j<m->size2;j++){
                std::cout<<m->data[index++]<<"\t";
            }
            std::cout<<std::endl;
        }
    }
    CTS(int dim):_dim(dim){
        _rand=QRandomGenerator::global();
        p= gsl_permutation_alloc (dim+1);
        _gslMat=gsl_matrix_alloc(dim+1,dim+1);
        _vs1=(VecD3d**)malloc(sizeof(VecD3d*)*(dim+1));
        _vs2=(VecD3d**)malloc(sizeof(VecD3d*)*(dim+1));
        double *mat=_gslMat->data;
        if(dim==3){
            mat[0]=mat[4]=mat[8]=mat[12]=1;
        }else{
            mat[0]=mat[3]=mat[6]=1;
        }

        _res=gsl_matrix_alloc(dim+1,dim+1);
        _coeffs=gsl_matrix_alloc(dim+1,dim+1);
        _temp1=gsl_vector_alloc(dim+1);
        _temp2=gsl_vector_alloc(dim+1);


        _Brows=dim+(dim*(dim-1))/2;
        _Bcolumns=dim*(dim+1);
        _B=gsl_matrix_alloc(_Brows,_Bcolumns);
        _d=gsl_vector_alloc(_Bcolumns);
        _strainV=gsl_vector_alloc(_Brows);
        std::fill(_B->data,_B->data+_Brows*_Bcolumns,0.0);
    }
    void setVecToVec3D(VecD3d src,gsl_vector *des){
        des->data[0]=src._coords[0];
        des->data[1]=src._coords[1];
        des->data[2]=src._coords[2];
    }
    void setVecToVec2D(VecD3d src,gsl_vector *des){
        des->data[0]=src._coords[0];
        des->data[1]=src._coords[2];

    }
    bool getNcoeffs3D(VecD3d **vs,gsl_matrix * coeffs){


        for(int i=0;i<4;i++){
            for(int j=1;j<4;j++){
                _gslMat->data[j+i*4]=vs[i]->_coords[j-1];
            }
        }
        if(fabs(gsl_linalg_LU_det(_gslMat,1))>0){
            int s;

            gsl_linalg_LU_decomp (_gslMat, p, &s);
            gsl_linalg_LU_invert(_gslMat,p,_res);
            for(int i=0;i<4;i++){
                std::fill(_temp1->data,_temp1->data+4,0.0);
                _temp1->data[i]=1;
                gsl_blas_dgemv(CblasNoTrans,1,_res,_temp1,0,_temp2);
                std::copy(_temp2->data,_temp2->data+4,coeffs->data+4*i);
            }
            return true;
        }else{
            return false;
        }
    }
    void compactifyVs3D(VecD3d *v1,VecD3d *v2,VecD3d *v3,VecD3d *v4,VecD3d **vs){
        vs[0]=v1;
        vs[1]=v2;
        vs[2]=v3;
        vs[3]=v4;
    }
    bool getStrainMatrix3D(VecD3d **vI,VecD3d **vF,gsl_matrix *res){
        if(getNcoeffs3D(vI,_coeffs)){
            for(int n=0;n<4;n++){
                _d->data[n*3]=vF[n]->_coords[0]-vI[n]->_coords[0];
                _d->data[n*3+1]=vF[n]->_coords[1]-vI[n]->_coords[1];
                _d->data[n*3+2]=vF[n]->_coords[2]-vI[n]->_coords[2];
            }

            for(int n=0;n<3;n++){


                ///
                _B->data[n+n*12]=_coeffs->data[1+n];
                _B->data[n+n*12+3]=_coeffs->data[5+n];
                _B->data[n+n*12+6]=_coeffs->data[9+n];
                _B->data[n+n*12+9]=_coeffs->data[13+n];
                ///
                ///
                _B->data[INDEX[n*2+1]-1+n*12+36]=_coeffs->data[INDEX[n*2]];
                _B->data[INDEX[n*2]-1+n*12+36]=_coeffs->data[INDEX[n*2+1]];

                _B->data[INDEX[n*2+1]-1+n*12+39]=_coeffs->data[INDEX[n*2]+4];
                _B->data[INDEX[n*2]-1+n*12+39]=_coeffs->data[INDEX[n*2+1]+4];

                _B->data[INDEX[n*2+1]-1+n*12+42]=_coeffs->data[INDEX[n*2]+8];
                _B->data[INDEX[n*2]-1+n*12+42]=_coeffs->data[INDEX[n*2+1]+8];

                _B->data[INDEX[n*2+1]-1+n*12+45]=_coeffs->data[INDEX[n*2]+12];
                _B->data[INDEX[n*2]-1+n*12+45]=_coeffs->data[INDEX[n*2+1]+12];
            }



            gsl_blas_dgemv(CblasNoTrans,1,_B,_d,0,_strainV);

            std::fill(res->data,res->data+9,0);


            for(int n=0;n<3;n++){
                res->data[n+3*n]=_strainV->data[n];
                int i=INDEX[n*2+1]-1;
                int j=INDEX[n*2]-1;
                res->data[j*3+i]=res->data[i*3+j]=_strainV->data[n+3];
            }
            return true;
        }else{
            return false;
        }
    }


    bool getNcoeffs2D(VecD3d **vs,gsl_matrix * coeffs){
        const int CORDS[]={0,2};

        for(int i=0;i<3;i++){
            for(int j=1;j<3;j++){
                _gslMat->data[j+i*3]=vs[i]->_coords[CORDS[j-1]];
            }
        }


        int s;

        gsl_linalg_LU_decomp (_gslMat, p, &s);
        gsl_linalg_LU_invert(_gslMat,p,_res);
        for(int i=0;i<3;i++){
            std::fill(_temp1->data,_temp1->data+3,0.0);
            _temp1->data[i]=1;
            gsl_blas_dgemv(CblasNoTrans,1,_res,_temp1,0,_temp2);
            std::copy(_temp2->data,_temp2->data+3,coeffs->data+3*i);
        }


        return true;
    }
    void compactifyVs2D(VecD3d *v1,VecD3d *v2,VecD3d *v3,VecD3d **vs){
        vs[0]=v1;
        vs[1]=v2;
        vs[2]=v3;

    }
    bool getStrainMatrix2D(VecD3d **vI,VecD3d **vF,gsl_matrix *res){
        if(getNcoeffs2D(vI,_coeffs)){

            for(int n=0;n<3;n++){
                _d->data[n*2]=vF[n]->_coords[0]-vI[n]->_coords[0];
                _d->data[n*2+1]=vF[n]->_coords[2]-vI[n]->_coords[2];

            }

            for(int n=0;n<2;n++){


                ///
                _B->data[n+n*6]=_coeffs->data[1+n];
                _B->data[n+n*6+2]=_coeffs->data[4+n];
                _B->data[n+n*6+4]=_coeffs->data[7+n];
            }
            ///
            ///
            _B->data[INDEX[1]-1+12]=_coeffs->data[INDEX[0]];
            _B->data[INDEX[0]-1+12]=_coeffs->data[INDEX[1]];

            _B->data[INDEX[1]-1+14]=_coeffs->data[INDEX[0]+3];
            _B->data[INDEX[0]-1+14]=_coeffs->data[INDEX[1]+3];

            _B->data[INDEX[1]-1+16]=_coeffs->data[INDEX[0]+6];
            _B->data[INDEX[0]-1+16]=_coeffs->data[INDEX[1]+6];






            gsl_blas_dgemv(CblasNoTrans,1,_B,_d,0,_strainV);

            std::fill(res->data,res->data+4,0);


            for(int n=0;n<2;n++){
                res->data[n+2*n]=_strainV->data[n];

            }

            int i=INDEX[1]-1;
            int j=INDEX[0]-1;
            res->data[j*2+i]=res->data[i*2+j]=_strainV->data[2];
            return true;
        }else{
            return false;
        }
    }



};
