#ifndef VECD3D_H
#define VECD3D_H
#define det(x1,x2,y1,y2) x1*y2-x2*y1
#include <math.h>
#include <QVector3D>
#include <iostream>


#define DR(a,b,i) a->_coords->_coords[i]-b->_coords->_coords[i]
#define DRi(a,b,i) a->_b->_coords->_coords[i]-b->_b->_coords->_coords[i]
#define getTheOtherBead(e,b) (e->_vid1==b?e->_vid2:e->_vid1)
#define Delta(i,j) (i==j?1:0)
#define deltaDelta(i,j,l,m) (Delta(i,j)-Delta(l,m));




QT_BEGIN_NAMESPACE
namespace Physics {
class VecD3d;
}
QT_END_NAMESPACE

class Physics::VecD3d
{

public:
    double _coords[3];
    VecD3d();
    VecD3d(double x,double y,double z);
    VecD3d(double[3]);
    VecD3d(VecD3d*);
    VecD3d(QVector3D*);
    void setValues(VecD3d*);
    void setValues(double,double,double);
    double len();
    double dot(VecD3d*);
    void cross(VecD3d*,VecD3d*);
    void zero();
    void add(double x,double y,double z);
    void add(VecD3d);
    void subVec(VecD3d*,VecD3d*);
    void sub(double x,double y,double z);
    void sub(VecD3d*);
    void multConst(double);
    void nomilize();
    void print();
    void debug();
    ~VecD3d();
};



#endif // VECD3D_H
