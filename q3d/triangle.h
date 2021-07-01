#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "utils/classes.h"
#include "q3d/edge.h"
#include "Phys/vecd3d.h"
#include "q3d/beadinfo.h"
#include <iostream>
#include <QObject>
QT_BEGIN_NAMESPACE
namespace Geometry {
class Triangle;
}
QT_END_NAMESPACE
class Geometry::Triangle:public QObject
{

private:

    bool isCommonEdge(Geometry::Edge *);
    Physics::VecD3d *_location;
double _area;
public:
    int _ID;
    double _curvater;
    Physics::VecD3d* _norm;

    Geometry::BeadInfo* _v[3];
    Geometry::Edge * _e[3];
    double _angles[3];
    bool _isObtuse;
    Triangle(QObject*,int, BeadInfo* i0,BeadInfo* i1,BeadInfo* i2,double);
    Triangle(QObject*,int,Geometry::Edge *e1,Geometry::Edge *e2,Geometry::Edge *e3,double);
    float* getPositions();
    Geometry::Edge* commonEdge(Triangle *);
    void updateVIndecies();
    Geometry::Edge ** edgeWithVertex(BeadInfo* vid);
    Geometry::Edge * edgeWithVertexExclude(BeadInfo* vid,Geometry::Edge*);
    Physics::VecD3d *getNormal();
    int includeVertex(BeadInfo* b);
    void addToEdge(Geometry::Edge *e);
    void printEdges();
    double area();
    Physics::VecD3d* getLocation();
    int getVertexIndex(BeadInfo*);
    bool contains(BeadInfo*);
};


#endif // TRIANGLE_H
