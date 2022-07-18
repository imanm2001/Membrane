#ifndef EDGE_H
#define EDGE_H
#include <QVector>
#include <iostream>
#include <QObject>
#include "utils/classes.h"
#include "Phys/bead.h"
#include "q3d/beadinfo.h"


QT_BEGIN_NAMESPACE
namespace Geometry {
class Edge;
}
QT_END_NAMESPACE

class Geometry::Edge:public QObject
{

public:
    double _restLength,_restLengthScale;
    BeadInfo * _vid1,*_vid2;
    BeadInfo* commonVertex(Edge *e2);
    Geometry::Triangle * _tris[2];
    Edge(QObject*,BeadInfo*,BeadInfo*);
    int triIndex(Geometry::Triangle*);
    void print();
    bool contain(BeadInfo* vid);
    bool equal(Geometry::Edge *);
    bool equal(BeadInfo * b1,BeadInfo * b2);
    inline void updateRL();
    void location(Physics::VecD3d*);
    ~Edge();
};



#endif // EDGE_H
