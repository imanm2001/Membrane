#ifndef PLANE_H
#define PLANE_H
#include "utils/classes.h"
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QGeometryRenderer>
#include <Qt3DRender/QBuffer>
#include <QVector3D>
#include <iostream>
#include <math.h>
#include "q3d/triangle.h"
#include <QVector>
#include "q3d/beadinfo.h"
#include "Phys/vecd3d.h"

#include <QObject>
#include <iostream>
#include <QSharedPointer>
#include <QRandomGenerator>

#include "utils/classes.h"
namespace Geometry{
class Plane: public QObject {

private:

    double _dw,_dh;
    QVector<Physics::VecD3d*> *_normals;
    float _triPoses[9];

    QRandomGenerator* _rand;



    Qt3DRender::QBuffer *_verBuff,*_indBuff;
    Qt3DRender::QGeometryRenderer * _mesh;
    void createPlan();
    void addQuadAt(int, int);
    void createTriangleByIndecies(int a,int b,int c,double,Geometry::Triangle*&);
    Geometry::Edge* createEdgeByIndecies(int a,int b);
    void rebuild();
    void updateNormals();
    bool updatable;
public:
    double _w,_h;
    int _rows,_columns;
    Physics::Membrane* _membrn;
    QVector<Geometry::Edge *> *_bb,*_bl,*_bt,*_br;
    QVector<Geometry::Edge *>* _edges;
    QVector<Physics::BeadInfo*>* _beads;
    QVector<Geometry::Triangle*> _tris;
    Plane(int rows,int columns,double w,double h);
    float* getTri(int i);

    Qt3DRender::QGeometryRenderer * mesh();
    void updateIndecies();
    void updateTopology();
    void updateVertecies();
    void updateVIN();
    QVector<Physics::BeadInfo*>* getBeads();
    QVector<Geometry::Edge*>* getEdges();
    void addNewEdge(Edge*);
    //Physics::Membrane* getMembraneInstance(double dt,double k);

    void expand();
    void swapEdge(Geometry::Triangle*,Geometry::Triangle*);
    ~Plane();
};
}



#endif // PLANE_H
