#ifndef WAVEFRONTOBJ_H
#define WAVEFRONTOBJ_H

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
#include <QFile>


QT_BEGIN_NAMESPACE
namespace Geometry {
class Plane;
class WaveFrontObj;
}

QT_END_NAMESPACE


class Geometry::WaveFrontObj:public QObject
{
public:
    QVector<Physics::VecD3d*> *_normals;
    WaveFrontObj(QString path);
    Qt3DRender::QGeometryRenderer * _mesh;
    QVector<Geometry::Edge *>* _edges;
    QVector<Physics::BeadInfo*>* _beads;
    QVector<Geometry::Triangle*> *_tris;
    Qt3DRender::QBuffer *_verBuff,*_indBuff;
    float* getTri(int i);
    Qt3DRender::QGeometryRenderer * mesh();
    Geometry::Edge* getEdge(Physics::BeadInfo*,Physics::BeadInfo*);
    void updateNormals();
    void updateIndecies();
    void updateVertecies();
    void buildMesh(Qt3DRender::QGeometry *customGeometry);
    void addNewEdge(Edge*);
    void saveToFile(QString *);
};

#endif // WAVEFRONTOBJ_H
