#ifndef MEMBRANE_H
#define MEMBRANE_H
#include <qvector.h>
#include <QFile>
#include <QTextStream>
#include "utils/classes.h"
#include "Phys/bendingenergy.h"
#include "q3d/edge.h"
#include "q3d/plane.h"
#include "q3d/triangle.h"
#include "q3d/beadinfo.h"
#include "Phys/springforce.h"
#include "Phys/tether.h"
#include "Phys/vecd3d.h"
#include "surfacewithphysics.h"

#include <random>
#include <QUuid>
QT_BEGIN_NAMESPACE
namespace Physics {
class Membrane;
}
QT_END_NAMESPACE
class Physics::Membrane:public SurfaceWithPhysics
{
private:
    std::default_random_engine _generator;
    std::normal_distribution<double> _norm;
    double _t;
    int _steps;
    QFile *file;
    QTextStream *out;
    VecD3d *_e1,*_e2,*_e3,*_E3,*_N1,*_N2,*_D,*_temp;
    Bead *_cb;

public:
    bool _capture;

    int _w,_h;
    double _dt,_lambda;
    Physics::SpringForce *_sf;
    Physics::Tether *_tether;
    Geometry::Plane * _plane;

    QVector<Geometry::Triangle*>* _tris;
    QVector<Geometry::Edge*>* _edges;
    QVector<Physics::BeadInfo*>* _beads,* _allBeads,*_centeral;
    void applyCurvutureOn(Geometry::Triangle*,Geometry::Triangle*,int);
    void applyCurvuture();

    bool containInBoundary(QVector<Geometry::Edge*> *ed,Physics::BeadInfo *b);
    Membrane(double dt,double k,Geometry::Plane*);
    void update();
    //void loadFromFile();
    //void capture();
    void saveObjToFile(QString *);

};

#endif // MEMBRANE_H
