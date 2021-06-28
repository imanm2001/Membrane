#ifndef BEAD_H
#define BEAD_H
#include "vecd3d.h"
#include <qobject.h>
#include "utils/classes.h"
#include <QDebug>
#define KbT 1

QT_BEGIN_NAMESPACE
namespace Physics {
class Bead;
}
QT_END_NAMESPACE


using Physics::VecD3d;
class Physics::Bead:public QObject
{

public:
    Bead(QObject*,VecD3d*,double D,int ID);
    void update(double dt);
    int ID;
    VecD3d *_coords;
    double _GAMMA;
    VecD3d *_force;

};


#endif // BEAD_H
