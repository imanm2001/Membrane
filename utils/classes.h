#ifndef CLASSES_H
#define CLASSES_H
#include <functional>
#include <cstdio>
#include <QObject>
#include <Qt3DRender/QGeometryRenderer>
#include <QFile>
#include <QVector>
#include <QString>
#define PI 3.14159265359
#define PI2 1.57079632679
/*
#include "Phys/beadinfo.h"
#include "Phys/bendingenergy.h"
#include "Phys/bendingparameters.h"
#include "Phys/tensor2.h"
#include "Phys/vecd3d.h"
#include "q3d/edge.h"
#include "q3d/triangle.h"*/


typedef std::function<bool()>  updateFunc;

enum Attributes{
    BENDINGENERGY=0,
    BORDER,
};


namespace Physics {
class Bead;
class Membrane;
class VecD3d;
class BendingParameters;
class BendingEnergy;
class Tensor2;
class Debuggable {
public:
    virtual void INIT()=0;
};

extern Debuggable *SIM;
}
#define IERROR(arg) \
    if ((arg)==0){ \
        std::cout<<"ERROR:\t"<<#arg<<std::endl;\
        std::cout<<__FILE__<<"\t:"<<__LINE__<<std::endl;\
        Physics::SIM->INIT();\
    }
namespace Geometry {
class QBoolean:QObject{

public:
    bool value;
    QBoolean(bool b):value(b){}
};

class BeadInfo;
class Edge;
class Triangle;
class Plane;
class Mesh{
public:
    virtual void addNewEdge(Edge*)=0;
};
}

#endif // CLASSES_H
