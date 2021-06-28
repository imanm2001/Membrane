#ifndef BEADINFO_H
#define BEADINFO_H
#include <QObject>
#include "utils/classes.h"
#include "Phys/bead.h"
#include "q3d/edge.h"
#include "Phys/bendingparameters.h"
#include <QMap>

QT_BEGIN_NAMESPACE

namespace Geometry {
class BeadInfo;
}
QT_END_NAMESPACE


class Geometry::BeadInfo:public Physics::Bead
{
private:
    QMap<Attributes,QObject*> *qmap;
public:   
    QVector<Geometry::Edge*> *_connections;
    QVector<Physics::BendingParameters*> *_bendingParameters;
    BeadInfo(QObject *, VecD3d *, double D, int ID);
    Geometry::Edge* createEdge(QObject *parent, BeadInfo* b);
    int findBeadIndexInTheConnection(Geometry::BeadInfo *);
    QObject* getAttribute(Attributes);
    void setAttribute(Attributes,QObject*);
};
static inline void DRALL(Physics::VecD3d *a ,Geometry::BeadInfo *b,Geometry::BeadInfo *c){
    a->_coords[0]=DR(b,c,0);
    a->_coords[1]=DR(b,c,1);
    a->_coords[2]=DR(b,c,2);
}
#endif // BEADINFO_H
