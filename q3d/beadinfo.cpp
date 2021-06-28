#include "beadinfo.h"

Geometry::BeadInfo::BeadInfo(QObject* p,VecD3d* loc,double D,int ID):Bead(p,loc,D,ID)
{
    qmap=new QMap<Attributes,QObject*>();
    _connections=new QVector<Geometry::Edge*>();
    _bendingParameters=new QVector<Physics::BendingParameters*>();

}

int Physics::BeadInfo::findBeadIndexInTheConnection(Geometry::BeadInfo * b){
    int ret=-1;
    for(int i=0;i<_connections->size();i++){
        auto c=_connections->at(i);
        if(c->_vid1==b||c->_vid2==b){
            ret=i;
            break;
        }
    }
    return  ret;
}
QObject* Physics::BeadInfo::getAttribute(Attributes attr){
    return qmap->value(attr,nullptr);
}
void Physics::BeadInfo::setAttribute(Attributes key,QObject* val){
    qmap->insert(key,val);
}
