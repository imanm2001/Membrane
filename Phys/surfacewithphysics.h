#ifndef SURFACEWITHPHYSICS_H
#define SURFACEWITHPHYSICS_H
#include <QObject>

#include <Qt3DRender/QGeometryRenderer>
#include <QRandomGenerator>
#include <QFile>
#include <QVector>
#include <QString>
#include "q3d/beadinfo.h"
#include "utils/classes.h"
#define _KBT 6.21
QT_BEGIN_NAMESPACE
namespace Physics {
class SurfaceWithPhysics;
}
QT_END_NAMESPACE
class Physics::SurfaceWithPhysics:public QObject{
protected:
    QString *_title;
    QString *_shape;
    QVector<Geometry::BeadInfo*> *_beads;
    QVector<Physics::VecD3d*> *_preLocs;
    QVector<Geometry::Triangle*> *_tris;
    bool _capture;
    int _accepted,_MCsteps;
    double _scale,_oldE;
    double _KBTN;
    void writeBeadCoordinates(QVector<Geometry::BeadInfo*> *beads,QTextStream *out){
        out->setRealNumberPrecision(16);
        for(int i=0;i<beads->size();i++){
            auto b=beads->at(i);

            out->operator<<(b->_coords->_coords[0]);
            out->operator<<("\t");
            out->operator<<(b->_coords->_coords[1]);
            out->operator<<("\t");
            out->operator<<(b->_coords->_coords[2]);
            out->operator<<("\r\n");
        }
        out->flush();
    }
    void loadFromFile(double *data,int len){

        auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape_Scaled\Shape_%1_%2.txt)").arg(*_shape, *_title);
        std::cout<<s.toStdString()<<std::endl;
        auto f=new QFile(s);
        if(f->exists()&&f->open(QIODevice::ReadOnly | QIODevice::Text)){

            auto fins=new QTextStream(f);
            int s;
            fins->operator>>(s);
            assert(s==_beads->size());

            double x,y,z;
            for(int i=0;i<_beads->size();i++){
                auto b=_beads->at(i);
                fins->operator>>(x);
                fins->operator>>(y);
                fins->operator>>(z);
                b->_coords->setValues(x,y,z);
            }
            for(int i=0;i<len&&!fins->atEnd();i++){
                double d;
                fins->operator>>(d);
                data[i]=d;
            }
        }
    }
    void updateTris(){
        for(int i=0;i<_tris->size();i++){
            _tris->at(i)->getNormal();
        }
    }
    void capture(double *data,int len){
        bool error=false;
        for(int i=0;i<_beads->size();i++){
            double l=_beads->at(i)->_coords->len();
            if(l!=l&&l>100){
                error=true;
                break;
            }
        }
        if(!error){
        //auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_force_%1.txt)").arg(QUuid::createUuid().toString());
        //auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_force_%1_%2.txt)").arg(QString::number(force),QUuid::createUuid().toString());
        auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape_Scaled\Shape_%1_%2.txt)").arg(*_shape, *_title);
        std::cout<<s.toStdString()<<std::endl;
        auto bak=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape_Scaled\Shape_%1_%2_bak.txt)").arg(*_shape,*_title);



        if (QFile::exists(s))
        {
            if(QFile::exists(bak)){
                QFile::remove(bak);
            }
            QFile::copy(s,bak);
        }
        auto file=new QFile(s);


        if (!file->open(QIODevice::WriteOnly | QIODevice::Text)){
            assert(0)   ;
            return;
        }
        auto out=new QTextStream(file);
        out->operator<<(_beads->size());
        out->operator<<("\r\n");


        writeBeadCoordinates(_beads,out);
        for(int i=0;i<len;i++){
        out->operator<<(data[i]);
        out->operator<<("\r\n");
        }
        file->close();
        delete out;
        }
        _capture=false;
    }
public:
    SurfaceWithPhysics():QObject(){
        _capture=false;
        _shape=new QString("NULL");
        _scale=1;
        _accepted=0;
        _MCsteps=0;
        _preLocs=new QVector<Physics::VecD3d*>();
    }
    QString getShape(){
        return *_shape;
    }
    virtual Qt3DRender::QGeometryRenderer * mesh()=0;
    virtual void update()=0;
    virtual void updateVIN()=0;
    virtual double calE()=0;
    virtual double calE(int)=0;
    virtual QString getTitle(){
        return *_title;
    }
    virtual void saveToFile(){
        _capture=true;
    };
    double volumne(){
        double ret=0.0;
        for(int i=0;i<_tris->size();i++){
            auto tri=_tris->at(i);

            double v=(tri->_v[0]->_coords->_coords[0])*(tri->_v[1]->_coords->_coords[1])*(tri->_v[2]->_coords->_coords[2])
                    +(tri->_v[1]->_coords->_coords[0])*(tri->_v[2]->_coords->_coords[1])*(tri->_v[0]->_coords->_coords[2])
                    +(tri->_v[2]->_coords->_coords[0])*(tri->_v[0]->_coords->_coords[1])*(tri->_v[1]->_coords->_coords[2])
                    -(tri->_v[1]->_coords->_coords[0])*(tri->_v[0]->_coords->_coords[1])*(tri->_v[2]->_coords->_coords[2])
                    -(tri->_v[0]->_coords->_coords[0])*(tri->_v[2]->_coords->_coords[1])*(tri->_v[1]->_coords->_coords[2])
                    -(tri->_v[2]->_coords->_coords[0])*(tri->_v[1]->_coords->_coords[1])*(tri->_v[0]->_coords->_coords[2]);
            v/=6.0;
            v=fabs(v);

            if(tri->getLocation()->dot(tri->_norm)<0){
                v*=-1;
            }
            ret+=v;
        }
        return ret;
    }
    void allocMem(){
        int size=_beads->size();
        for(int i=0;i<size;i++){
            _preLocs->append(new Physics::VecD3d(_beads->at(i)->_coords));
        }
        _oldE=0;
        _KBTN=_KBT*_beads->size();
    }

    void MC(){
        _MCsteps++;

        if(_MCsteps%100==0){
            double noE=calE();
            _accepted=(_accepted*2)/_beads->size();
            std::cout<<_accepted<<"\t"<<_scale<<"\t"<<noE-_oldE<<std::endl;

            if(_accepted>55){
                _scale*=1.05;
            }else if(_accepted<45){
                _scale*=0.95;
            }
            _accepted=0;
            _oldE=noE;
        }
        auto rand=QRandomGenerator::global();
        for(int ii=0;ii<_beads->size();ii++){
            auto b=_beads->at(ii);
            double oE=calE(ii);
            for(int j=0;j<3;j++){
                double s=(rand->generateDouble()-0.5)*_scale;
                b->_coords->_coords[j]+=s;
            }

            double nE=calE(ii);
            //double dE=nE-_oldE;
            double dE=nE-oE;
            //if(rand->generateDouble()>std::exp((nE-_oldE)/_KBT)){
            //if(dE<0||rand->generateDouble()<std::exp(-dE/_KBT)){
            if(dE<0){
                _preLocs->at(ii)->setValues(b->_coords);
                _accepted++;
                _oldE=calE();
            }else{
                b->_coords->setValues(_preLocs->at(ii));
            }
        }
    }

/*
    void MC(){
        _MCsteps++;
        auto rand=QRandomGenerator::global();
        for(int i=0;i<_beads->size();i++){
            auto b=_beads->at(i);
            for(int j=0;j<3;j++){
                double s=(rand->generateDouble()-0.5)*_scale;
                b->_coords->_coords[j]+=s;
            }
        }
        double nE=calE();
        double dE=nE-_oldE;


        if(_MCsteps%100==0){
            std::cout<<_accepted<<"\t"<<_scale<<"\t"<<(nE-_oldE)/_KBTN<<"\t"<<(nE-_oldE)<<"\t"<<_beads->size() <<std::endl;
            if(_accepted>60){
                //_scale*=1.25;
            }else if(_accepted<40){
                //_scale*=0.75;
            }
            _accepted=0;
        }

        if(dE<0||rand->generateDouble()<std::exp(-dE/_KBT)){
            for(int i=0;i<_beads->size();i++){
                auto b=_beads->at(i);
                auto v=_preLocs->at(i);
                v->setValues(b->_coords);
            }
            _accepted++;
            _oldE=nE;
        }else{
            for(int i=0;i<_beads->size();i++){
                auto b=_beads->at(i);
                auto v=_preLocs->at(i);
                b->_coords->setValues(v);
            }
        }

    }
*/
};

#endif // SURFACEWITHPHYSICS_H
