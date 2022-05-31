#include "membrane.h"
#define PI 3.14159265358979323846
#define force 0.08
#define _K 50.01
#define _kappa 1.6

inline double threshold(double f){
    //return std::fabs(f)>(1e-3)?f:0;
    return f;

}
#define u14(le2,N,Ns,i) (le2*N->_coords[i]/Ns)
#define u23(dot1,dot2,N1,N2,le2,Ns1,Ns2,i) (dot1*N1->_coords[i]/(le2*Ns1)+dot2*N2->_coords[i]/(le2*Ns2))
#define ANG 0.001
/*
inline void curvV1ALL(Physics::Bead * b1,Physics::Bead * b2,Physics::Bead * b3,Physics::Bead * b4,
                      Physics::VecD3d *e2,Physics::VecD3d *N1,Physics::VecD3d *N2,double le2,
                      double Ns1,double Ns2,Physics::VecD3d *temp, double lambda){
    b1->_force->_coords[0]+=threshold(lambda*u14(le2,N1,Ns1,0));
    b1->_force->_coords[1]+=threshold(lambda*u14(le2,N1,Ns1,1));
    b1->_force._coords[2]+=threshold(lambda*u14(le2,N1,Ns1,2));


}

inline void curvV4ALL(Physics::Bead * b1,Physics::Bead * b2,Physics::Bead * b3,Physics::Bead * b4,
                      Physics::VecD3d *e2,Physics::VecD3d *N1,Physics::VecD3d *N2,double le2,
                      double Ns1,double Ns2,Physics::VecD3d *temp, double lambda){
    b4->_force._coords[0]+=threshold(lambda*u14(le2,N2,Ns2,0));
    b4->_force._coords[1]+=threshold(lambda*u14(le2,N2,Ns2,1));
    b4->_force._coords[2]+=threshold(lambda*u14(le2,N2,Ns2,2));

}

inline void curvV2ALL(Physics::Bead * b1,Physics::Bead * b2,Physics::Bead * b3,Physics::Bead * b4,
                      Physics::VecD3d *e2,Physics::VecD3d *N1,Physics::VecD3d *N2,double le2,
                      double Ns1,double Ns2,Physics::VecD3d *temp, double lambda){
    b1->_coords.subVec(&b3->_coords,temp);
    double dot1=temp->dot(e2);
    b4->_coords.subVec(&b3->_coords,temp);
    double dot2=temp->dot(e2);

    b2->_force._coords[0]+=threshold(lambda*u23(dot1,dot2,N1,N2,le2,Ns1,Ns2,0));
    b2->_force._coords[1]+=threshold(lambda*u23(dot1,dot2,N1,N2,le2,Ns1,Ns2,1));
    b2->_force._coords[2]+=threshold(lambda*u23(dot1,dot2,N1,N2,le2,Ns1,Ns2,2));

}
inline void curvV3ALL(Physics::Bead * b1,Physics::Bead * b2,Physics::Bead * b3,Physics::Bead * b4,
                      Physics::VecD3d *e2,Physics::VecD3d *N1,Physics::VecD3d *N2,double le2,
                      double Ns1,double Ns2,Physics::VecD3d *temp, double lambda){
    b1->_coords.subVec(&b2->_coords,temp);
    double dot1=temp->dot(e2);
    b4->_coords.subVec(&b2->_coords,temp);
    double dot2=temp->dot(e2);

    b3->_force._coords[0]-=threshold(lambda*u23(dot1,dot2,N1,N2,le2,Ns1,Ns2,0));
    b3->_force._coords[1]-=threshold(lambda*u23(dot1,dot2,N1,N2,le2,Ns1,Ns2,1));
    b3->_force._coords[2]-=threshold(lambda*u23(dot1,dot2,N1,N2,le2,Ns1,Ns2,2));

}
*/
Physics::Membrane::Membrane(double dt,double k,Geometry::Plane * p):_dt(dt)
  ,_norm(0,1*KbT),_t(0),_steps(0),_capture(false)

{
    _plane=p;
    _w=p->_columns;
    _h=p->_rows,

    _title=new QString(QString("%1_%2_%3").arg(QString::number(force),QString::number(_K),QString::number(_kappa)));
#ifdef OUT
    auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Diffusion\%1.txt)").arg(QUuid::createUuid().toString());
    s=s.replace("{","");
    s=s.replace("}","");

    file=new QFile(s);
    if (!file->open(QIODevice::WriteOnly | QIODevice::Text)){
        assert(0)   ;
        return;
    }
    out=new QTextStream(file);
#endif
    _lambda=0.5;
    _e1=new VecD3d();
    _e2=new VecD3d();
    _e3=new VecD3d();
    _E3=new VecD3d();
    _N1=new VecD3d();
    _N2=new VecD3d();
    _D=new VecD3d();
    _temp=new VecD3d();

    _dt*=0.025;
    _edges=_plane->_edges;

    _tris=&_plane->_tris;
    //_beads=p->_beads;
    _beads=new QVector<Physics::BeadInfo*>();
    _centeral=new QVector<Physics::BeadInfo*>();
    auto lb=_plane->_bb->first()->commonVertex(_plane->_bl->first());
    auto lt=_plane->_bt->first()->commonVertex(_plane->_bl->last());

    auto rb=_plane->_bb->last()->commonVertex(_plane->_br->first());
    auto rt=_plane->_bt->last()->commonVertex(_plane->_br->last());
    assert(lb!=nullptr);
    assert(lt!=nullptr);
    assert(rb!=nullptr);
    assert(rt!=nullptr);
    _allBeads=_plane->_beads;
    int l=std::sqrt(_allBeads->size());

    for(int i=0;i<_allBeads->size();i++){
        Physics::BeadInfo * b=_allBeads->at(i);
        Physics::BendingEnergy *be=new Physics::BendingEnergy(b);
        be->orderConnections();
        //if(b!=lb&&b!=rb&&b!=lt&&b!=rt)
        if(fabs(b->_coords->_coords[0])<1e-3){
            _centeral->append(b);
            if(fabs(b->_coords->_coords[1])<1e-3){
                _cb=b;

            }
        }

        if(!(containInBoundary(_plane->_bb,b)||containInBoundary(_plane->_bl,b)||containInBoundary(_plane->_br,b)||containInBoundary(_plane->_bt,b)))
        {



            // b->_coords._coords[2]=_norm(_generator)*0.05*b->_curv;
            //      b->_coords._coords[0]+=_norm(_generator)*0.01;
            //            b->_coords._coords[1]+=_norm(_generator)*0.01;
            _beads->append(b);
        }
    }
    _sf=new SpringForce(_K);
    _tether=new Tether(_K);



    double d[0];
    loadFromFile(d,1);
    applyCurvuture();
}
bool Physics::Membrane::containInBoundary(QVector<Geometry::Edge*> *ed,Physics::BeadInfo *b ){
    bool ret=false;
    for(int i=0;i<ed->size();i++){
        auto e=ed->at(i);
        if(e->contain(b)){
            ret=true;
            break;
        }
    }
    return ret;
}

void Physics::Membrane::update(){

    if(_steps%10000==0||_capture){

        //capture(_data,2);
    }



    applyCurvuture();


    for(auto b: *_beads){

        //   double nx=_norm(_generator),ny=_norm(_generator),nz=_norm(_generator);

        //b->_force.add(nx,ny,nz);
        _temp->zero();




        //b->dEdx(-1,_temp);

        ((Physics::BendingEnergy*)b->getAttribute(BENDINGENERGY))->Ep(-1,_temp);
        std::cout<<_temp->len()<<std::endl;
        _temp->multConst(-_kappa);

        if(b==_cb){

            _temp->_coords[2]+=force;

            if(_steps%1000==0){
                std::cout<<b->_coords->_coords[2]<<std::endl;

            }
        }
        //_temp->debug();
        b->_force->setValues(_temp);
        //b->dAdx(_K,_temp);
        //b->_force.add(_temp);

    }

    for(auto e:*_edges){
        _tether->eval(e);
    }
    for(auto b: *_beads){
        b->update(_dt);
    }
    /*
    for(auto b: *_beads){
    double phx=(b->_coords._coords[0])*10;
    double phy=b->_coords._coords[1]*10;
    b->_coords._coords[2]=0.1*(std::cos(phx+_t)+std::sin(phy+phx+_t));
    }
_t+=_dt*0.0002;*/
    /*

    if(_steps%100==0){
        for(int i=0;i<1;i++){
            _plane->updateTopology();
        }
    }
*/
    _steps++;
    for(auto tri:*_tris){
        tri->getNormal();

    }
#define NUM 5
#ifdef OUT
    if(_steps%1000==0){
        for(int i=0;i<NUM*NUM;i++){
            int index=(60-NUM)/2+(i/NUM)*60+i%NUM;
            auto bead=_beads->at(index);
            for(int j=0;j<3;j++)
                if(j==2&&i==NUM*NUM-1){

                    out->operator<<(bead->_coords._coords[j]);
                    out->operator<<("\r\n");
                }else{
                    out->operator<<(bead->_coords._coords[j]);
                    out->operator<<("\t");
                }
        }

        out->flush();

        std::cout<<_steps<<std::endl;
        assert(_steps<1200000);
        if(_steps>=1200000){
            exit(0);
        }
    }
#endif
}
/*
void Physics::Membrane::applyCurvutureOn(Geometry::Triangle* t1,Geometry::Triangle* t2,int eInd){
    auto e=t1->_e[eInd];
    Physics::BeadInfo * vids[4];
    int id=1;
    bool adj=false;
    for(int i=0;i<3;i++){
        auto v=t1->_v[i];
        if(e->contain(v)){
            if(id==2){
                if(adj){
                    vids[id]=v;
                }else{
                    vids[2]=vids[1];
                    vids[1]=v;
                }
            }else{
                vids[1]=v;
                adj=true;
            }

            id++;

        }else{
            vids[0]=v;
            adj=false;
        }
    }

    for(auto v:t2->_v){
        if(!e->contain(v)){

            vids[3]=v;

        }
    }

    DRALL(_e1,vids[0],vids[1]);
    DRALL(_E3,vids[0],vids[2]);
    _e1->cross(_E3,_N1);

    DRALL(_e1,vids[3],vids[2]);
    DRALL(_E3,vids[3],vids[1]);
    _e1->cross(_E3,_N2);

    DRALL(_e2,vids[2],vids[1]);


    double Ns1=_N1->dot(_N1),Ns2=_N2->dot(_N2);

    double N1=std::sqrt(Ns1),N2=std::sqrt(Ns2);
    double e2s=_e2->dot(_e2),le2=std::sqrt(e2s);

    _N1->cross(_N2,_temp);
    double sint2=std::sqrt(std::fmax(1-_N1->dot(_N2)/(N1*N2),0)/2)*(_temp->dot(_e2)>=0?1:-1);
    //std::cout<<std::sin(tetha/2)<<","<<sint2<<std::endl;
    //sint2=std::sin(std::asin(_temp->dot(_e2)/(N1*N2))/2);

    double K=_lambda*e2s/(N1+N2)*(sint2-0.05);
    if(N1>0&&N2>0&&le2>0){
        curvV1ALL(vids[0],vids[1],vids[2],vids[3],_e2,_N1,_N2,le2,Ns1,Ns2,_temp,K);

        curvV2ALL(vids[0],vids[1],vids[2],vids[3],_e2,_N1,_N2,le2,Ns1,Ns2,_temp,K);

        curvV3ALL(vids[0],vids[1],vids[2],vids[3],_e2,_N1,_N2,le2,Ns1,Ns2,_temp,K);
        curvV4ALL(vids[0],vids[1],vids[2],vids[3],_e2,_N1,_N2,le2,Ns1,Ns2,_temp,K);
    }
}*/
/*
void Physics::Membrane::applyCurvuture(){

    for(auto tri:*_tris){
        for(int i=0;i<3;i++){
            auto e=tri->_e[i];
            auto t2=e->_tris[1-e->triIndex(tri)];

            if(t2!=nullptr&&t2->_ID>tri->_ID){
                applyCurvutureOn(tri,t2,i);

            }
        }
    }
    */
void Physics::Membrane::applyCurvuture(){

    for(auto b:*_allBeads){
        ((BendingEnergy*)b->getAttribute(BENDINGENERGY))-> updateBendingParameters();
    }
}
/*
void Physics::Membrane::loadFromFile(){
    auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_force_%1.txt)").arg(_title);
    auto f=new QFile(s);
    if(f->exists()&&f->open(QIODevice::ReadOnly | QIODevice::Text)){

        auto fins=new QTextStream(f);
        int w,h;
        fins->operator>>(w);
        fins->operator>>(h);
        assert(w==_w);
        assert(h==_h);
        double x,y,z;
        for(int i=0;i<_beads->size();i++){
            auto b=_beads->at(i);
            fins->operator>>(x);
            fins->operator>>(y);
            fins->operator>>(z);
            b->_coords->setValues(x,y,z);
        }
    }
}
void Physics::Membrane::capture(){
    //auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_force_%1.txt)").arg(QUuid::createUuid().toString());
    //auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_force_%1_%2.txt)").arg(QString::number(force),QUuid::createUuid().toString());
    auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_force_%1.txt)").arg(_title);
    auto bak=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_force_%1_bak.txt)").arg(_title);
    s=s.replace("{","");
    s=s.replace("}","");
    std::cout<<s.toStdString()<<"\t"<<_cb->_coords->_coords[2]<<std::endl;


    if (QFile::exists(s))
    {
        if(QFile::exists(bak)){
            QFile::remove(bak);
        }
        QFile::copy(s,bak);
    }
    file=new QFile(s);


    if (!file->open(QIODevice::WriteOnly | QIODevice::Text)){
        assert(0)   ;
        return;
    }
    out=new QTextStream(file);
    out->operator<<(_w);
    out->operator<<("\t");
    out->operator<<(_h);
    out->operator<<("\r\n");
    for(int i=0;i<_beads->size();i++){
        auto b=_beads->at(i);
        out->operator<<(b->_coords->_coords[0]);
        out->operator<<("\t");
        out->operator<<(b->_coords->_coords[1]);
        out->operator<<("\t");
        out->operator<<(b->_coords->_coords[2]);
        out->operator<<("\r\n");
    }
    out->flush();
    file->close();
    _capture=false;
    delete out;
}
*/
