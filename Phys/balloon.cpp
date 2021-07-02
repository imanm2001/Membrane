#include "balloon.h"
#define _P 0
#define _K 0
#define _kappa 10

Physics::Balloon::Balloon(double dt):SurfaceWithPhysics(),_dt(dt)
{
    generator=std::default_random_engine();
    distribution=std::normal_distribution<double>(0.0,1.0);
    _rand=_rand=QRandomGenerator::global();
    _dt*=1e-6;
    _sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\oval3_r1.obj)"));
    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\sphere2_r1.obj)"));
    double el=0;
    for(int i=0;i<_sphere->_edges->size();i++){
        auto e=_sphere->_edges->at(i);
        el+=e->_restLength;
        e->_restLength=0.39528;
    }
    std::cout<<">>>"<<el/_sphere->_edges->size()<<std::endl;
    _tris=_sphere->_tris;
    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\hex.obj)"));
    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\cap.obj)"));
    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\pyramid.obj)"));
    for(int i=0;i<_sphere->_tris->size();i++){
        _sphere->_tris->at(i)->_curvater=1/1.0;
    }
    _temp=new Physics::VecD3d();

    for(int i=0;i<_sphere->_beads->size();i++){
        auto b=_sphere->_beads->at(i);
        BendingEnergy *be=new BendingEnergy(b);
        be->orderConnections();
        be->updateBendingParameters();

    }
    _sf=new SpringForce(_K);
    _step=0;
    _l=_lsq=0;
    _title=new QString(QString("%1_%2_%3").arg(QString::number(_P),QString::number(_K),QString::number(_kappa)));
    _totalA=0;
    _scale=1e-3;
    _beads=_sphere->_beads;
    allocMem();
}
Qt3DRender::QGeometryRenderer * Physics::Balloon::mesh(){
    return _sphere->mesh();
}
void Physics::Balloon::update(){
    for(int i=0;i<_sphere->_tris->size();i++){
        auto tri=_sphere->_tris->at(i);
        tri->getNormal();
    }

    double H=111110,MH=0;
    _l=0;
    double A=0;
    for(int i=0;i<_sphere->_beads->size();i++){
        auto b=_sphere->_beads->at(i);

        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        be->_curv=1;

        be->orderConnections();
        be->updateBendingParameters();
        A+=be->_Av;
        if(be->_curv>0){
            H=std::fmin(H,be->_H);
            MH=std::fmax(MH,be->_H);

        }
        double len=b->_coords->len();
        _l+=len;
    }
    _l/=(double)_beads->size();

    if(_step<40000){
        _dt=1e-8;
    }else{
        _dt=(1e-8-2.5e-7)*std::exp(-(_step-40000)/40000.0)+2.5e-7;
    }
    _dt=1e-6;
    //double p=_step<100000?50:(50-_P)*std::exp(-(_step-100000)/100000.0)+_P;
    double p=_P;
    double kappa=_kappa*(1-std::exp(-_step/10.0));
    if(_step%500==0){
        std::cout<<H<<"\t"<<MH<<"\t"<<_l<<"\t"<<p<<"\t"<<kappa<<"\t"<<A<<"\t"<<_dt<<std::endl;
        //std::exit(-1);
    }
    //MC();
    _step++;


    for(int i=0;i<_sphere->_beads->size();i++){
        auto b=_sphere->_beads->at(i);
        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        _temp->zero();


        be->Ep2(-1,_temp);
        _temp->multConst(-kappa);
        if(be->_curv>0){
            //    _temp->print();
            b->_force->add(_temp);
        }



        _temp->setValues(be->_n);
        _temp->nomilize();
        _temp->multConst(p*be->_Av);
        b->_force->add(_temp);


    }

    //std::exit(-1);


    for(int i=0;i<_sphere->_edges->size();i++){
        auto edge=_sphere->_edges->at(i);
        edge->_restLength=0.3952;
        _sf->eval(edge);
    }
    double s=_sphere->_beads->size();
    bool upd=_step%100==0;
    _l=0;
    _lsq=0;
    for(int i=0;i<s;i++){
        auto b=_sphere->_beads->at(i);
        if(b->_GAMMA>0){
            b->update(_dt);

//            b->_coords->_coords[0]+=1.5e-3*(_rand->generateDouble()-0.5);
//            b->_coords->_coords[1]+=1.5e-3*(_rand->generateDouble()-0.5);
//            b->_coords->_coords[2]+=1.5e-3*(_rand->generateDouble()-0.5);

//            b->_coords->_coords[0]+=1e-4*(distribution(generator));
//            b->_coords->_coords[1]+=1e-4*(distribution(generator));
//            b->_coords->_coords[2]+=1e-4*(distribution(generator));
        }

        if(upd){
            double len=b->_coords->len();
            _l+=len/s;
            _lsq+=len*len/s;
        }
    }

    _step++;
}
double Physics::Balloon::calE(){
    double ret=0;

    for(int i=0;i<_beads->size();i++){
        auto b=_beads->at(i);
        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        //  be->orderConnections();
        be->updateBendingParameters();
    }
    for(int i=0;i<_beads->size();i++){
        auto b=_beads->at(i);
        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        ret+=be->calE();
    }
    ret*=_kappa;
    ret-=_P*volumne();
    return ret;
}

double Physics::Balloon::calE(int i){
    double ret=0;


    auto b=_beads->at(i);
    BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
    //be->orderConnections();
    be->updateBendingParameters();

    for(int j=0;j<b->_connections->size();j++){
        auto bj=getTheOtherBead(b->_connections->at(j),b);
        BendingEnergy *bej=(BendingEnergy *)bj->getAttribute(BENDINGENERGY);
        //bej->orderConnections();
        bej->updateBendingParameters();
    }
    for(int j=0;j<b->_connections->size();j++){
        auto bj=getTheOtherBead(b->_connections->at(j),b);
        BendingEnergy *bej=(BendingEnergy *)bj->getAttribute(BENDINGENERGY);
        ret+=bej->calE();

    }

    ret+=be->calE();
    ret*=_kappa;
    ret-=_P*volumne();
    return ret;
}
void Physics::Balloon:: updateVIN(){
    _sphere->updateVertecies();
    _sphere->updateIndecies();
    _sphere->updateNormals();
}
