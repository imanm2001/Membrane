#include "membranefromobj.h"
#define _P 1
#define _K 50
#define _kappa 18522
#define _T 0
#define _E 2*_K/1.73205081
#define _THRESHOLD 42
#define _F 4010
Physics::MembraneFromObj::MembraneFromObj(double dt):SurfaceWithPhysics(),_dt(dt),_appliedF(_F),_py(500),_frad(55),_alpha(0)
{
    std::cout<<"Load From Files"<<std::endl;
    _scale=1;
    //Mesh Project
    _Rind=8;
    _cb=nullptr;
    _kappaFactor=1;
    _radiusFactor=1;
    std::cout<<"Load output"<<std::endl;
    _disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\disc_r44_d60_relaxed.obj)"));
    std::cout<<"Load Destination"<<std::endl;
    _discDest=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\disc_r44_d50_relaxed.obj)"));
    std::cout<<"All files are loaded"<<std::endl;
    _tris=_disc->_tris;

    _rand=_rand=QRandomGenerator::global();
    _border=new QVector<Geometry::BeadInfo*>();
    _xprofile=new QVector<Geometry::BeadInfo*>();
    _beads=_disc->_beads;
    _beadsDestI=_discDest->_beads;
    _beadsDestF=new QVector<Geometry::BeadInfo*>();
    _beadsDestTransiation=new QVector<Physics::Bead*>();
    _updatable=new QVector<Geometry::BeadInfo*>();
    _quads=new QVector<Physics::Quad*>();
    SurfaceWithPhysics::_shape=new QString("Disc");

    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\hex.obj)"));
    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\cap.obj)"));
    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\pyramid.obj)"));
    for(int i=0;i<_disc->_tris->size();i++){
        auto tri=_disc->_tris->at(i);
        if(tri->getLocation()->len()>_THRESHOLD){
            tri->_curvater=1;
        }else{
            tri->_curvater=1;
        }
    }
    _temp=new Physics::VecD3d();
    _temp2=new Physics::VecD3d();
    double mRdis=0;
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        if(fabs(b->_coords->_coords[2])<1e-3){
            double l=b->_coords->_coords[0];
            if(_xprofile->size()==0||std::fabs(l-_THRESHOLD)<std::fabs(mRdis-_THRESHOLD)){
                _Rind=_xprofile->size();
                mRdis=l;
            }
            _xprofile->append(b);
            if(fabs(b->_coords->_coords[0])<1e-3){
                _cb=b;
            }
        }

        //        b->_coords->_coords[0]+=1e-4*(_rand->generateDouble()-0.5);
        //        b->_coords->_coords[1]+=1e-4*(_rand->generateDouble()-0.5);
        //        b->_coords->_coords[2]+=1e-4*(_rand->generateDouble()-0.5);

    }
    _initalArea=0;
    _APV=0;
    updateTris();
    _TIE=0;
    int N=0;
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        BendingEnergy *be=new BendingEnergy(b);

        be->orderConnections();
        be->updateBendingParameters();
        be->_curv=0.048;
        if(b->_coords->len()<_THRESHOLD){
            _TIE+=_kappa*be->calE();
        }
        double r=b->_coords->len();

        if(b->getAttribute(BORDER)!=nullptr){
            _border->append(b);
            be->_border=true;
        }else {
            if(r<_THRESHOLD){
                N++;
                _initalArea+=be->_Av;
                double d=r/_THRESHOLD;
                //b->_coords->_coords[1]=10*(1-d*d);
            }

            _updatable->append(b);
        }
        _APV=_initalArea/(double)N;

    }
    std::cout<<"TOTAL::"<<_initalArea<<std::endl;
    //    _sf=new SpringForce(_K);
    _sf=new SpringForce(10);
    _tet=new Tether(100000);
    _step=0;

    SurfaceWithPhysics::_title=new QString(QString("%1_%2_%3_%4").arg(QString::number(_P),QString::number(_K),QString::number(_kappa),QString::number(_F)));
    _frad=1;
    for(int i=0;i<_border->size();i++){
        auto b=_border->at(i);
        double l=b->_coords->len();
        if(i==0||l>_frad){
            _frad=l;
        }
        //b->_coords->multConst(10.0/b->_coords->len());
    }



    _pBE=0;
    allocMem();
    std::cout<<volumne()<<std::endl;



    auto bb=_xprofile->at(8);
    _kappaFactor=bb->_coords->_coords[0]/2;
    _kappaFactor=2*_P*_kappaFactor*_kappaFactor*_kappaFactor/_kappa;
    _radiusFactor=(bb->_coords->_coords[0]-1e-2)/_THRESHOLD;
    _radiusFactor=1;
    std::cout<<"---"<<_kappaFactor<<"\t"<<_radiusFactor<<std::endl;
    for(int i=0;i<_beads->size();i++){
        auto b=_beads->at(i);
        auto q=new Quad(b,1e3,&_scale);
        _quads->append(q);
    }

    for(int i=0;i<_beadsDestI->size();i++){
        auto b=new Physics::Bead((Physics::Bead*) _beadsDestI->at(i));
        _beadsDestTransiation->append(b);
        auto bi=new Geometry::BeadInfo(_beadsDestI->at(i));
        _beadsDestF->append(bi);
    }

    loadFromFile(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape_Scaled_fixed_radial_force_r44_d50\Shape_%1_%2.txt)").arg(*_shape, *_title),_beadsDestF);
    std::cout<<"Quads"<<std::endl;
    for(int i=0;i<_quads->size();i++){
        auto q=_quads->at(i);
        if(i%1000==0){
            std::cout<<std::round(100*i/(double)_beads->size())<<std::endl;
        }
        for(int j=0;j<_beadsDestTransiation->size();j++){
            auto b=_beadsDestTransiation->at(j);
            q->addNewCandidate(b);
        }
    }
    std::cout<<"100"<<std::endl;
    loadFromFile();
}

void Physics::MembraneFromObj::updateBeads(QVector<Geometry::BeadInfo*> *beads,double P,double kappa,double dt){




    _temp->zero();
    _temp2->zero();

    int j=_step/1000;
    for(int i=0;i<_quads->size();i++){
        auto q=_quads->at(i);
        q->eval();
    }
    for(int i=0;i<_disc->_edges->size();i++){
        auto e=_disc->_edges->at(i);
        _sf->eval(e->_vid1,e->_vid2,0);
    }
    for(int i=0;i<beads->size();i++){
        auto b=beads->at(i);
        b->update(dt);

        if(j%10==0){
            b->_coords->_coords[0]+=1e-2*(_rand->generateDouble()-0.5);
            b->_coords->_coords[1]+=1e-2*(_rand->generateDouble()-0.5);
            b->_coords->_coords[2]+=1e-2*(_rand->generateDouble()-0.5);
        }
    }

}
void Physics::MembraneFromObj::update(){


    updateTris();
    updateBeads(_updatable,0,0,5e-5);
    if(_step%1000==0){
        capture();
        std::cout<<_alpha<<"\t"<<_scale<<"\t"<<_step<<std::endl;

    }

    if(_step%10==0||_capture){

        int mRdis=0;
        _scale=(1-_alpha)*0.025+0.925;
        _alpha=std::fmin(1,_alpha+0.001);


        for(int i=0;i<_beadsDestI->size();i++){
            auto bi=_beadsDestI->at(i);
            auto bf=_beadsDestF->at(i);
            auto bt=_beadsDestTransiation->at(i);
            _temp->setValues(bi->_coords);
            _temp->multConst(1-_alpha);
            bt->_coords->setValues(_temp);
            _temp->setValues(bf->_coords);
            _temp->multConst(_alpha);
            bt->_coords->add(_temp);
        }

    }

    _step++;
}

Qt3DRender::QGeometryRenderer* Physics::MembraneFromObj::mesh(){
    return _disc->mesh();
}
void Physics::MembraneFromObj::capture(){
    SurfaceWithPhysics::capture();
    auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape_ScaledHD\Shape_%1_profile_%2.txt)").arg(*_shape,*_title);

    auto file=new QFile(s);


    if (!file->open(QIODevice::WriteOnly | QIODevice::Text)){
        assert(0)   ;
        return;
    }
    auto out=new QTextStream(file);
    writeBeadCoordinates(_xprofile,out);
    file->close();
    delete out;

}
void Physics::MembraneFromObj::updateVIN(){
    _disc->updateVertecies();
    _disc->updateIndecies();
    _disc->updateNormals();
}
double Physics::MembraneFromObj::calE(){
    return 0;
}
double Physics::MembraneFromObj::calE(int i){
    double ret=0;


    auto b=_beads->at(i);
    double l=b->_coords->len();
    if(l<_THRESHOLD+5){
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


    }else{
        ret+=10*b->_coords->_coords[1]*b->_coords->_coords[1];
    }
    ret-=_P*volumne();



    return ret;
}

void Physics::MembraneFromObj::saveObjToFile(QString *path){
    _disc->saveToFile(path);
}

