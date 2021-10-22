#include "membranefromobj.h"
#define _P 0
#define _K 1000
#define _kappa 0
#define _T 0
#define _E 2*_K/1.73205081
#define _THRESHOLD 42
#define _F 0
Physics::MembraneFromObj::MembraneFromObj(double dt):SurfaceWithPhysics(),_dt(dt),_appliedF(_F),_py(500),_frad(55)
{
    //Mesh Project
    _Rind=8;
    _cb=nullptr;
    _kappaFactor=1;
    _radiusFactor=1;
    _disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\disc_r44_d60.obj)"));
    //_disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\oval3_r1_3scaled_smaller.obj)"));
    _tris=_disc->_tris;
    _scale=1e-5;
    _rand=_rand=QRandomGenerator::global();
    _border=new QVector<Geometry::BeadInfo*>();
    _xprofile=new QVector<Geometry::BeadInfo*>();
    _beads=_disc->_beads;
    _beadsSrc=new QVector<Geometry::BeadInfo*>();
    _updatable=new QVector<Geometry::BeadInfo*>();
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

        b->_coords->_coords[0]+=1e-4*(_rand->generateDouble()-0.5);
        b->_coords->_coords[1]+=1e-4*(_rand->generateDouble()-0.5);
        b->_coords->_coords[2]+=1e-4*(_rand->generateDouble()-0.5);

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
    _sf=new SpringForce(10000);
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
    loadFromFile();
    auto bb=_xprofile->at(8);
    _kappaFactor=bb->_coords->_coords[0]/2;
    _kappaFactor=2*_P*_kappaFactor*_kappaFactor*_kappaFactor/_kappa;
    _radiusFactor=(bb->_coords->_coords[0]-1e-2)/_THRESHOLD;
    _radiusFactor=1;
    std::cout<<"---"<<_kappaFactor<<"\t"<<_radiusFactor<<std::endl;
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        BendingEnergy *be=new BendingEnergy(b);
        be->orderConnections();
        be->updateBendingParameters();

    }

    std::cout<<_TIE<<std::endl;
}

void Physics::MembraneFromObj::updateBeads(QVector<Geometry::BeadInfo*> *beads,double P,double kappa,double dt){
    double tA=0;
    double BE=0,SE=0;


    double totalBE=0;



    int N=0;
    for(int i=0;i<_disc->_edges->size();i++){
        auto edge=_disc->_edges->at(i);
        edge->location(_temp);
        double r=_temp->len();
        double e=_sf->eval(edge);
        //SE+=e;


        //if(edge->_vid1->_coords->len()<_THRESHOLD||edge->_vid2->_coords->len()<_THRESHOLD){
        //
        //}
    }



    _temp->zero();
    _temp2->zero();
    double ten=0,ten2=0;
    for(int i=0;i<beads->size();i++){
        auto b=beads->at(i);
        b->_coords->_coords[1]=0;
        b->update(dt);


#if _T>0
        double dts=std::sqrt(dt);
        b->_coords->_coords[0]+=1e-1*(_rand->generateDouble()-0.5)*_T*dts;
        b->_coords->_coords[1]+=1e-1*(_rand->generateDouble()-0.5)*_T*dts;
        b->_coords->_coords[2]+=1e-1*(_rand->generateDouble()-0.5)*_T*dts;
#endif
    }
    double minF=0;
    if(_step%100==0){
        std::cout<<ten<<"\t"<<ten2<<std::endl;
    }
    _temp->zero();
    for(int i=0;i<_border->size();i++){
        auto b=_border->at(i);
        double r=b->_coords->len();
        b->_force->_coords[1]=0;
        if(r>42){
            double f=b->_force->len();
            minF=std::fmax(f,minF);
            //b->_force->_coords[2]+=10000*(_frad-r)*b->_coords->_coords[2]/r;
            //b->_force->_coords[0]+=10000*(_frad-r)*b->_coords->_coords[0]/r;

            b->_coords->_coords[0]=b->_coords->_coords[0]*_frad/r;
            b->_coords->_coords[1]=0;
            b->_coords->_coords[2]=b->_coords->_coords[2]*_frad/r;

            //;
            // if(f>100){
            //   if(std::fabs(_tension)>0.042){
            b->update(dt);
            //  }
        }
    }
    //   std::cout<<minF<<std::endl;
    //    std::cout<<tA<<std::endl;
    //    _tension=_E*(tA-296.397)/296.397;
    //_tension=_E*(tA-_initalArea)/_initalArea;
    //TENP
    //std::cout<<(_sf->_k/_K)*_E*(tA-_initalArea)<<"\t"<<BE<<"\t"<<volumne()*_P<<"\t"<<_F*_cb->_coords->_coords[1]<<"\t"<<N<<std::endl;

    if(_cb!=nullptr){

        totalBE=(-(-_F*_cb->_coords->_coords[1]+volumne()*_P));
    }
    // _tension=(_TIE-4*(_sf->_k/_K)*_E*(tA-_initalArea))/(2*_initalArea);
    _tension=-(0.5*(-_TIE+BE)-totalBE+0*(_sf->_k/_K)*_E*(tA-_initalArea))/((_THRESHOLD*_radiusFactor)*(tA-_initalArea));
    //  std::cout<<tA-_initalArea<<std::endl;
    // _tension=(totalBE-_pBE)/(tA-_initalArea);
    // std::cout<<_tension<<std::endl;
    // _tension=(totalBE-3*N)/(2*_initalArea);
    /*
    if(_tension>1e-4){
        _sf->_k/=10;
    }
    if(_tension<1e-5){
        _sf->_k*=10;
    }
*/
    _pBE=totalBE;
    //    _initalArea=tA;


}
void Physics::MembraneFromObj::update(){

    double p=_P;
    double kappa=_kappa;
    //_tet->_scale=0.2*std::exp(-_step/100.0)+0.9;
    _tet->_scale=1.1;
    //_tet->_k=_K*std::exp(-_step/1000.0)+0.01;
    _tet->_k =10000;
    // _sf->_k=100.0;
    /*  for(int i=0;i<_disc->_tris->size();i++){
        _disc->_tris->at(i)->getNormal();
    }*/
    for(int i=0;i<1;i++){
        updateTris();
        updateBeads(_updatable,p,kappa,1e-6);
    }
    auto bb=_xprofile->at(_Rind);
    _radiusFactor=((bb->_coords->_coords[0]-1e-2)/_THRESHOLD);//*(1-std::exp(-_step/2000.0))+1;
    if(_step%1000==0||_capture){
        capture();
        int mRdis=0;
        for(int i=0;i<_xprofile->size();i++){
            auto b=_xprofile->at(i);
            double l=b->_coords->_coords[0];
            if(i==0||std::fabs(l-_THRESHOLD)<std::fabs(mRdis-_THRESHOLD)){
                _Rind=i;
                mRdis=l;
            }
        }

        _kappaFactor=mRdis/2;

        _kappaFactor=2*_P*_kappaFactor*_kappaFactor*_kappaFactor/_kappa;

        if(_cb!=nullptr){
            if(std::fabs(_cb->_coords->_coords[1]-_py)<0.01 && _tension>1e-5){
                //_sf->_k*=0.99;
            }
            _cb->_coords->print();
            if(_cb->_coords->_coords[1]-_py<1&&_appliedF<_F){
                _appliedF+=10;
                _appliedF=std::fmin(_appliedF,_F);
            }
            if(_step>1000&&std::fabs((_tension-_ptension))<1e-3&&std::fabs(_cb->_coords->_coords[1]-_py)<1e-4){
                if(_tension>1){
                    _frad*=0.999;
                }else if(_tension<0){
                    //   _frad*=1.001;
                }
            }
            _ptension=_tension;


            if(_cb->_coords->_coords[1]<-1e-1&&_sf->_k>1&&_cb->_coords->_coords[1]-_py<0.1){
                //   _sf->_k-=10;
            }
            auto be=(BendingEnergy *)_cb->getAttribute(BENDINGENERGY);

        }
        /*
        if(_tension>1e-4){
            _sf->_k=std::fmin(500,_sf->_k/10);
        }
        if(_tension<1e-5){
            _sf->_k=std::fmax(5000,_sf->_k*10);
        }*/
        _py=_cb->_coords->_coords[1];
        std::cout<<"TEN:"<<_tension<<"\t"<<_initalArea<<std::endl;
        std::cout<<_appliedF<<"\t"<<_sf->_k<<"\t"<<_initalArea<<"\t"<<_frad<<"\t"<<_kappaFactor<<"\t"<<bb->_coords->_coords[0]<<"\t"<<_Rind<<std::endl;
        std::cout<<"STEP:"<<_step<<std::endl;
        std::cout.flush();
    }
    //std::exit(-1);
    /*
    for(int i=0;i<_beads->size();i++){
        auto b=_beads->at(i);
        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        be->orderConnections();
        be->updateBendingParameters();
    }*/

    //MC();

    _step++;
}

Qt3DRender::QGeometryRenderer* Physics::MembraneFromObj::mesh(){
    return _disc->mesh();
}
void Physics::MembraneFromObj::capture(){
    SurfaceWithPhysics::capture();
    auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape_Scaled\Shape_%1_profile_%2.txt)").arg(*_shape,*_title);

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


    for(int i=0;i<_disc->_edges->size();i++){
        //        ret+=_tet->calE(_disc->_edges->at(i));
    }
    for(int i=0;i<_tris->size();i++){
        auto tri=_tris->at(i);
        double min=100,max=-1;
        for(int j=0;j<3;j++){
            double l=tri->_v[j]->_coords->len();
            min=(min,l);
            max=fmax(max,l);
        }
        if(min>_THRESHOLD){
            double y=fabs(tri->getNormal()->_coords[1])-1;
            ret+=y*y*10;

        }
    }
    return ret;
}

void Physics::MembraneFromObj::saveObjToFile(QString *path){
    _disc->saveToFile(path);
}

