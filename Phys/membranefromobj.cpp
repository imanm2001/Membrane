#include "membranefromobj.h"
#define _P 0
#define _K 0.001
#define _kappa 1234
#define _T 0
#define _E 2*_K/1.73205081
#define _THRESHOLD 30
#define _F 000
Physics::MembraneFromObj::MembraneFromObj(double dt):SurfaceWithPhysics(),_dt(dt)
{

    _disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\disc_r40_d30.obj)"));
    _tris=_disc->_tris;
    _scale=1e-5;
    _rand=_rand=QRandomGenerator::global();
    _border=new QVector<Geometry::BeadInfo*>();
    _xprofile=new QVector<Geometry::BeadInfo*>();
    _beads=new QVector<Geometry::BeadInfo*>();
    SurfaceWithPhysics::_shape=new QString("Disc");

    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\hex.obj)"));
    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\cap.obj)"));
    //_sphere=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\pyramid.obj)"));
    for(int i=0;i<_disc->_tris->size();i++){
        auto tri=_disc->_tris->at(i);
        if(tri->getLocation()->len()>_THRESHOLD){
            tri->_curvater=0;
        }else{
            tri->_curvater=1;
        }
    }
    _temp=new Physics::VecD3d();
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        BendingEnergy *be=new BendingEnergy(b);
        be->orderConnections();

        if(fabs(b->_coords->_coords[2])<1e-3){
            _xprofile->append(b);
            if(fabs(b->_coords->_coords[0])<1e-3){
                _cb=b;
            }
        }

        b->_coords->_coords[0]+=1e-1*(_rand->generateDouble()-0.5);
        b->_coords->_coords[1]+=1e-1*(_rand->generateDouble()-0.5);
        b->_coords->_coords[2]+=1e-1*(_rand->generateDouble()-0.5);

    }
    _initalArea=0;
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        BendingEnergy *be=new BendingEnergy(b);

        be->orderConnections();
        be->updateBendingParameters();

        double r=b->_coords->len();

        if(b->getAttribute(BORDER)!=nullptr){
            _border->append(b);
            be->_border=true;
        }else {
            if(r<_THRESHOLD){
                _initalArea+=be->_Av;
                double d=r/_THRESHOLD;
                //b->_coords->_coords[1]=10*(1-d*d);
            }

            _beads->append(b);
        }


    }
    std::cout<<"TOTAL::"<<_initalArea<<std::endl;
    _sf=new SpringForce(_K);
    _tet=new Tether(_K);
    _step=0;

    SurfaceWithPhysics::_title=new QString(QString("%1_%2_%3_%4").arg(QString::number(_P),QString::number(_K),QString::number(_kappa),QString::number(_F)));
    loadFromFile();
    for(int i=0;i<_border->size();i++){
        auto b=_border->at(i);
        //b->_coords->multConst(10.0/b->_coords->len());
    }
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        BendingEnergy *be=new BendingEnergy(b);
        if(b->_coords->len()<_THRESHOLD){
            be->_curv=1;
        }else{
            be->_curv=0;
        }
        be->orderConnections();
        be->updateBendingParameters();
    }
    _pBE=0;
    allocMem();
    std::cout<<volumne()<<std::endl;
}

void Physics::MembraneFromObj::updateBeads(QVector<Geometry::BeadInfo*> *beads,double P,double kappa,double dt){
    double tA=0;
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        double l=b->_coords->len();
        if(l<_THRESHOLD+2){
            _temp->zero();
            if( l<_THRESHOLD){
                be->_curv=0.067*(1-std::exp((_step-0)/100.0));
                be->_curv=0.033;

            }else{
                be->_curv=0.0;
            }
            be->Ep2(-1,_temp);
            _temp->multConst(-kappa);

            //    _temp->print();
            b->_force->add(_temp);
        }

        //be->orderConnections();
        be->updateBendingParameters();
    }
    double totalBE=0;
    int N=0;

    for(int i=0;i<beads->size();i++){
        auto b=beads->at(i);
        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);




        if(b->_coords->len()<_THRESHOLD){
            _temp->setValues(be->_n);
            _temp->multConst((-P+1)*be->_Av);
            b->_force->add(_temp);
            tA+=be->_Av;
            N++;
            totalBE+=_kappa*be->calE();
        }


    }



    for(int i=0;i<_disc->_edges->size();i++){
        auto edge=_disc->_edges->at(i);
        edge->location(_temp);
        double r=_temp->len();

        if(r<_THRESHOLD){
            BendingEnergy *be=(BendingEnergy *)edge->_vid1->getAttribute(BENDINGENERGY);
            double a=be->_Av;
            be=(BendingEnergy *)edge->_vid2->getAttribute(BENDINGENERGY);
            a+=be->_Av;
            totalBE+=_tet->calE(edge);
        }

        _tet->eval(edge);

        //double e=_sf->eval(edge);
        //if(edge->_vid1->_coords->len()<_THRESHOLD||edge->_vid2->_coords->len()<_THRESHOLD){
        //  totalBE+=e;
        //}
    }

    _cb->_force->_coords[1]+=_F+1000;
    for(int i=0;i<beads->size();i++){
        auto b=beads->at(i);
        double r=b->_coords->len();
        if(r>_THRESHOLD){

            b->_force->_coords[1]=-0.1*b->_coords->_coords[1];
            b->_force->_coords[1]=0;

        }else{

        }
        b->update(dt);

#if _T>0
        double dts=std::sqrt(dt);
        b->_coords->_coords[0]+=1e0*(_rand->generateDouble()-0.5)*_T*dts;
        b->_coords->_coords[1]+=1e0*(_rand->generateDouble()-0.5)*_T*dts;
        b->_coords->_coords[2]+=1e0*(_rand->generateDouble()-0.5)*_T*dts;
#endif
    }
    /*
    for(int i=0;i<_border->size();i++){
        auto b=_border->at(i);
        double r=b->_coords->len();
        b->_force->_coords[1]=0;
        b->_force->_coords[2]=1*(100-r)*b->_coords->_coords[2]/r;
        b->_force->_coords[0]=1*(100-r)*b->_coords->_coords[0]/r;


        b->update(dt);

    }*/
    //    std::cout<<tA<<std::endl;
    //    _tension=_E*(tA-296.397)/296.397;
    //_tension=_E*(tA-_initalArea)/_initalArea;
    totalBE-=3*volumne()*_P;
    totalBE-=_F*_cb->_coords->_coords[1]*1;
    _tension=(totalBE-_pBE)/(tA-_initalArea);
    _tension=(totalBE)/(2*_initalArea);
    _pBE=totalBE;
    //_initalArea=tA;


}
void Physics::MembraneFromObj::update(){

    double p=_P;
    double kappa=_kappa;
    //_sf->_k=_K*(1-std::exp(-_step/2000.0));
    _tet->_k =_K;
    updateTris();
    for(int i=0;i<5;i++){

        updateBeads(_beads,p,kappa,1e-7);
    }
    if(_step%1000==0||_capture){
            capture();
        _cb->_coords->print();
        std::cout<<_tension<<std::endl;
        std::cout<<volumne()<<std::endl;
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
    auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_%1_profile_%2.txt)").arg(*_shape,*_title);

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
    if(l<_THRESHOLD+1){
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
        ret+=_tet->calE(_disc->_edges->at(i));
    }
    for(int i=0;i<_tris->size();i++){
        auto tri=_tris->at(i);
        double min=100,max=-1;
        for(int j=0;j<3;j++){
            double l=tri->_v[j]->_coords->len();
            min=fmin(min,l);
            max=fmax(max,l);
        }
        if(min>_THRESHOLD){
            double y=fabs(tri->getNormal()->_coords[1])-1;
            ret+=y*y*10;

        }
    }
    return ret;
}

