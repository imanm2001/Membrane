#include "membranefromobj.h"
#define _P 1
#define _K 50
#define _kappa 18522
#define _T 1
#define _E 2*_K/1.73205081
#define _THRESHOLD 42
#define _F 5010
#define _DT 6e-6

Physics::MembraneFromObj::MembraneFromObj(double dt):SurfaceWithPhysics(),_dt(dt),_appliedF(_F),_py(500),_frad(55)
{
    _pA=-1;
    //Amixed
    _FSign=-1;
    int res=_F%10;
    if(res>5){
        _radialForce=25*(5-res);
    }else{
        _radialForce=25*res;
    }

    _Rind=4;
    _cb=nullptr;
    _kappaFactor=1;
    _radiusFactor=1;
    _disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\disc_r44_d50_relaxed.obj)"));
    //_disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\disc_r44_d60_relaxed.obj)"));
    //_disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\oval3_r1_3scaled_smaller.obj)"));
    _tris=_disc->_tris;
    _scale=1e-5;
    _rand=_rand=QRandomGenerator::global();
    _border=new QVector<Geometry::BeadInfo*>();
    _xprofile=new QVector<Geometry::BeadInfo*>();
    _beads=_disc->_beads;
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
            int k=0;
            //            for(;k<_xprofile->size();k++){
            //                auto bb=_xprofile->at(k);
            //                if(bb->_coords->_coords[0]>b->_coords->_coords[0]){
            //                    break;
            //                }
            //            }
            //            _xprofile->insert(k,b);
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
    std::cout<<"TOTAL::"<<_initalArea<<"\t"<<_APV<<std::endl;
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
    auto bb=_xprofile->at(9);
    _kappaFactor=bb->_coords->_coords[0]/2;
    _kappaFactor=2*_P*_kappaFactor*_kappaFactor*_kappaFactor/_kappa;
    _radiusFactor=(bb->_coords->_coords[0]-1e-2)/_THRESHOLD;


    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        BendingEnergy *be=new BendingEnergy(b);
        be->orderConnections();
        be->updateBendingParameters();

    }

    std::cout<<":::"<<_TIE<<std::endl;
}

void Physics::MembraneFromObj::updateBeads(QVector<Geometry::BeadInfo*> *beads,double P,double kappa,double dt){
    double tA=0;
    double BE=0,SE=0;
    double NR=_THRESHOLD*_radiusFactor;
    bool thermal=0&(_step%20000)<100;
    double dts=0;

    if(thermal){
        std::sqrt(dt);
    }
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);

        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        double l=b->_coords->len();
        be->_curv=0.0;
        if( l< NR){
            //be->_curv=0.048*(1-std::exp(-(_step-0)/1000.0));
            be->_curv=2.0/(NR);

        }else{
            //  be->_curv=std::fmax(0.0,0.048*(1-(l-_THRESHOLD)/4));
            be->_curv=0.048;
        }
        be->_curv=0.0;
        be->orderConnections();
        be->updateBendingParameters();

    }

    double totalBE=0,KE=0;
    for(int i=0;_step>0&&i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        double l=b->_coords->len();
        _temp->zero();
        if(l<NR){


            //be->_curv=0.0;
            be->Ep2(-1,_temp);
            if(l>NR){
                _temp->multConst(-0*kappa);
            }else{
                _temp->multConst(-_kappa*_kappaFactor);
            }

            //    _temp->print();
            b->_force->add(_temp);
            BE+=_kappaFactor*_kappa*be->calE();
            //BE+=kappa*be->calcSig();
        }
    }



    int N=0;
    for(int i=0;i<_disc->_edges->size();i++){
        auto edge=_disc->_edges->at(i);
        edge->location(_temp);
        double r=_temp->len();
        double e=_sf->eval(edge);

        KE+=e;
        if(r<NR){

            BendingEnergy *be=(BendingEnergy *)edge->_vid1->getAttribute(BENDINGENERGY);
            double a=be->_Av;
            be=(BendingEnergy *)edge->_vid2->getAttribute(BENDINGENERGY);

            a+=be->_Av;

        }else{
            // edge->_vid1->_force->zero();
            //edge->_vid2->_force->zero();
        }



        //if(edge->_vid1->_coords->len()<_THRESHOLD||edge->_vid2->_coords->len()<_THRESHOLD){
        //
        //}
    }
    N=beads->size();

    double fixedBoundary=0;
    for(int i=0;i<beads->size();i++){
        auto b=beads->at(i);
        BendingEnergy *be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
        double l=b->_coords->len();



        if(l<NR+2){

            if(l<NR+1){
                _temp->setValues(be->_n);
                _temp->nomilize();
                _temp->multConst(-(_P)*be->_Av);
                b->_force->add(_temp);
                tA+=be->_Av;
            }


            //totalBE+=_kappa*be->calE();



        }
        if(l>NR-0.1){

            double f=-300000*b->_coords->_coords[1];
            double e=-f*b->_coords->_coords[1];
            fixedBoundary+=e;
            b->_force->_coords[1]+=f;
            //b->_coords->_coords[1]=b->_force->_coords[1]=0;
        }

    }




    if(_cb!=nullptr){
        //_cb->_force->_coords[1]+=std::fmin(_F,_appliedF);
        auto be=(BendingEnergy *)_cb->getAttribute(BENDINGENERGY);

        _cb->_force->_coords[1]+=_appliedF;
        _cb->_force->_coords[0]+=-1*_cb->_coords->_coords[0];
        _cb->_force->_coords[2]+=-1*_cb->_coords->_coords[2];

        //std::cout<<_APV/be->_Av<<"\t"<<std::endl;
        //_cb->_force->_coords[0]=0;
        //_cb->_force->_coords[2]=0;
    }
    //for(int i=0;i<_disc->_beads->size();i++){
    //  auto b=_disc->_beads->at(i);

    //_cfy=_cb->_force->_coords[1];
    _temp->zero();
    _temp2->zero();
    double ten=0,ten2=0,ten3=0;
    for(int i=0;i<beads->size();i++){
        auto b=beads->at(i);
        double r=b->_coords->len();
        if(b!=_cb&&r<NR){
            _temp2->setValues(b->_coords);
            _temp2->nomilize();
            auto be=(BendingEnergy *)b->getAttribute(BENDINGENERGY);
            _temp->setValues(be->_n);
            _temp->nomilize();
            double d=_temp->dot(_temp2);


            _temp->multConst(d);
            _temp2->sub(_temp);
            _temp2->nomilize();
            ten+=_temp2->dot(b->_force)/be->_Av;
            double dA=be->_Av-_APV;
            ten2+=(dA*std::fabs(dA))/_APV;
            ten3+=(dA*dA)/_APV;

        }
        b->update(dt);


#if _T>0
        if(thermal){

        b->_coords->_coords[0]+=0.75*(_rand->generateDouble()-0.5)*_T*dts;
        b->_coords->_coords[1]+=0.75*(_rand->generateDouble()-0.5)*_T*dts;
        b->_coords->_coords[2]+=0.75*(_rand->generateDouble()-0.5)*_T*dts;
        }
#endif
    }
    double minF=0;

    _temp->zero();
    for(int i=0;i<_border->size();i++){
        auto b=_border->at(i);
        double r=b->_coords->len();
        b->_coords->_coords[1]=b->_force->_coords[1]=0;

        if(r>NR+0.1){
            double f=b->_force->len();
            double dr=(_frad-r);
            minF=std::fmax(f,minF);
            //b->_force->_coords[2]+=10*dr*b->_coords->_coords[2]/r;
            //b->_force->_coords[0]+=10*dr*b->_coords->_coords[0]/r;

            // SE+=dr*dr*5;

            b->_force->_coords[0]+=_radialForce*b->_coords->_coords[0]/r;
            b->_force->_coords[2]+=_radialForce*b->_coords->_coords[2]/r;

            SE-=dr*_radialForce;
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
    //_tension=-(0.0*(-_TIE+BE)-0*totalBE+0*(_sf->_k/_K)*_E*(tA-_initalArea))/((_THRESHOLD*_radiusFactor)*(tA-_initalArea));

    _tension=(-_TIE+BE)*(0.04977764316091744)+(ten2)*(-0.001294612012995185)+(_F*_cb->_coords->_coords[1])*(0.0011541393111540513)+(volumne()*_P)*(0.016966224106566494)+(KE)*(0.024632178160253912)+(SE)*(-0.009050863839029583)+145.84154811364436;
    if(_step%100==0){
        double W=_F*_cb->_coords->_coords[1];
        double PV=volumne()*_P;
        std::cout<<-_TIE+BE<<"\t"<<ten2<<"\t"<<W<<"\t"<<PV<<"\t"<<KE<<"\t"<<SE<<"\t"<<fixedBoundary<<"\t"<<ten3<<"\t"<<tA<<std::endl;
        std::cout<<"DIFF"<<std::endl;
        std::cout<<BE-_pE[0]<<"\t"<<ten2-_pE[1]<<"\t"<<W-_pE[2]<<"\t"<<PV-_pE[3]<<"\t"<<KE-_pE[4]<<"\t"<<SE-_pE[5]<<"\t"<<fixedBoundary-_pE[6]<<"\t"<<ten3-_pE[7]<<"\t"<<tA-_pA<<std::endl;
        _pE[0]=BE;
        _pE[1]=ten2;
        _pE[2]=W;
        _pE[3]=PV;
        _pE[4]=KE;
        _pE[5]=SE;
        _pE[6]=fixedBoundary;
        _pE[7]=ten3;
        _pA=tA;

    }
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
        updateBeads(_updatable,p,kappa,_DT);
    }
    auto bb=_xprofile->at(_Rind);
    _radiusFactor=((bb->_coords->_coords[0]-5e-2)/_THRESHOLD);//*(1-std::exp(-_step/2000.0))+1;
    if(_step%1000==0||_capture){
        if(_step==1000){
            _ptension=_tension;
        }
        capture();
        double mRdis=0;
        for(int i=0;i<_xprofile->size();i++){
            auto b=_xprofile->at(i);
            double l=b->_coords->_coords[0];
            if(i==0||std::fabs(l-_THRESHOLD)<std::fabs(mRdis-_THRESHOLD)){
                _Rind=i;
                mRdis=l;
            }
        }
        std::cout<<_Rind<<std::endl;

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
            if(_step>1000&&std::fabs(_ptension2-_tension)<0.1){
                if(_tension>0&&_ptension<_tension){
                    _FSign*=-1;
                }
                if(_tension<0&&_ptension>_tension){
                    _FSign*=-1;
                }
                if(_tension>5){
                    // _radialForce+=0.1*_FSign;
                    //_frad*=0.999;
                }else if(_tension<0){
                    //  _radialForce+=0.5*_FSign;
                    //   _frad*=1.001;
                }
                _ptension=_tension;
            }
            _ptension2=_tension;



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
        std::cout<<_appliedF<<"\t"<<_sf->_k<<"\t"<<_initalArea<<"\t"<<_frad<<"\t"<<_border->at(0)->_coords->len()<<"\t"<< _kappaFactor<<"_"<<mRdis<<"\t"<<bb->_coords->_coords[0]<<"\t"<<_Rind<<"\t"<<_radialForce<<std::endl;
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

