#include "membranefromobj.h"

#define _P 1
#define _K 50006
#define _kappa 18522
//#define _kappa 1
#define _T 1
#define _E 2*_K/1.73205081
#define _THRESHOLD 42
#define _F 437

#define _DT 3e-6

Physics::MembraneFromObj::MembraneFromObj(double dt):SurfaceWithPhysics(),_dt(dt),_appliedF(0),_frad(55)
{
    _dtF=1;
    _tempTri=new Geometry::Triangle(this,12345,new Geometry::BeadInfo(this,new VecD3d(1,1,0),1,0),new Geometry::BeadInfo(this,new VecD3d(1,0,0),1,1),new Geometry::BeadInfo(this,new VecD3d(0,1,0),1,2),0);
    _cts=new CTS(2);

    _CTSvs1=(VecD3d**)malloc(sizeof(VecD3d*)*4);
    _CTSvs2=(VecD3d**)malloc(sizeof(VecD3d*)*4);

    for(int i=0;i<4;i++){
        _CTSvs1[i]=new Physics::VecD3d();
        _CTSvs2[i]=new Physics::VecD3d();
    }
    _CTSbeads=(Geometry::BeadInfo**)malloc(sizeof(BeadInfo*)*4);
    _strainMatrix=gsl_matrix_alloc(2,2);
    _strainDirection=gsl_vector_alloc(2);
    _strainDirection2=gsl_vector_alloc(2);

    _pA=-1;
    //Amixed
    _FSign=-1;
    int res=_F%10;
    if(res>5){
        _radialForce=25*(5-res);
    }else{
        _radialForce=25*res;
    }

    _appliedF=0;
    _radialForce=-0;
    _step=000;
    _saveData[0]=_appliedF;
    _saveData[1]=_radialForce;
    _saveData[2]=_step;

    _Rind=4;
    _cb=nullptr;
    _kappaFactor=1;
    _radiusFactor=1;
    //_disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\sphere2_r1.obj)"));
    //_disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\disc_r55_d60_relaxed.obj)"));
    _disc=new Geometry::WaveFrontObj(QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\disc_r80_d90_relaxed.obj)"));
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

    _temp=new Physics::VecD3d();
    _temp2=new Physics::VecD3d();
    _temp3=new Physics::VecD3d();
    double mRdis=0;
    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        if(fabs(b->_coords->_coords[2])<1e-3){
            double l=b->_coords->len();
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



    }

    for(int i=0;i<_disc->_tris->size();i++){
        auto tri=_disc->_tris->at(i);
        for(int j=0;j<3;j++){
            tri->_v[j]->_coords->_coords[1]=0;
        }
        tri->setOrginals();

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
        //be->_curv=0.048;
        be->_curv=1;
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
    _sf=new SpringForce(_K);

    _tet=new Tether(100000);


    SurfaceWithPhysics::_title=new QString(QString("%1_%2_%3_%4").arg(QString::number(_P),QString::number(_K),QString::number(_kappa),QString::number(_F)));
    _frad=1;
    double R=0;
    for(int i=0;i<_border->size();i++){
        auto b=_border->at(i);
        double l=b->_coords->len();
        R+=l;
        if(i==0||l>_frad){
            _frad=l;
        }
        //b->_coords->multConst(10.0/b->_coords->len());
    }
    R/=_border->size();


    _pBE=0;
    allocMem();
    std::cout<<volumne()<<"R:\t"<<R<<std::endl;
    double r=0;
    for(int i=0;i<_beads->size();i++){
        auto b=_beads->at(i);
        r+=b->_coords->len();

    }
    std::cout<<r/_beads->size()<<std::endl;
    loadFromFile(_saveData,_NUMSAVEDATA);
    _appliedF=_saveData[0];
    _radialForce=_saveData[1];
    _pstep=_step=_saveData[2];

    if(_xprofile->size()>0){
        auto bb=_xprofile->at(9);

        _kappaFactor=bb->_coords->_coords[0]/2;
        _kappaFactor=2*_P*_kappaFactor*_kappaFactor*_kappaFactor/_kappa;
        _radiusFactor=(bb->_coords->_coords[0]-1e-2)/_THRESHOLD;
    }else{
        _kappaFactor=1;
        _radiusFactor=1;
    }

    for(int i=0;i<_disc->_beads->size();i++){
        auto b=_disc->_beads->at(i);
        b->_coords->_coords[0]+=5e-3*(_rand->generateDouble()-0.5);
        b->_coords->_coords[1]+=5e-3*(_rand->generateDouble()-0.5);
        b->_coords->_coords[2]+=5e-3*(_rand->generateDouble()-0.5);

        BendingEnergy *be=new BendingEnergy(b);
        be->orderConnections();
        be->updateBendingParameters();

    }
    if(_cb!=nullptr){
        _py=_cb->_coords->_coords[1];
    }
    std::cout<<":::"<<_TIE<<std::endl;
    _RFA=_APFA=0;
}
void Physics::MembraneFromObj::findThefourthVertex(Geometry::Triangle* t1,Geometry::Triangle* t2,Geometry::BeadInfo** bi){
    std::copy(t1->_v,t1->_v+3,bi);
    bi[3]=nullptr;
    bool found=false;
    for(int i=0;i<3&&!found;i++){

        auto b=t2->_v[i];
        for(int j=0;!found&&j<3;j++){
            found=b==t1->_v[j];
        }
        if(!found){
            bi[3]=b;

        }else{
            found=false;
        }
    }
    assert(bi[3]!=nullptr);
}
double Physics::MembraneFromObj::calStrain(){
    double ret=0,area=0;
    /*
    for(int i=0;i<_beads->size();i++){
        r+=_beads->at(i)->_coords->len();
    }*/
    double NR=_THRESHOLD*_radiusFactor;
    for(int i=0;i<_disc->_tris->size();i++){
        auto tri=_disc->_tris->at(i);

        if(tri->getLocation()->len()<NR){

            for(int j=0;j<3;j++){
                auto tri2=tri->_e[j]->_tris[0]!=tri?tri->_e[j]->_tris[0]:tri->_e[j]->_tris[1];
                if(tri2!=nullptr&&tri2!=tri){
                    findThefourthVertex(tri,tri2,_CTSbeads);
                    break;
                }
            }
            for(int j=0;j<4;j++){
                _CTSvs1[j]=_CTSbeads[j]->_originalLocations;
                _CTSvs2[j]=_CTSbeads[j]->_coords;


            }

            if(_cts->getStrainMatrix2D(_CTSvs1,_CTSvs2,_strainMatrix)){
                _temp->setValues(tri->getNormal());
                //        _temp->nomilize();
                //get Phi Vector
                _temp2->setValues(-_temp->_coords[2],0,_temp->_coords[0]);
                _temp2->nomilize();
                _temp2->cross(tri->getNormal(),_temp);
                _temp2->nomilize();
                _cts->setVecToVec2D(_temp,_strainDirection);
                gsl_blas_dgemv(CblasNoTrans,1,_strainMatrix,_strainDirection,0,_strainDirection2);
                double strainT;
                gsl_blas_ddot(_strainDirection,_strainDirection2,&strainT);
                double a=tri->area();
                ret+=strainT*a;
                /*
                _cts->setVecToVec(_temp2,_strainDirection);
                gsl_blas_dgemv(CblasNoTrans,1,_strainMatrix,_strainDirection,0,_strainDirection2);

                gsl_blas_ddot(_strainDirection,_strainDirection2,&strainT);
                ret+=strainT*a;*/
                /*
        _temp2->setValues(tri->getLocation());
        _temp2->nomilize();
        _cts->setVecToVec(_temp2,_strainDirection);
        gsl_blas_dgemv(CblasNoTrans,1,_strainMatrix,_strainDirection,0,_strainDirection2);

        gsl_blas_ddot(_strainDirection,_strainDirection2,&strainT);
        ret+=strainT*a;*/

                area+=a;
            }
        }

    }
    std::cout<<area<<std::endl;
    return ret*area;
}
void Physics::MembraneFromObj::Project2D(VecD3d** res,VecD3d* e1,VecD3d* e2,VecD3d* temp){

    double edots[6]={0,0,0,0,0,0};
    for(int i=0;i<6;i++){
        edots[i]=0;
    }
    double ddot=fabs(e1->dot(e2));
    if(ddot!=ddot||ddot>=1e-4||e1->len()<0.9||e2->len()<0.9){
        std::cout<<"---->>><<EE\t"<<e1->dot(e2)<<std::endl;
        e1->print();
        e2->print();

        assert(0);
    }
    temp->zero();
    for(int i=0;i<3;i++){
        temp->add(res[i]);
    }
    temp->multConst(1/3.0);
    for(int i=0;i<3;i++){
        auto v=res[i];
        v->sub(temp);


        edots[i*2]=e1->dot(v);
        edots[i*2+1]=e2->dot(v);

    }

    for(int i=0;i<3;i++){
        res[i]->setValues(edots[i*2],0,edots[i*2+1]);

    }
    /*
    temp->zero();
    for(int i=0;i<3;i++){
        temp->add(res[i]);
    }
    temp->multConst(1/3.0);*/
    /*for(int i=0;i<3;i++){
       res[i]->sub(temp);
    }*/

}
double Physics::MembraneFromObj::calStrain2D2(){
    double ret=0,area=0,r=0;

    for(int i=0;i<_beads->size();i++){
        r+=_beads->at(i)->_coords->len();
    }
    double NR=_THRESHOLD*_radiusFactor;
    for(int n=0;n<_disc->_tris->size();n++){
        auto tri=_disc->_tris->at(n);
        _temp->zero();
        _temp2->zero();
        _temp3->zero();
        if(tri->getLocation()->len()<NR){
            for(int j=0;j<3;j++){

                _CTSvs2[j]->setValues(tri->_v[j]->_coords);
            }

            double dr=0;
            for(int j=0;j<3;j++){
                _temp->setValues(tri->_e[j]->_vid1->_coords);
                _temp2->setValues(tri->_e[j]->_vid2->_coords);
                _temp2->sub(_temp);
                dr+=(_temp2->len()-tri->_e[j]->_restLength)/tri->_e[j]->_restLength;
            }
            dr/=3;
            _temp->setValues(tri->getNormal());
            if(fabs(_temp->_coords[1])==1){

                _temp->nomilize();
                bool b=_temp->_coords[1]!=_temp->_coords[1]||fabs(_temp->_coords[1])==1;

                for(int i=0;i<3&&b==1;i++){
                    _temp->setValues(tri->_v[i]->_coords);
                    _temp->nomilize();
                    b=_temp->_coords[1]!=_temp->_coords[1]||fabs(_temp->_coords[1])==1;
                }

            }
            //temp3=normal
            _temp3->setValues(_temp);
            _temp3->nomilize();
            _temp->setValues(tri->getLocation());
            //        _temp->nomilize();
            //temp2= Phi Vector


            _temp2->setValues(_temp->_coords[2],0,-_temp->_coords[0]);

            _temp2->nomilize();

            //_temp= theta
            _temp2->cross(_temp3,_temp);

            _temp->nomilize();


            Project2D(_CTSvs2,_temp,_temp2,_temp3);

            double mY=_CTSvs2[0]->_coords[0];
            double MY=mY;
            for(int j=1;j<2;j++){
                mY=fmin(mY,_CTSvs2[j]->_coords[0]);
                MY=fmax(MY,_CTSvs2[j]->_coords[0]);
            }
            ret+=(MY-mY)*dr;
        }

    }


    std::cout<<area<<"\t"<<r/(double)_beads->size()<<std::endl<<std::endl<<"---2:"<<std::endl;
    return ret;
}
void Physics::MembraneFromObj::calRect(VecD3d** vis,double* w,double* h){
    double mY=vis[0]->_coords[0];
    double MY=mY;

    double mX=vis[0]->_coords[2];
    double MX=mX;
    for(int j=1;j<3;j++){
        mY=fmin(mY,vis[j]->_coords[0]);
        MY=fmax(MY,vis[j]->_coords[0]);

        mX=fmin(mX,vis[j]->_coords[2]);
        MX=fmax(MX,vis[j]->_coords[2]);
    }
    h[0]=MY-mY;
    w[0]=MX-mX;
}
double Physics::MembraneFromObj::calStrain2D(){
    double ret=0,area=0,r=0;

    for(int i=0;i<_border->size();i++){
        r+=_border->at(i)->_coords->len();
    }
    r/=_border->size();
    double NR=_THRESHOLD*_radiusFactor;
    for(int n=0;n<_disc->_tris->size();n++){
        auto tri=_disc->_tris->at(n);
        bool b=tri->getLocation()->len()<NR+0.8;
        if(b){
            _temp->zero();
            _temp2->zero();
            _temp3->zero();

            for(int j=0;j<3;j++){
                _CTSvs1[j]->setValues(tri->_v[j]->_originalLocations);
                _CTSvs2[j]->setValues(tri->_v[j]->_coords);

            }


            if(true){
                _temp->setValues(tri->_oLocation);
                _temp->_coords[1]=0;
                bool bb=_temp->len()<1e-5;

                for(int i=0;i<3&&bb==1;i++){
                    _temp->setValues(tri->_v[i]->_originalLocations);
                    _temp->_coords[1]=0;
                     bb=_temp->len()<1e-5;
                    _temp->nomilize();
                }
                _temp->_coords[1]=0;
                _temp->nomilize();
                //temp2=phi
                _temp2->_coords[0]=_temp->_coords[2];
                _temp2->_coords[1]=0;
                _temp2->_coords[2]=-_temp->_coords[0];
                _temp2->nomilize();

            }else{
                _temp->setValues(tri->_oNorm);

                if(fabs(_temp->_coords[1])==1){
                    _temp->setValues(tri->_oLocation);
                    _temp->nomilize();
                    bool b=_temp->_coords[1]!=_temp->_coords[1]||fabs(_temp->_coords[1])==1;
                    for(int i=0;i<3&&b;i++){
                        _temp->setValues(tri->_v[i]->_originalLocations);
                        _temp->nomilize();
                        b=_temp->_coords[1]!=_temp->_coords[1]||fabs(_temp->_coords[1])==1;
                    }

                }

                //Normal Direction
                _temp->nomilize();
                _temp3->setValues(_temp);

                //temp2= Phi Vector
                _temp2->setValues(_temp3->_coords[2],0,-_temp3->_coords[0]);
                _temp2->nomilize();

                //_temp= theta
                _temp2->cross(_temp3,_temp);
                _temp->nomilize();

            }
            _temp->nomilize();
            _temp2->nomilize();
            /*
            _temp->setValues(1,0,0);
            _temp2->setValues(0,0,1);*/
            //_temp=r
            Project2D(_CTSvs1,_temp,_temp2,_temp3);

            _temp->setValues(tri->getNormal());
            if(fabs(_temp->_coords[1])==1){

                _temp->nomilize();
                bool b=_temp->_coords[1]!=_temp->_coords[1]||((1-fabs(_temp->_coords[1]))<1e-6);

                for(int i=0;i<3&&b==1;i++){
                    _temp->setValues(tri->_v[i]->_coords);
                    _temp->nomilize();
                    b=_temp->_coords[1]!=_temp->_coords[1]||((1-fabs(_temp->_coords[1]))<1e-6);
                }

            }
            //temp3=normal
            _temp3->setValues(_temp);
            _temp3->nomilize();
            _temp->setValues(tri->getLocation());
            //        _temp->nomilize();
            //temp2= Phi Vector


            _temp2->setValues(_temp->_coords[2],0,-_temp->_coords[0]);

            _temp2->nomilize();

            //_temp= theta
            _temp2->cross(_temp3,_temp);

            _temp->nomilize();


            Project2D(_CTSvs2,_temp,_temp2,_temp3);
            if(_cts->getStrainMatrix2D(_CTSvs1,_CTSvs2,_strainMatrix)){
                _temp2->setValues(1,0,0);
                _cts->setVecToVec2D(_temp2,_strainDirection);
                _strainDirection2->data[0]=_strainDirection2->data[1]=0;
                gsl_blas_dgemv(CblasNoTrans,1,_strainMatrix,_strainDirection,0,_strainDirection2);
                double strainT=0;

                gsl_blas_ddot(_strainDirection,_strainDirection2,&strainT);
                double a=tri->_area;

/*
                _temp2->setValues(0,0,1);
                _temp2->nomilize();
                _cts->setVecToVec2D(_temp2,_strainDirection);
                gsl_blas_dgemv(CblasNoTrans,1,_strainMatrix,_strainDirection,0,_strainDirection2);
                double strainT2=0;

                gsl_blas_ddot(_strainDirection,_strainDirection2,&strainT2);

                //ret+=strainT/strainT;
                double dA=a-tri->_oA;
                //ret+=(dA-strainT2)/tri->_oA;
                double w1=0,w2=0,h1=0,h2=0;
                calRect(_CTSvs1,&w1,&h1);
                calRect(_CTSvs2,&w2,&h2);

                //ret+=strainT;

                if(h1<1e-2){
                    for(int i2=0;i2<3;i2++){
                        std::cout<<"---"<<std::endl;
                        _CTSvs1[i2]->print();
                        tri->_v[i2]->_originalLocations->print();

                    }
                    std::cout<<_tempTri->area()<<"<<\t>>"<<tri->_oA<<std::endl;

                    std::cout<<":::AA   "<<w1<<"\t"<<h1<<"\t"<<a<<std::endl;
                    tri->_oNorm->print();
                    ret+=strainT;
                    assert(0);
                }
*/
                ret+=strainT*a*(b?1:-1);

                /*
                _cts->setVecToVec(_temp2,_strainDirection);
                gsl_blas_dgemv(CblasNoTrans,1,_strainMatrix,_strainDirection,0,_strainDirection2);

                gsl_blas_ddot(_strainDirection,_strainDirection2,&strainT);
                ret+=strainT*a;*/
                /*
        _temp2->setValues(tri->getLocation());
        _temp2->nomilize();
        _cts->setVecToVec(_temp2,_strainDirection);
        gsl_blas_dgemv(CblasNoTrans,1,_strainMatrix,_strainDirection,0,_strainDirection2);

        gsl_blas_ddot(_strainDirection,_strainDirection2,&strainT);
        ret+=strainT*a;*/

                area+=a;
            }
        }

    }


  //  std::cout<<area<<"\tR:\t"<<r<<std::endl<<std::endl<<"---:"<<std::endl;
    _maxR=r;
    return ret;
}

void Physics::MembraneFromObj::testStrain(){
    std::cout<<":::CTS"<<std::endl;
    Geometry::BeadInfo **beads;
    beads=(BeadInfo**)malloc(sizeof(BeadInfo*)*4);
    VecD3d** vs=(VecD3d**)malloc(sizeof(VecD3d*)*4);
    VecD3d** vs2=(VecD3d**)malloc(sizeof(VecD3d*)*4);

    vs[0]=new VecD3d(0,1e-2,0);
    vs[1]=new VecD3d(1,4e-3,0);
    vs[2]=new VecD3d(0,5e-3,1);
    vs[3]=new VecD3d(1,1e-2,1);


    vs2[0]=new VecD3d(0,1e-2,0);
    vs2[1]=new VecD3d(1,4e-3,0);
    vs2[2]=new VecD3d(0,5e-3,2);
    vs2[3]=new VecD3d(1,1e-2,2);

    /*
    vs[0]=new VecD3d(8.66025404e-01, -5.00000000e-01,1e-2);
    vs[1]=new VecD3d(-8.66025404e-01, -5.00000000e-01,4e-3);
    vs[2]=new VecD3d(6.12323400e-17,  1.00000000e+00,5e-3);
    vs[3]=new VecD3d(1.73205081e+00, 1.00000000e+00,1e-2);


    vs2[0]=new VecD3d(0.96592583,0.25881905,1e-2);
    vs2[1]=new VecD3d(-0.25881905,-0.96592583,4e-3);
    vs2[2]=new VecD3d(-0.70710678,0.70710678,5e-3);
    vs2[3]=new VecD3d(0.51763809,1.93185165,1e-2);*/

    beads[0]=new Geometry::BeadInfo(this,vs2[0],0,1);
    beads[1]=new BeadInfo(this,vs2[1],0,1);
    beads[2]=new BeadInfo(this,vs2[2],0,2);
    beads[3]=new BeadInfo(this,vs2[3],0,3);

    Geometry::Triangle *t1,*t2;
    t1=new Geometry::Triangle(this,0,beads[0],beads[1],beads[2],0);
    t2=new Geometry::Triangle(this,0,beads[1],beads[3],beads[2],0);

    gsl_matrix *m=gsl_matrix_alloc(2,2);

    _cts->getStrainMatrix2D(vs,vs2,m);
    std::cout<<std::endl;
    std::cout<<std::endl;
    _cts->printM(m);

    std::cout<<std::endl;
    double nd[3]={0,1};

    double nd2[3]={1,0};
    //double nd[3]={1,0,0};

    auto n=gsl_vector_view_array(nd,2);
    auto n2=gsl_vector_view_array(nd2,2);

    auto temp=gsl_vector_alloc(2);
    gsl_blas_dgemv(CblasNoTrans,1,m,&n.vector,0,temp);
    double res=123;
    gsl_blas_ddot(&n.vector,temp,&res);
    std::cout<<"----"<<std::endl;
    gsl_vector_fprintf(stdout,temp,"%g");
    std::cout<<res<<std::endl;




}

void Physics::MembraneFromObj:: updateBeads(QVector<Geometry::BeadInfo*> *beads,double P,double kappa,double dt){
    double tA=0;
    double BE=0,SE=0;
    double NR=_THRESHOLD*_radiusFactor;
    int t=(_step%1000);
    bool thermal=t>100&&t<800;;
    thermal=1;
    double dts=0;

    if(thermal){
        dts=std::sqrt(dt);
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
        //be->_curv= 1;
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
                //_temp->multConst(-_kappa);
            }

            //    _temp->print();
            b->_force->add(_temp);
            //BE+=_kappaFactor*_kappa*be->calE();
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
                _temp->multConst((-_P)*be->_Av);
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
        if(b->_coords->_coords[1]<0){
            b->_force->_coords[1]-=10000*b->_coords->_coords[1];

        }



#if _T>0
        if(thermal){

            b->_coords->_coords[0]+=1e0*(_rand->generateDouble()-0.5)*_T*dts;
            b->_coords->_coords[1]+=1e0*(_rand->generateDouble()-0.5)*_T*dts;
            b->_coords->_coords[2]+=1e0*(_rand->generateDouble()-0.5)*_T*dts;
        }

#endif
          b->update(dt);
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

    //_tension=(-_TIE+BE)*(0.04977764316091744)+(ten2)*(-0.001294612012995185)+(_F*_cb->_coords->_coords[1])*(0.0011541393111540513)+(volumne()*_P)*(0.016966224106566494)+(KE)*(0.024632178160253912)+(SE)*(-0.009050863839029583)+145.84154811364436;
    if(_step%1000==0){
        std::cout<<NR<<"\t"<<tA<<"\t"<<_sf->_k<<std::endl;
    }
    if(_step%100==0&&_cb!=nullptr){
        double W=_F*_cb->_coords->_coords[1];
        double PV=volumne()*_P;
        /*
        std::cout<<-_TIE+BE<<"\t"<<ten2<<"\t"<<W<<"\t"<<PV<<"\t"<<KE<<"\t"<<SE<<"\t"<<fixedBoundary<<"\t"<<ten3<<"\t"<<tA<<std::endl;
        std::cout<<"DIFF"<<std::endl;
        std::cout<<BE-_pE[0]<<"\t"<<ten2-_pE[1]<<"\t"<<W-_pE[2]<<"\t"<<PV-_pE[3]<<"\t"<<KE-_pE[4]<<"\t"<<SE-_pE[5]<<"\t"<<fixedBoundary-_pE[6]<<"\t"<<ten3-_pE[7]<<"\t"<<tA-_pA<<std::endl;*/
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
    _appliedF=fmin(fmax(_F, _appliedF+_APFA),800);
    _radialForce=fmax(-4000,fmin(_radialForce+_RFA,0));

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
    Geometry::BeadInfo *bb=nullptr;
    if(_Rind>0&&_Rind<_xprofile->size()){
        auto bb=_xprofile->at(_Rind);
        _radiusFactor=((bb->_coords->len()+1e-3)/_THRESHOLD);//*(1-std::exp(-_step/2000.0))+1;
    }else{
        _radiusFactor=1;
    }

    if(_step%1000==0||_capture){
        double cR=-1;
        for(int i=0;i<_xprofile->size();i++){
            auto b=_xprofile->at(i);

            if(b->_coords->_coords[1]<0){
                double l=b->_coords->len();
                if (l<cR||cR==-1){
                    cR=l;
                }
            }
        }
        if(_step==1000){
            _ptension=_tension;
        }
        _saveData[0]=_appliedF;
        _saveData[1]=_radialForce;
        _saveData[2]=_step;
        capture(_saveData,_NUMSAVEDATA);
        double mRdis=-1;
        for(int i=0;i<_xprofile->size();i++){
            auto b=_xprofile->at(i);
            double l=b->_coords->len();
            if(i==0||std::fabs(l-_THRESHOLD)<std::fabs(mRdis-_THRESHOLD)){
                _Rind=i;
                mRdis=l;
            }
        }
        if(cR==-1){
            cR=mRdis;
        }
        std::cout<<_Rind<<std::endl;

        _kappaFactor=mRdis/2;

        _kappaFactor=2*_P*_kappaFactor*_kappaFactor*_kappaFactor/_kappa;
        double strain=calStrain2D();
        double r=_THRESHOLD*_radiusFactor;
        //double param=(_maxR*_maxR-54.3153*54.3153)-(r*r-42.064*42.064);
        double param=(_maxR*_maxR-79.6341*79.6341)-(r*r-42.0564*42.0564);
        double alpha=cR/(mRdis-0.8);
        std::cout<<std::endl<<alpha<<"\t"<<cR<<"\t"<<mRdis<<std::endl;
        //_tension=(strain*(14.607392476665238)+param*(2.611268641852894)+(464.58341269938956))+0*27086.6*(1-alpha);
        //_tension=strain*(14.5093297530095)+param*(3.25999986286536)+(521.642529910544); //r55  d60
        //_tension=strain*(14.125051896565202)+param*(1.3863713270604137)+(66.91552217921635);
        _tension=strain*(14.965913222422378)+param*(0.757068516982182)+(65.57094366299579);
        if(_cb!=nullptr){

            _cb->_coords->print();

            if(_step>100){

                if(1||_tension>20){
                   // _radialForce-=((_tension-_ptension)*0.1+(2e-2)*(_tension+800))*(_step-_pstep)*_DT*1e4;
                   // _RFA=-((_tension-_ptension)*0.1+(5e-2)*_tension)*(_step-_pstep)*_DT*1e-1;
                    //_frad*=0.999;
                }
                else if(_tension>10){
                    _radialForce-=10;
                }
                else if(_tension>5){
                    _radialForce-=1;
                }
                else if(_tension<-20){
                    _radialForce-=((_tension-_ptension)*0.1+(2e-2)*_tension)*(_step-_pstep)*1e-3;
                    //   _frad*=1.001;
                }
                else if(_tension<-10){
                    _radialForce+=5;
                    //   _frad*=1.001;
                }
                else if(_tension<-5){
                    _radialForce+=0.5;
                    //   _frad*=1.001;
                }
                _ptension=_tension;
                //_radialForce=fmax(_radialForce,-1000);
            }
            _ptension2=_tension;





            double dy=_cb->_coords->_coords[1]-_py;
            if(_step%1000==0){
                auto fileN=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape_Scaled\Shape_%1_%2_FvT.txt)").arg(*_shape,*_title);
                if(_step==1000){

/*
                    if (QFile::exists(fileN))
                    {
                        QFile::remove(fileN);
                    }
*/
                }
                auto file=new QFile(fileN);

                if (!file->open(QIODevice::Append| QIODevice::WriteOnly | QIODevice::Text)){
                    assert(0)   ;
                    return;
                }
                auto out=new QTextStream(file);
                out->operator<<(_step*_DT);
                out->operator<<("\t");
                out->operator<<(_appliedF);

                out->operator<<("\t");
                out->operator<<(_radialForce);
                out->operator<<("\t");
                out->operator<<(_tension);
                out->operator<<("\t");
                out->operator<<(_cb->_coords->_coords[1]);
                out->operator<<("\r\n");
                out->flush();
                file->close();
                //double a=fmax(0,fmin(1,22-_cb->_coords->_coords[1]));
                double a=1;
                //_appliedF+=(20*(_step-_pstep)*_DT -dy)*20000*a;
                //_APFA=0;
                _APFA=(12*(_step-_pstep)*_DT -dy)*20*a;
                _RFA=-(12*(_step-_pstep)*_DT -dy)*10*a;
                _py=_cb->_coords->_coords[1];
                _pstep=_step;
            }
            //_appliedF=404;
/*
            _appliedF=0;
            _radialForce=0;*/
        }

        /*
        if(_tension>1e-4){
            _sf->_k=std::fmin(500,_sf->_k/10);
        }
        if(_tension<1e-5){
            _sf->_k=std::fmax(5000,_sf->_k*10);
        }*/
        //_py=_cb->_coords->_coords[1];
        /*
         if(_step%3000==0){
             _appliedF+=(_rand->generateDouble()-0.5)*0.1*_appliedF;
         }*/
        //_tension=(strain*(14.607392476665238)+param*(2.611268641852894)+(464.58341269938956));
        std::cout<<"TEN:"<<_THRESHOLD*_radiusFactor<<"\t"<<cR<<"\t"<<_initalArea<<"...\tSS:  "<<strain<<"\t"<<param<<" ten\t"<<_tension<<"\t"<<_kappaFactor<<"\t"<<_kappa<<"\t"<<_radiusFactor<<std::endl;
        std::cout<<"rForce"<<_radialForce<<"\t"<<_appliedF<<std::endl;
        //std::cout<<_appliedF<<"\t"<<_sf->_k<<"\t"<<_initalArea<<"\t"<<_frad<<"\t"<<_border->at(0)->_coords->len()<<"\t"<< _kappaFactor<<"_"<<mRdis<<"\t"<<bb->_coords->_coords[0]<<"\t"<<_Rind<<"\t"<<_radialForce<<std::endl;
        //std::cout<<_appliedF<<"\t"<<_sf->_k<<"\t"<<_initalArea<<"\t"<<_frad<<"\t"<<_border->at(0)->_coords->len()<<"\t"<< _kappaFactor<<"_"<<mRdis<<"\t"<<_Rind<<"\t"<<_radialForce<<std::endl;
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
void Physics::MembraneFromObj::capture(double *data,int len){
    SurfaceWithPhysics::capture(data,len);
    QString numb=QString::number(_step/1000);
    auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape_Scaled\Shape_%1_profile_%2_%3.txt)").arg(*_shape,*_title,numb);

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
    for(int i=0;i<_disc->_tris->size();i++){
        auto tri=_disc->_tris->at(i);
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

