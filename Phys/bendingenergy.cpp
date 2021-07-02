#include "bendingenergy.h"
#define SIGN(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)
Physics::BendingEnergy::BendingEnergy(BeadInfo *bi):_curv(1)
{
    _border=false;
    _bi=bi;
    _bi->setAttribute(BENDINGENERGY,this);
    _temp1=new VecD3d();
    _temp2=new VecD3d();
    _temp3=new VecD3d();
    _temps[0]=new VecD3d();
    _temps[1]=new VecD3d();
    _temps[2]=new VecD3d();
    _temps[3]=new VecD3d();
    _tempT[0]=new Tensor2();
    _tempT[1]=new Tensor2();
    _tempT[2]=new Tensor2();
    _tempT[3]=new Tensor2();
    _A=new VecD3d();
    _AvP=new VecD3d();

    _v1=new VecD3d();
    _v2=new VecD3d();
    _LB=new VecD3d();

    _n=new VecD3d();

}
double Physics::BendingEnergy::getChi(int j1,int j2){


    auto bij=_bi->_bendingParameters->at(j1);
    auto bijpm=_bi->_bendingParameters->at(j2);

    _temp2->setValues(bijpm->_dxP);
    _temp2->sub(bij->_dxP);
    _temp2->nomilize();
    double tt=bijpm->_dxP->dot(_temp2)/(bijpm->_l);
    if(tt==1||tt!=tt){
        std::cout<<"----"<<_border<<std::endl;
        _temp2->print();
        bijpm->_dxP->print();
        std::cout<<"----"<<std::endl;
        std::cout.flush();
        assert(0);

    }
    assert(tt==tt);
    assert(!isinf(tt));

    return tt;
}

void Physics::BendingEnergy::UP(int l,int j,int j2,Physics::VecD3d* t1,
                                Physics::VecD3d* t2,Physics::VecD3d* t3,Physics::Tensor2*tensor){


    auto bp1=_bi->_bendingParameters->at(j);
    auto bp2=_bi->_bendingParameters->at(j2);
    double dnum=1.0/(bp1->_lsq*bp2->_lsq);
    tensor->zero();
    if(l==-1){
        t1->setValues(bp2->_dxP);
        t1->cross(bp1->_dxP,t2);

        for(int k=0;k<3;k++){
            for(int c=0;c<3;c++){
                double d=0;
                for(int b=0;b<3;b++){
                    d+=LeviC[k][b][c]*(-bp1->_dxP->_coords[b]+bp2->_dxP->_coords[b]);
                }
                d*=bp1->_l*bp2->_l;
                double d2=bp1->_dxP->_coords[k]*bp2->_l/bp1->_l;
                d2+=bp2->_dxP->_coords[k]*bp1->_l/bp2->_l;
                d2*=t2->_coords[c];
                tensor->_vecs[c]->_coords[k]=((d+d2)*dnum);
            }
        }
    }else if(l==j){
        t1->setValues(bp1->_dxP);
        t1->cross(bp2->_dxP,t2);

        for(int k=0;k<3;k++){
            for(int c=0;c<3;c++){
                double d=0;
                for(int b=0;b<3;b++){
                    d-=LeviC[k][b][c]*bp2->_dxP->_coords[b];
                }
                d*=bp1->_l*bp2->_l;
                double d2=bp1->_dxP->_coords[k]*bp2->_l/bp1->_l;
                d2*=t2->_coords[c];
                tensor->_vecs[c]->_coords[k]=((d+d2)*dnum);
            }
        }
    }else if(l==j+1){
        for(int k=0;k<3;k++){
            for(int c=0;c<3;c++){
                double d=0;
                for(int b=0;b<3;b++){
                    d=LeviC[k][b][c]*bp1->_dxP->_coords[b];
                }
                d*=bp1->_l*bp2->_l;
                double d2=bp2->_dxP->_coords[k]*bp1->_l/bp2->_l;
                d2*=t2->_coords[c];
                tensor->_vecs[c]->_coords[k]=((d+d2)*dnum);
            }
        }
    }

}
void Physics::BendingEnergy::updateBendingParameters(){

    int size=_bi->_connections->size();
    assert(size>2);
    //    repetitiveConnections();
    _v1->zero();
    _v2->zero();
    _LB->zero();
    _Av=0;

    double Av2=0;
    for(int j=0;j<size;j++){

        int jp1=fixIndex(j+1,_bi->_connections),jp2=fixIndex(j-1,_bi->_connections);

        auto e=_bi->_connections->at(j);
        auto param=_bi->_bendingParameters->at(j);
        param->update(getTheOtherBead(e,_bi));
    }

    for(int j=0;j<size;j++){

        int jp1=fixIndex(j+1,_bi->_connections),jp2=fixIndex(j-1,_bi->_connections);

        auto e=_bi->_connections->at(j);
        auto param=_bi->_bendingParameters->at(j);
        auto paramp=_bi->_bendingParameters->at(jp1);
        _temp1->setValues(param->_dxP);
        _temp1->cross(paramp->_dxP,_temp2);
        param->_triA=_temp2->len()/2.0;

        double chi1=getChi(j,jp1), chi2=getChi(j,jp2);
        if(chi1*chi1>=1){
            std::cout<<chi1<<"\t"<<chi2<<"\t"<<_border<<std::endl;
        }
        assert(chi1*chi1<1);
        assert(chi2*chi2<1);
        /*
        double cd1=std::sqrt(1-chi1*chi1);
        double cd2=std::sqrt(1-chi2*chi2);
*/
        param->chiP=chi1;
        param->chiM=chi2;
        param->_cot1=chi1/std::sqrt(1-chi1*chi1);
        double cos2=getChi(jp1,j);
        param->_cot2=cos2/std::sqrt(1-cos2*cos2);
        param->_T=param->_cot1+chi2/std::sqrt(1-chi2*chi2);


        _temp1->setValues(param->_dxP);
        _temp1->multConst(param->_T);

        _temp1->debug();
        _v1->add(_temp1);

        auto tri=param->_tri;

        if(tri->_isObtuse){
            if(tri->_angles[tri->getVertexIndex(_bi)]<PI2){
                _Av+=param->_triA/4.0;
            }else{
                _Av+=param->_triA/2.0;
            }
        }else{
            double t=(param->_lsq*param->_cot1+paramp->_lsq*param->_cot2)/8.0;
            if(t<0){
              //  std::cout<<chi1<<"__\t"<<cos2<<std::endl;
            }
            //assert(t>0);
            _Av+=t;
        }

//        _Av+=param->_lsq*param->_T;

    }
    //_Av=_Av/8.0;
    //Av2/=8.0;


    assert(_Av==_Av);
    assert(_Av!=0);
    _Avs=_Av*_Av;
    _LB->setValues(_v1);
    _v1->debug();
    _v1->multConst(1.0/_Av);
    _v1->debug();
    double dn=_v1->dot(_v1);
    if(dn!=dn){
        qDebug()<<_Av;
    }
    assert(dn==dn);
    _A->zero();
    double Kg=2*3.14159265359;
    for(int j=0;j<size;j++){
        int jp=(j+1)%size;
        auto bp1=_bi->_bendingParameters->at(j);
        auto bp2=_bi->_bendingParameters->at(jp);
        bp1->_phi=std::acos(bp1->_dxP->dot(bp2->_dxP)/(bp1->_l*bp2->_l));
        Kg-=bp1->_phi;
        assert(bp1->_phi!=0);
        assert(bp1->_normal!=nullptr);
        _temp1->setValues(bp1->_normal);
        _temp1->multConst(bp1->_phi);
        _A->add(_temp1);
        UP(-1,j,fixIndex(j+1,_bi->_connections),_temp1,_temp2,_temp3,bp1->_tensor);
    }
    Kg/=_Av;

    _absAs=_A->dot(_A);
    _absA=std::sqrt(_absAs);
    if(_A0==-1){
        _A0=_absA;
    }

    _n->setValues(_A);
    _n->nomilize();
    _n->debug();
    //_triA/=3.0;

    _LBdotN=_LB->dot(_n)/2;
    _SIGN=SIGN(_LBdotN);



    //_Av=_triA/3.0;
    _H=_SIGN*_LB->len()/(4*_Av);

    //_H= _LB->len()/(4*_Av);
    double delta=std::sqrt(_H*_H-Kg);
    //std::cout<<_H+delta<<"\t"<<_H-delta<<"\t"<<(_H*_H-Kg)/Kg<<std::endl;
    //_H=3*_LB->len()/(4*_triA);
    Avp(-1,_temp1,_temp2,_temp3,_AvP);

}



void Physics::BendingEnergy::orderConnections(){
    int size=_bi->_connections->size();
    if(_bi->_bendingParameters->size()<size){
        int s=_bi->_bendingParameters->size();
        while(s<size){
            _bi->_bendingParameters->append(new Physics::BendingParameters(_bi));
            s++;
        }
    }
    int csize=_bi->_connections->size();
    if(csize>2){
        int numnull=0;


        //bool zero=false;
        bool switched=false;
        for(int i=0;i<_bi->_connections->size();i++){
            auto ed=_bi->_connections->at(i);
            if(ed->_tris[1]==nullptr){
                if(i>0){
                    if(!switched){
                        auto t=_bi->_connections->at(0);
                        _bi->_connections->replace(0,ed);
                        _bi->_connections->replace(i,t);

                        auto bp=_bi->_bendingParameters->at(0);
                        auto bp2=_bi->_bendingParameters->at(i);
                        _bi->_bendingParameters->replace(0,bp2);
                        _bi->_bendingParameters->replace(i,bp);
                    }else{
                        auto t=_bi->_connections->at(csize-1);
                        _bi->_connections->replace(0,ed);
                        _bi->_connections->replace(csize-1,t);

                        auto bp=_bi->_bendingParameters->at(csize-1);
                        auto bp2=_bi->_bendingParameters->at(i);
                        _bi->_bendingParameters->replace(csize-1,bp2);
                        _bi->_bendingParameters->replace(i,bp);
                    }
                }
                numnull++;
            }
        }
        if(numnull!=0&&numnull!=2){
            assert(0);
        }

        Geometry::Edge* firstE=_bi->_connections->at(0);

        Geometry::Triangle *tri=firstE->_tris[0];
        Geometry::Edge* e2=tri->edgeWithVertexExclude(_bi,firstE);

        _bi->_bendingParameters->at(0)->_normal=tri->_norm;
        _bi->_bendingParameters->at(0)->_tri=tri;
        int index=1;

        //zero|=tri->_curvater==0;

        while(e2!=firstE){

            Geometry::Edge *temp=_bi->_connections->at(index);
            int find=_bi->_connections->indexOf(e2);

            _bi->_connections->replace(index,e2);
            _bi->_connections->replace(find,temp);

            tri=e2->_tris[0]==tri?e2->_tris[1]:e2->_tris[0];

            if(tri!=nullptr){
                e2=tri->edgeWithVertexExclude(_bi,e2);

                //zero|=tri->_curvater==0;
                 auto bp=_bi->_bendingParameters->at(index);
                bp->_normal=tri->_norm;
                bp->_tri=tri;

            }else{
                _bi->_bendingParameters->at(index)->_normal=_bi->_bendingParameters->at(index-1)->_normal;
                assert(index==(_bi->_connections->size()-1));
                index++;
                break;
            }
            index++;
        }
        assert(index==(_bi->_connections->size()));
        if(numnull==2){
            _curv=0;
            if(!_border){
                _bi->setAttribute(BORDER,(QObject*) new Geometry::QBoolean(true));
                _border=true;
            }
        }else{
            _curv=tri->_curvater;
            ;
            //_curv/=_bi->_connections->size();
        }
    }else{
        assert(0);
        _curv=0;
        if(!_border){
            _bi->setAttribute(BORDER,(QObject*) new Geometry::QBoolean(true));
            _border=true;
        }
    }
}

void Physics::BendingEnergy::LsqP(int l,int j,Physics::VecD3d *ret){
    ret->setValues(_bi->_bendingParameters->at(j)->_dxP);

    if(j==l){
        ret->multConst(-2);

    }else if(l==-1){
        ret->multConst(2);

    }else{
        ret->multConst(0);
    }

}

void Physics::BendingEnergy::TP(int l,int j,Physics::VecD3d *temp,Physics::VecD3d *temp2,Physics::VecD3d *ret){
    auto bp= _bi->_bendingParameters->at(j);

    double t1=(1-bp->chiM*bp->chiM);
    assert(t1>0);
    t1=std::sqrt(t1*t1*t1);
    chiP(l,j,fixIndex(j-1,_bi->_connections),temp,temp2);
    temp2->multConst(1.0/t1);
    temp2->debug();
    ret->setValues(temp2);

    t1=(1-bp->chiP*bp->chiP);
    t1=std::sqrt(t1*t1*t1);
    chiP(l,j,fixIndex(j+1,_bi->_connections),temp,temp2);
    temp2->multConst(1.0/t1);
    temp2->debug();
    ret->add(temp2);

}
void Physics::BendingEnergy::TP1(int l,int j,Physics::VecD3d *temp,Physics::VecD3d *ret){
    auto bp= _bi->_bendingParameters->at(j);

    double t1=(1-bp->chiP*bp->chiP);
    assert(t1>0);
    t1=std::sqrt(t1*t1*t1);
    chiP(l,j,fixIndex(j+1,_bi->_connections),temp,ret);
    ret->multConst(1.0/t1);
    ret->debug();


}

void Physics::BendingEnergy::TP2(int l,int j,Physics::VecD3d *temp,Physics::VecD3d *ret){
    int jp=fixIndex(j+1,_bi->_connections);
    auto bp= _bi->_bendingParameters->at(jp);

    double t1=(1-bp->chiM*bp->chiM);
    assert(t1>0);
    t1=std::sqrt(t1*t1*t1);
    chiP(l,jp,j,temp,ret);
    ret->multConst(1.0/t1);
    ret->debug();


}

void Physics::BendingEnergy::chiP(int l,int j,int m,Physics::VecD3d *temp,Physics::VecD3d *ret){

    auto bj=_bi->_bendingParameters->at(j);
    auto bm=_bi->_bendingParameters->at(m);

    double chi=getChi(j,m);

    temp->setValues(bm->_dxP);
    temp->sub(bj->_dxP);
    double lim=bm->_l;
    double ljm=temp->len();
    ret->zero();

    if(l==-1){
        ret->setValues(bm->_dxP);
        ret->multConst(-ljm*chi/lim);
        ret->add(temp);
    }else if(l==j){
        ret->setValues(temp);
        ret->multConst(-lim*chi/ljm);
        ret->add(bm->_dxP);
    }else if(l==m){
        ret->setValues(temp);
        ret->add(bm->_dxP);

        temp->multConst(-lim*chi/ljm);
        ret->add(temp);

        temp->setValues(bm->_dxP);
        temp->multConst(-ljm*chi/lim);

        ret->add(temp);
        ret->multConst(-1);

    }
    ret->multConst(1.0/(ljm*lim));

    assert(ljm!=0);
    assert(lim!=0);
    assert(lim==lim);
    assert(ljm==ljm);
    ret->debug();
}
void Physics::BendingEnergy::Ep(int l,VecD3d * ret){
    ret->zero();

    if(!_border){
        double HH=_H-_curv;
        HH=_H;
        if(l==-1){
            Hp(-1,_temps,_tempT,_temp1);
            _temp1->multConst(2*_Av*(HH));
            _temp1->debug();
            _temp2->setValues(_AvP);

            _temp2->multConst(HH*HH);
            _temp2->debug();
            ret->add(_temp1);
            ret->add(_temp2);
            ret->debug();
            for(int j=0;j<_bi->_connections->size();j++){
                auto b=getTheOtherBead(_bi->_connections->at(j),_bi);
                int jp=b->findBeadIndexInTheConnection(this->_bi);
                BendingEnergy* be=(BendingEnergy*)b->getAttribute(BENDINGENERGY);
                be->Ep(jp,_temp3);
                _temp3->debug();

                ret->add(_temp3);
                ret->debug();
            }
            ret->debug();



        }else{
            Hp(l,_temps,_tempT,_temp1);
            _temp1->debug();
            _temp1->multConst(2*_Av*HH);

            _temp1->debug();
            Avp(l,_temps[0],_temps[1],_temps[2],_temp2);
            _temp2->multConst(HH*HH);
            _temp2->debug();
            ret->add(_temp1);
            ret->add(_temp2);
        }
    }
    ret->debug();
}

void Physics::BendingEnergy::absLp(int l,int j,VecD3d* ret){

    if(l==-1||l==j){
        ret->setValues(this->_bi->_bendingParameters->at(j)->_dxP);
        if(l==j){
            ret->multConst(-1);
        }
    }else{
        ret->zero();
    }

}
void Physics::BendingEnergy::Ep2(int l,VecD3d * ret){
    ret->zero();




    if(l==-1){


        H2p(-1,_temps,_tempT,_temp1);


        ret->setValues(_temp1);
        //  HpSV(-1,_temps,_temp1);





        ret->debug();
        for(int j=0;j<_bi->_connections->size();j++){
            auto b=getTheOtherBead(_bi->_connections->at(j),_bi);
            int jp=b->findBeadIndexInTheConnection(this->_bi);
            BendingEnergy* be=(BendingEnergy*)b->getAttribute(BENDINGENERGY);

            be->Ep2(jp,_temp3);

            _temp3->debug();

            ret->add(_temp3);
            ret->debug();
        }
        ret->debug();
        ret->multConst(2.0);


    }else{
        H2p(l,_temps,_tempT,ret);
        _temp1->debug();


        //HpSV(l,_temps,_temp1);

    }

    ret->debug();
}


void Physics::BendingEnergy::Atp(int l,Physics::VecD3d* t1,Physics::VecD3d* t2,Physics::VecD3d*t3,Physics::VecD3d* ret){
    ret->zero();
    int st=0,end=_bi->_connections->size();

    if(l!=-1){
        st=l-1;
        end=l+2;
    }

    for(int jp=st;jp<end;jp++){
        int j=fixIndex(jp,_bi->_connections);
        auto A=_bi->_bendingParameters->at(j);
        int j1=fixIndex(jp+1,_bi->_connections);
        double dA=deltaDelta(l,j,l,-1);
        double dB=deltaDelta(l,j1,l,-1);
        auto B=_bi->_bendingParameters->at(j1);
        t1->setValues(A->_dxP);
        t1->cross(B->_dxP,t2);
        double triA=t2->len()/2.0;
        double A2=A->_dxP->dot(A->_dxP);
        double B2=B->_dxP->dot(B->_dxP);
        double AB=A->_dxP->dot(B->_dxP);
        for(int i=0;i<3;i++){
            t3->_coords[i]=A->_dxP->_coords[i]*B2*dA+B->_dxP->_coords[i]*A2*dB-AB*(B->_dxP->_coords[i]*dA+A->_dxP->_coords[i]*dB);
        }
        t3->multConst(-1.0/(4*triA));
        ret->add(t3);

    }

    ret->debug();
}

//old
void Physics::BendingEnergy::Ap(int l,Physics::VecD3d* t1,Physics::VecD3d* t2,Physics::VecD3d*t3,Physics::VecD3d* ret){
    ret->zero();
    int st=0,end=_bi->_connections->size();

    if(l!=-1){
        st=l-1;
        end=l+2;
    }

    for(int jp=st;jp<end;jp++){
        int j=fixIndex(jp,_bi->_connections);
        int j1=fixIndex(jp+1,_bi->_connections);
        auto bp=_bi->_bendingParameters->at(j);

        auto tri=bp->_tri;

        if(tri->_isObtuse){
            double cte=2;
            if(tri->_angles[tri->getVertexIndex(_bi)]<PI2){
                cte=4;
            }
            double dA=deltaDelta(l,-1,l,j);
            double dB=deltaDelta(l,-1,l,j1);
            auto B=_bi->_bendingParameters->at(j1);
            t1->setValues(bp->_dxP);
            t1->cross(B->_dxP,t2);
            double triA=t2->len()/2.0;
            double A2=bp->_dxP->dot(bp->_dxP);
            double B2=B->_dxP->dot(B->_dxP);
            double AB=bp->_dxP->dot(B->_dxP);
            for(int i=0;i<3;i++){
                t3->_coords[i]=bp->_dxP->_coords[i]*B2*dA+B->_dxP->_coords[i]*A2*dB-AB*(B->_dxP->_coords[i]*dA+bp->_dxP->_coords[i]*dB);
            }
            t3->multConst(1.0/(cte*4*triA));
            ret->add(t3);
        }else{
            TP1(l,j,t1,t3);
            t3->multConst(bp->_lsq);
            t3->debug();
            auto bp2=_bi->_bendingParameters->at(j1);
            TP2(l,j,t1,t2);
            t2->multConst(bp2->_lsq);
            t2->debug();
            t3->add(t2);

            LsqP(l,j,t1);
            t1->multConst(bp->_cot1);
            t1->debug();
            t3->add(t1);

            LsqP(l,j1,t1);
            t1->multConst(bp->_cot2);
            t1->debug();
            t3->add(t1);
            t3->multConst(1/8.0);
            ret->add(t3);

        }


    }

    ret->debug();
}
void Physics::BendingEnergy::Avp(int l,Physics::VecD3d* t1,Physics::VecD3d* t2,Physics::VecD3d*t3,Physics::VecD3d* ret){
    ret->zero();
    int st=0,end=_bi->_connections->size();

    if(l!=-1){
        st=l-1;
        end=l+2;
    }

    for(int jp=st;jp<end;jp++){
        int j=fixIndex(jp,_bi->_connections);
        auto bp=_bi->_bendingParameters->at(j);

        TP(l,j,t1,t2,t3);

        t3->multConst(bp->_lsq);
        t3->debug();
        LsqP(l,j,t1);
        t1->multConst(bp->_T);
        t1->debug();
        ret->add(t1);
        ret->add(t3);

    }
    ret->multConst(1/8.0);
    ret->debug();
}

void Physics::BendingEnergy::Hp(int l, VecD3d **tmps,Tensor2 **tn, VecD3d *ret){
    if(!_border){
        tn[0]->zero();
        tn[1]->zero();
        tn[2]->zero();
        tn[3]->zero();
        ret->zero();
        int st=0,end=_bi->_connections->size();
        if(l!=-1){
            st=l-1;
            end=l+2;
        }
        nP(l,tmps[0],tmps[1],tmps[2],tn[0],tn[1],tn[2]);
        _temp1->zero();

        for(int jp=st;jp<end;jp++){
            int j=fixIndex(jp,_bi->_connections);
            auto bp=_bi->_bendingParameters->at(j);
            //First Term
            TP(l,j,tmps[0],tmps[1],tmps[2]);
            tn[0]->produc(bp->_dxP,tmps[2]);
            double t=bp->_T*deltaDelta(l,-1,j,l);
            for(int n=0;n<3;n++){
                tn[0]->_vecs[n]->_coords[n]+=t;
            }

            tn[3]->add(tn[0]);
            _temp2->setValues(bp->_dxP);
            _temp2->multConst(bp->_T/(2*_Av));
            _temp1->add(_temp2);
            //Second Term
            tn[2]->dotC(bp->_dxP,tmps[0]);
            tmps[0]->debug();
            tmps[0]->multConst(bp->_T);
            ret->add(tmps[0]);
        }

        ret->debug();
        tn[3]->dotC(_n,tmps[0]);
        ret->add(tmps[0]);
        ret->multConst(_Av);
        if(l==-1){
            tmps[3]->setValues(_AvP);
        }else{
            Avp(l,tmps[0],tmps[1],tmps[2],tmps[3]);
        }
        tmps[3]->multConst(-_LBdotN);
        tmps[3]->debug();
        ret->add(tmps[3]);
        ret->multConst(1.0/_Avs);
        ret->debug();
    }

}
void Physics::BendingEnergy::HpSV(int l, VecD3d **tmps,VecD3d *ret){
    tmps[0]->zero();
    tmps[1]->zero();
    tmps[2]->zero();
    tmps[3]->zero();
    ret->zero();
    double x,y,z;
    x=y=z=0;
    for(int jp=0;jp<_bi->_connections->size();jp++){
        int j=fixIndex(jp,_bi->_connections);
        auto bp=_bi->_bendingParameters->at(j);
        //First Term
        TP(l,j,tmps[0],tmps[1],tmps[2]);
        //*** tmps[2]=Tp


        //*** tn[0] = (x(i)-x(j) x Tp)
        double t=bp->_T*deltaDelta(l,-1,j,l);
        //tn[0]->produc(bp->_dxP,tmps[2]);
        double a=_LB->dot(bp->_dxP);
        tmps[0]->setValues(tmps[2]);
        tmps[0]->multConst(a);

        x+=tmps[0]->_coords[0]+t*_LB->_coords[0];
        y+=tmps[0]->_coords[1]+t*_LB->_coords[1];
        z+=tmps[0]->_coords[2]+t*_LB->_coords[2];
    }


    ret->setValues(x,y,z);
    ret->multConst(1/(_LB->len()));
    //ret->add(tmps[3]);

}

//void Physics::BendingEnergy::H2p(int l, VecD3d **tmps,Tensor2 **tn, VecD3d *ret){
//    tn[0]->zero();
//    tn[1]->zero();
//    tn[2]->zero();
//    tn[3]->zero();
//    ret->zero();


//    int st=0,end=_bi->_connections->size();
////    if(l!=-1){
////        st=l-2;
////        end=l+2;
////    }

//    //Second term -2[\{sum [x(i)-x(j)Tij]^2/Av}*sum[Tijlsqp+lsq*Tp]
//    tmps[3]->zero();

//    double x,y,z;
//    x=y=z=0;
//    for(int jp=st;jp<end;jp++){
//        int j=fixIndex(jp,_bi->_connections);
//        auto bp=_bi->_bendingParameters->at(j);
//        //First Term
//        TP(l,j,tmps[0],tmps[1],tmps[2]);
//        //*** tmps[2]=Tp


//        //*** tn[0] = (x(i)-x(j) x Tp)
//        double t=bp->_T*deltaDelta(l,-1,j,l);
//        //tn[0]->produc(bp->_dxP,tmps[2]);
//        double a=_LB->dot(bp->_dxP);
//        tmps[0]->setValues(tmps[2]);
//        tmps[0]->multConst(a);

//        x+=tmps[0]->_coords[0]+t*_LB->_coords[0];
//        y+=tmps[0]->_coords[1]+t*_LB->_coords[1];
//        z+=tmps[0]->_coords[2]+t*_LB->_coords[2];

//        //***(x(i)-x(j) x Tp + t.ek = (x(i)-x(j)) x Tp +tI
//        /*
//            for(int n=0;n<3;n++){
//                tn[0]->_vecs[n]->_coords[n]+=t;
//            }*/


//        //            tn[3]->add(tn[0]);


//        //Second Term
//        //*** Tij x lsqp
//        LsqP(l,j,tmps[0]);
//        tmps[0]->multConst(bp->_T);
//        tmps[3]->add(tmps[0]);
//        //*** lsq x Tp
//        tmps[0]->setValues(tmps[2]);
//        tmps[0]->multConst(bp->_lsq);
//        tmps[3]->add(tmps[0]);



//    }
//    double C=-_LB->dot(_LB)/(32*_Av*_Av);


//    tmps[3]->multConst(C);
//    ret->setValues(tmps[3]);

//    ret->debug();

//    //tn[3]->dotK(_LB,tmps[0]);
//    tmps[0]->setValues(x,y,z);
//    tmps[0]->multConst(1/(2*_Av));
//    ret->add(tmps[0]);
//    ret->debug();



//}



//absLBp
void Physics::BendingEnergy::LBp(int l,VecD3d ** tmps,VecD3d *ret){
    tmps[0]->zero();
    tmps[1]->zero();
    tmps[2]->zero();
    tmps[3]->zero();
    ret->zero();
    double x,y,z;
    x=y=z=0;
    int st=0,end=_bi->_connections->size();

    if(l!=-1){
        st=l-1;
        end=l+2;
    }

    for(int jp=st;jp<end;jp++){
        int j=fixIndex(jp,_bi->_connections);
        auto bp=_bi->_bendingParameters->at(j);
        TP(l,j,tmps[0],tmps[1],tmps[2]);
        double t=bp->_T*deltaDelta(l,-1,j,l);
        //tn[0]->produc(bp->_dxP,tmps[2]);
        double a=_LB->dot(bp->_dxP);
        tmps[0]->setValues(tmps[2]);
        tmps[0]->multConst(a);

        x+=tmps[0]->_coords[0]+t*_LB->_coords[0];
        y+=tmps[0]->_coords[1]+t*_LB->_coords[1];
        z+=tmps[0]->_coords[2]+t*_LB->_coords[2];
    }
    ret->setValues(x,y,z);
    ret->multConst(1/(_LB->len()));
}
void Physics::BendingEnergy::H2p(int l, VecD3d **tmps,Tensor2 **tn, VecD3d *ret){

    double HH0=_H-_curv;

    tn[0]->zero();
    tn[1]->zero();
    tn[2]->zero();
    tn[3]->zero();
    ret->zero();




    //Second term -2[\{sum [x(i)-x(j)Tij]^2/Av}*sum[Tij lsqp+lsq*Tp]
    tmps[3]->zero();

    //tmps[2] <-LBp
    //tmps[3] <-Avp
    LBp(l,tmps,ret);
    //ret->multConst(SIGN(_LBdotN));
    //std::cout<<"----"<<std::endl;
    ///Avp(l,tmps[0],tmps[1],tmps[2],tmps[3]);
    //tmps[3]->print();
    Ap(l,tmps[0],tmps[1],tmps[2],tmps[3]);
    //tmps[3]->print();


    tmps[2]->setValues(ret);




    //tmp[0] <- LBp
    tmps[0]->setValues(tmps[2]);

    //tmp[0]<- LBp*_AV (moshtaqe soorat dar makhraj
    tmps[0]->multConst(_Av);



    //tmps[1] <- Avp
    tmps[1]->setValues(tmps[3]);
    //tmps[1] <- -LB*Avp (-moshtaqe makhraj dar soorat)
    tmps[1]->multConst(-_LB->len());

    tmps[0]->add(tmps[1]);
    //tmps[0] <- moshtaqe sorat*makhraj-moshtaqe makhraj dar soorat
    tmps[0]->multConst( HH0*_SIGN/(2*_Av));
    //tmps[0] <- moshtaqe (sorat*makhraj-moshtaqe makhraj dar soorat)*(H-H0)*Av/(makhraj^2)
    ret->setValues(tmps[0]);


    //Derivative of surface

    tmps[3]->multConst(HH0*HH0);
    ret->add(tmps[3]);



}

/* OLD2
void Physics::BendingEnergy::H2p(int l, VecD3d **tmps,Tensor2 **tn, VecD3d *ret){

    double HH0=_H-_curv;
    tn[0]->zero();
    tn[1]->zero();
    tn[2]->zero();
    tn[3]->zero();
    ret->zero();


    int st=0,end=_bi->_connections->size();
    //    if(l!=-1){
    //        st=l-2;
    //        end=l+2;
    //    }

    //Second term -2[\{sum [x(i)-x(j)Tij]^2/Av}*sum[Tij lsqp+lsq*Tp]
    tmps[3]->zero();
    LBp(l,tmps,ret);
    Avp(l,tmps[0],tmps[1],tmps[2],tmps[3]);
    tmps[2]->setValues(ret);
    //tmps[2] <-LBp
    //tmps[3] <-Avp
    tmps[0]->setValues(tmps[2]);
    tmps[0]->multConst(_Av);



    tmps[1]->setValues(tmps[3]);
    tmps[1]->multConst(-_LB->len());
    tmps[0]->add(tmps[1]);
    tmps[0]->multConst(HH0/(2*_Av));
    ret->setValues(tmps[0]);
    //Derivative of surface
//    tmps[3]->multConst(HH0*HH0);
  //  ret->add(tmps[3]);


}
*/
void Physics::BendingEnergy::nP(int l,Physics::VecD3d* t1,Physics::VecD3d* t2,Physics::VecD3d* t3,Physics::Tensor2 *tt,Physics::Tensor2 *tt2,Physics::Tensor2 *ret){
    ret->zero();

    t1->zero();
    t2->zero();
    t3->zero();
    int size=_bi->_connections->size();
    int sp=0,ep=size;
    if(l!=-1){
        sp=l-1;
        ep=l+2;
    }
    tt->zero();
    for(int jp=sp;jp<ep;jp++){
        int j=fixIndex(jp,_bi->_connections);
        auto bp=_bi->_bendingParameters->at(j);
        int j2=(j+1)%size;
        PhiP(l,j,j2,t1,t2,t3);
        t3->debug();
        tt2->produc(bp->_normal,t3);
        tt2->debug();
        tt->add(tt2);
        tt->debug();
        UP(l,j,j2,t1,t2,t3,tt2);
        tt->add(tt2);
        tt->debug();
    }
    tt->dotC(_A,t1);
    t1->debug();
    /*
    tt2->produc(_A,t1);
    tt2->multConst(1.0/(_absA*_absAs));*/
    ret->produc(_A,t1);
    ret->multConst(-1.0/(_absA*_absAs));

    tt->multConst(1.0/_absA);
    ret->add(tt);
    //    ret->zero();
}
void Physics::BendingEnergy::PhiP(int l,int j,int j2,Physics::VecD3d* t1,Physics::VecD3d* t2,Physics::VecD3d* ret){
    auto bp1=_bi->_bendingParameters->at(j);
    auto bp2=_bi->_bendingParameters->at(j2);
    double dot=bp1->_dxP->dot(bp2->_dxP);
    double t=dot/(bp1->_l*bp2->_l);

    t*=t;
    assert(t==t);
    assert(t<1);
    t=std::sqrt(1-t);
    double denum=bp1->_lsq*bp2->_lsq*t;
    if(denum==0){
        bp1->_dxP->print();
        bp2->_dxP->print();
    }
    assert(denum==denum);
    assert(denum>0);
    ret->zero();
    if(l==-1){
        t1->setValues(bp1->_dxP);
        t1->add(bp2->_dxP);
        t1->multConst(-bp1->_l*bp2->_l);
        ret->setValues(t1);
        t1->setValues(bp1->_dxP);
        t1->multConst(dot*bp2->_l/bp1->_l);
        t2->setValues(bp1->_dxP);
        t2->multConst(dot*bp1->_l/bp2->_l);
        t1->add(t2);
        ret->add(t1);
        ret->multConst(1.0/denum);
        ret->debug();
    }else if(l==j){
        t1->add(bp2->_dxP);
        t1->multConst(bp1->_l*bp2->_l);
        ret->setValues(t1);
        t1->setValues(bp1->_dxP);
        t1->multConst(-dot*bp2->_l/bp1->_l);
        ret->add(t1);
        ret->multConst(1.0/denum);
        ret->debug();
    }else if(l==j+1){
        t1->add(bp1->_dxP);
        t1->multConst(bp1->_l*bp2->_l);
        ret->setValues(t1);
        t1->setValues(bp2->_dxP);
        t1->multConst(-dot*bp1->_l/bp2->_l);
        ret->add(t1);
        ret->multConst(1.0/denum);
        ret->debug();

    }

    ret->debug();
}

double Physics::BendingEnergy::calE(){
    double ret=0;
    ret=_H-_curv;
    ret*=ret;
    ret*=2*_Av;
    return ret;
}
