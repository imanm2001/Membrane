#ifndef BENDINGENERGY_H
#define BENDINGENERGY_H
#include <QObject>

#include "utils/classes.h"
#include "Phys/bead.h"
#include "q3d/beadinfo.h"
#include "Phys/bendingparameters.h"
#include "q3d/triangle.h"
#include "vecd3d.h"
namespace Physics {

using Geometry::BeadInfo;
class BendingEnergy:public QObject
{
private:


    VecD3d *_temps[4];
    Tensor2 *_tempT[4];
    double _Avs;
    double _SIGN;
    double _A0,_absA,_absAs,_LBdotN;
    Physics::VecD3d *_v1,*_v2,*_LB,*_temp1,*_temp2,*_temp3,*_AvP,*_A;

    void Ap(int,Physics::VecD3d*,Physics::VecD3d*,Physics::VecD3d*,Physics::VecD3d*);
    void Avp(int,Physics::VecD3d*,Physics::VecD3d*,Physics::VecD3d*,Physics::VecD3d*);
    void Atp(int,Physics::VecD3d*,Physics::VecD3d*,Physics::VecD3d*,Physics::VecD3d*);

    double getChi(int j1,int j2);
    void UP(int,int,int,Physics::VecD3d*,Physics::VecD3d*,Physics::VecD3d*,Physics::Tensor2* );
    void LsqP(int l,int j,Physics::VecD3d *ret);
    void TP(int,int,Physics::VecD3d *,Physics::VecD3d *,Physics::VecD3d *);
    void TP1(int,int,Physics::VecD3d *,Physics::VecD3d *);
    void TP2(int,int,Physics::VecD3d *,Physics::VecD3d *);
    void chiP(int,int,int,Physics::VecD3d *,Physics::VecD3d *);


    void Hp(int, VecD3d **,Tensor2 **,VecD3d *);
    void H2p(int, VecD3d **,Tensor2 **,VecD3d *);
    void HpSV(int, VecD3d **,VecD3d *);
    void nP(int l,Physics::VecD3d* t1,Physics::VecD3d* t2,Physics::VecD3d* t3,Physics::Tensor2 *tt,Physics::Tensor2 *tt2,Physics::Tensor2 *ret);
    void PhiP(int,int,int,Physics::VecD3d*,Physics::VecD3d*,Physics::VecD3d*);
    void absLp(int,int,VecD3d*);
    void LBp(int,VecD3d **,VecD3d *);

public:
    bool _border;
    double _curv,_Av,_H;
    BeadInfo *_bi;
    Physics::VecD3d *_n;
    BendingEnergy(BeadInfo *);
    void Ep(int,VecD3d*);
    void Ep2(int,VecD3d*);
    void updateBendingParameters();
    void orderConnections();
    double calE();
    double calcSig();

};
}

#endif // BENDINGENERGY_H
