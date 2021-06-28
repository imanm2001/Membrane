#ifndef TENSOR2_H
#define TENSOR2_H
#include "utils/classes.h"
#include "Phys/vecd3d.h"
namespace Physics {


class Tensor2
{
public:
    Physics::VecD3d *_vecs[3];
    Tensor2();
    void zero();
    void produc(Physics::VecD3d *v1,Physics::VecD3d *vk);
    void add(Tensor2*);
    void addConst(double);
    void multConst(double);
    void dotC(VecD3d*,VecD3d*);
    void dotK(VecD3d*,VecD3d*);
    void print();
    void I();
    void debug();
};
}

#endif // TENSOR2_H
