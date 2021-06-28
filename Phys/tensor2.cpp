#include "tensor2.h"

Physics::Tensor2::Tensor2()
{
    _vecs[0]=new Physics::VecD3d();
    _vecs[1]=new Physics::VecD3d();
    _vecs[2]=new Physics::VecD3d();
}

void Physics::Tensor2::zero(){
    _vecs[0]->zero();
    _vecs[1]->zero();
    _vecs[2]->zero();

}
void Physics::Tensor2::produc(Physics::VecD3d *v1,Physics::VecD3d *vk){
    for(int c=0;c<3;c++){
        for(int k=0;k<3;k++){
            _vecs[c]->_coords[k]=v1->_coords[c]*vk->_coords[k];
        }
    }
}

void Physics::Tensor2::add(Tensor2* t){
    _vecs[0]->add(t->_vecs[0]);
    _vecs[1]->add(t->_vecs[1]);
    _vecs[2]->add(t->_vecs[2]);
}
void Physics::Tensor2::multConst(double c){
    _vecs[0]->multConst(c);
    _vecs[1]->multConst(c);
    _vecs[2]->multConst(c);
}
void Physics::Tensor2::addConst(double c){
    _vecs[0]->add(c,c,c);
    _vecs[1]->add(c,c,c);
    _vecs[2]->add(c,c,c);
}

void Physics::Tensor2::dotK(VecD3d* v,VecD3d* ret){
    ret->zero();
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ret->_coords[i]+=_vecs[j]->_coords[i]*v->_coords[j];
        }
    }
}
void Physics::Tensor2::dotC(VecD3d* v,VecD3d* ret){
    ret->zero();
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ret->_coords[i]+=_vecs[i]->_coords[j]*v->_coords[j];
        }
    }
}
void Physics::Tensor2::print(){
    _vecs[0]->print();
    _vecs[1]->print();
    _vecs[2]->print();
}
void Physics::Tensor2::I(){
    _vecs[0]->setValues(1,0,0);
    _vecs[1]->setValues(0,1,0);
    _vecs[2]->setValues(0,0,1);
}
void Physics::Tensor2::debug(){
    _vecs[0]->debug();
    _vecs[1]->debug();
    _vecs[2]->debug();
}
