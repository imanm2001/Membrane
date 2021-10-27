#include "quad.h"

Physics::Quad::Quad(Physics::Bead* mb,double k)
{
    _mb=mb;
    _sf=new Physics::SpringForce(k);
    _beads=(Physics::Bead**)malloc(sizeof(Physics::Bead*)*4);
    _beads[0]=nullptr;
    _beads[1]=nullptr;
    _beads[2]=nullptr;
    _beads[3]=nullptr;
    _temp=new Physics::VecD3d();

}

void Physics::Quad::addNewCandidate(Physics::Bead* b){
    _temp->setValues(_mb->_coords);
    _temp->sub(b->_coords);
    double len=_temp->len();
    for(int i=0;i<4;i++){

        if(_beads[0]!=nullptr){
            if(len<_lens[i]){
                for(int j=3;j>i;j--){
                    _beads[j]=_beads[j-1];
                    _lens[j]=_lens[j-1];
                }
                _beads[i]=b;
                _lens[i]=len;
                break;

            }

        }else{
            assert(i==0);
            _beads[0]=_beads[1]=_beads[2]=_beads[3]=b;
          _lens[0]=_lens[1]=_lens[2]=_lens[3]=len;
        }

    }
}

void Physics::Quad::eval(){
    for(int i=0;i<4;i++){
        _sf->eval(_mb,_beads[i],_lens[i]);

        _beads[i]->_force->zero();
    }
}
