#include "mthread.h"
#include <iostream>
MThread::MThread(updateFunc u):sleep(0),update(u)
{
_threadIsActive=false;
}
bool MThread::isActive(){
    return _threadIsActive;
}
void MThread::run()
{
    _threadIsActive=true;
   forever{
        if(!update()){
            break;
        }
        if(sleep>0){
            msleep(sleep);
        }

    }
_threadIsActive=false;
}

