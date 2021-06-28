#ifndef MTHREAD_H
#define MTHREAD_H

#include <QThread>
#include "utils/classes.h"
typedef void tupdate (int);

class MThread : public QThread
{
    Q_OBJECT
private:
updateFunc  update;
bool _threadIsActive;
public:
    MThread(updateFunc u);
    int sleep;
    bool isActive();
protected:
    void run();



};

#endif // MTHREAD_H
