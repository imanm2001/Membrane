#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow),_mutex(),_init(false)
{
    std::cout<<"INIT"<<std::endl;
    _tri1Ind=-1;
    ui->setupUi(this);


    connect(ui->focusPushButton,&QPushButton::clicked,this,&MainWindow::focusCamera);
    connect(ui->listWidget,&QListWidget::itemSelectionChanged,this,&MainWindow::listChanged);
    connect(ui->pushButton,&QPushButton::clicked,this,&MainWindow::pushButtClicked);
    connect(ui->swapPushButton,&QPushButton::clicked,this,&MainWindow::swapPushButtClicked);
    connect(ui->randomPushButton,&QPushButton::clicked,this,&MainWindow::randomPushButtClicked);
    connect(ui->captureButton,&QPushButton::clicked,this,&MainWindow::captureButtClicked);
    _running=true;
    _mt=_mt2=nullptr;
}
void MainWindow::focusCamera(bool b){
    ui->widget3D->lookAtCenter();
}
void MainWindow::showEvent(QShowEvent *)
{
    if(!_init){

        ui->widget3D->init();
        Physics::SIM=(Physics::Debuggable*) ui->widget3D->_mesh;
        IERROR(0);


        /*for(int i=0;i<ui->widget3D->_plan->_tris.size();i++){
            ui->listWidget->addItem(QString::number(i));
        }
        ui->widget3D->_plan->updateIndecies();

        _memb=ui->widget3D->_memb;
        */
        std::cout<<"INIT\r\n";

        auto u2=std::bind(&MainWindow::updatePlan,this);
        _mt=new MThread(u2);
        _mt->sleep=1000;
        _mt->start();


        auto u=std::bind(&MainWindow::updateMembrane,this);
        _mt2=new MThread(u);
        _mt2->sleep=100;
        _mt2->start();
        _init=true;
        this->setWindowTitle(ui->widget3D->_mesh->getTitle());
    }

}

MainWindow::~MainWindow()
{

    std::cout<<"terminate"<<std::endl;
    _running=false;
    if(_mt!=nullptr){
        while(_mt->isActive()){
            QThread::currentThread()->msleep(100);
        }
        delete  _mt;
    }
    if(_mt2!=nullptr){
        while(_mt2->isActive()){
            QThread::currentThread()->msleep(100);
        }
        delete  _mt2;
    }
    //_mt->terminate();
    //_mt2->terminate();
    std::cout<<"AAAA"<<std::endl;
    delete ui;

}
void MainWindow::listChanged(){
    /*
    int i=ui->listWidget->currentRow();

    ui->widget3D->drawTriangle(i,_tri1Ind==-1?ui->widget3D->_triGem:ui->widget3D->_triGem2);
    if(_tri1Ind!=-1){
        int tri2Ind=ui->listWidget->currentRow();
        if(tri2Ind!=-1){
            Geometry::Triangle *tri1=ui->widget3D->_plan->_tris.at(_tri1Ind);
            Geometry::Triangle *tri2=ui->widget3D->_plan->_tris.at(tri2Ind);
            ui->swapPushButton->setEnabled(_tri1Ind!=tri2Ind&& tri1->commonEdge(tri2)!=nullptr);
        }
    }*/
}
void MainWindow::swapPushButtClicked(bool b){
    /*
    int tri2Ind=ui->listWidget->currentRow();
    if(tri2Ind!=-1){

        Geometry::Triangle *tri1=ui->widget3D->_plan->_tris.at(_tri1Ind);
        Geometry::Triangle *tri2=ui->widget3D->_plan->_tris.at(tri2Ind);
        ui->widget3D->_plan->swapEdge(tri1,tri2);
        ui->widget3D->drawTriangle(_tri1Ind,ui->widget3D->_triGem);
        ui->widget3D->drawTriangle(tri2Ind,ui->widget3D->_triGem2);
        ui->widget3D->_plan->updateIndecies();
        ui->widget3D->_plan->updateVertecies();


    }
*/


}
bool MainWindow::updatePlan(){
    //std::cout<<"1PlanL"<<std::endl;
    //QMutexLocker locker(&_mutex);

    if(QThread::currentThread()!=QApplication::instance()->thread()){

        _mutex.lock();
        //ui->widget3D->_plan->updateVIN();
        ui->widget3D->_mesh->updateVIN();
        _mutex.unlock();
        QMetaObject::invokeMethod(this,
                                  "updatePlan",
                                  Qt::QueuedConnection);     // val1F
    }else{
        ui->widget3D->update();
        if(ui->captureCheckBox->isChecked()){
            ui->widget3D->snapshot();
        }
        if(_mutex.tryLock()){
            ui->widget3D->_view->requestUpdate();

            _mutex.unlock();
        }
    }





    //std::cout<<"3PlanUL\r\n"<<std::endl;
    return _running;
}

bool MainWindow::updateMembrane(){

    if( ui->checkBox->isChecked()){
        _mt2->sleep=0;
        //      std::cout<<"1MemL"<<std::endl;
        //QMutexLocker locker(&_mutex);

        //        std::cout<<"2Mem"<<std::endl;
        if(_mutex.tryLock()){

            for(int i=0;i<200&&_running;i++){

                //_memb->update();
                ui->widget3D->_mesh->update();


            }
            /*double l=ui->widget3D->_balloon->_l;
            double lsq=ui->widget3D->_balloon->_lsq;
            double A=ui->widget3D->_balloon->_totalA;

            ui->lineEdit->setText(QString("%1 %2").arg(QString::number(l),QString::number(std::sqrt(lsq-l*l))));*/
            //ui->lineEdit->setText(QString("%1 %2").arg(QString::number(l),QString::number(A)));
            _mutex.unlock();
        }





        //   std::cout<<"3MemUL"<<std::endl;
    }else{
        _mt2->sleep=500;
    }


    /*
    if(QThread::currentThread()!=QApplication::instance()->thread()){
        QMetaObject::invokeMethod(this,         // obj
                                  "tupdate",         // member: don't put parameters
                                  Qt::QueuedConnection,     // connection type
                                  Q_ARG(int, i));     // val1
    }else{
        ui->progressBar->setValue(i);
    }

*/
    return _running;
}

void MainWindow::pushButtClicked(bool b){
    _tri1Ind=ui->listWidget->currentRow();
    ui->widget3D->drawTriangle(_tri1Ind,ui->widget3D->_triGem);
    ui->triNameLabel->setText(QString::number(_tri1Ind));
    ui->widget3D->_triGem2->parentNode()->setEnabled(true);
}

void MainWindow::randomPushButtClicked(bool b){

}
void MainWindow::captureButtClicked(bool){
    ui->widget3D->_mesh->saveToFile();
}
