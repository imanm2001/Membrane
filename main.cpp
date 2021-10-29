#include "mainwindow.h"
#include <QApplication>
#include <QCoreApplication>
#include "Phys/bead.h"
#include "Phys/vecd3d.h"
#include "iostream"

#include "utils/classes.h"
#include "Phys/bead.h"
#include "q3d/beadinfo.h"
#include "Phys/bendingparameters.h"
#include "q3d/triangle.h"
#include "Phys/vecd3d.h"
#include "Phys/tensor2.h"

int main(int argc, char *argv[])
{

    QApplication app(argc, argv);
    double a=10;
    double *b=&a;
    std::cout<<*b<<std::endl;

    MainWindow w;
    w.show();

    return app.exec();

}

//#include <QCoreApplication>
//#include "Phys/bead.h"
//#include "Phys/vecd3d.h"
//#include "Phys/springforce2.h"
//#include "iostream"
//#include <random>
//#include <fstream>

//#include <iostream>
//#include <iomanip>
//#include <string>
//#include <map>
//#include <random>
//#include <cmath>
//#include <math.h>

////using namespace std;

//int main(int argc, char *argv[])
//{
//    QCoreApplication a(argc, argv);
//    double D=10,gamma=KbT/D;

//    Physics::Bead *b1=new Physics::Bead(nullptr,new Physics::VecD3d(),gamma,0);
//    Physics::Bead *b2=new Physics::Bead(nullptr,new Physics::VecD3d(10,0,0),gamma,1);//D=?, ID=2
//    Physics::SpringForce2 *spring=new Physics::SpringForce2(10,10);
//    Physics::VecD3d *res=new Physics::VecD3d(); //res for b1-b2


//    // Write your code here





//    std::ofstream springfile1;
//    springfile1.open ("spring_dist.txt");


//    int i=1; //iterations
//    double p=1; //number of particles
//    double dt = 1e-3, t=0; //dt=1e-2 for msd


//    std::random_device rd{};
//    std::mt19937 gen{rd()};
//    std::normal_distribution<> d{0,1};
//    double noiseConstant=KbT*sqrt(2/(D*dt));
//    while (p<=100){ //for week 1 msd
//        t=0;i=0;
//        while (i<=20000){ //100 time steps for msd





//            double b1_randnum_x = d(gen);
//            double b1_randnum_y = d(gen);
//            double b1_randnum_z = d(gen);

//            double b1_fx=noiseConstant*b1_randnum_x;
//            double b1_fy=noiseConstant*b1_randnum_y*0;
//            double b1_fz=noiseConstant*b1_randnum_z*0;


//            b1->_force.setValues(b1_fx,b1_fy,b1_fz);




//            double b2_randnum_x = d(gen);
//            double b2_randnum_y = d(gen);
//            double b2_randnum_z = d(gen);

//            double b2_fx=noiseConstant*b2_randnum_x;
//            double b2_fy=noiseConstant*b2_randnum_y*0;
//            double b2_fz=noiseConstant*b2_randnum_z*0;

//            b2->_force.setValues(b2_fx,b2_fy,b2_fz);

//            spring->eval(b1,b2);

//            b1->update(dt);
//            b2->update(dt);


//            b1->_coords.subVec(&b2->_coords,res);//'&' turns reference into a pointer
//            double distance_b1b2=res->len();
//            if(i%200==0){
//                springfile1 << t << '\t' << distance_b1b2 << std::endl;

//            }

//            i++;
//            t+=dt;
//        }
//        p++;
//        springfile1.flush();
//    }
//    springfile1.flush();

//    springfile1.close();

//    std::cout << "Membrane simulation" << std::endl;

//    return a.exec();
//}
