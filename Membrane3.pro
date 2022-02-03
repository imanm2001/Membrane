QT += core gui qml
QT += 3dcore 3drender 3dinput 3dlogic 3dextras 3danimation
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11 console
CONFIG -= app_bundle

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        Phys/balloon.cpp \
        Phys/bead.cpp \
        Phys/bendingenergy.cpp \
        Phys/bendingparameters.cpp \
        Phys/membrane.cpp \
        Phys/membranefromobj.cpp \
        Phys/springforce.cpp \
        Phys/springforce2.cpp \
        Phys/tensor2.cpp \
        Phys/tether.cpp \
        Phys/vecd3d.cpp \
        main.cpp \
        mainwindow.cpp \
        q3d/beadinfo.cpp \
        q3d/edge.cpp \
        q3d/plane.cpp \
        q3d/q3dwidget.cpp \
        q3d/triangle.cpp \
        q3d/wavefrontobj.cpp \
        utils/mthread.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    Phys/CTS.h \
    Phys/balloon.h \
    Phys/bead.h \
    Phys/bendingenergy.h \
    Phys/bendingparameters.h \
    Phys/membrane.h \
    Phys/membranefromobj.h \
    Phys/springforce.h \
    Phys/springforce2.h \
    Phys/surfacewithphysics.h \
    Phys/tensor2.h \
    Phys/tether.h \
    Phys/vecd3d.h \
    mainwindow.h \
    q3d/beadinfo.h \
    q3d/edge.h \
    q3d/plane.h \
    q3d/q3dwidget.h \
    q3d/triangle.h \
    q3d/wavefrontobj.h \
    utils/classes.h \
    utils/mthread.h

FORMS += \
    mainwindow.ui

RESOURCES += \
    wireframe.qrc

DISTFILES += \
    WireframeEffect2.qml


INCLUDEPATH +="C:/Users/sm2983/Downloads/gsl-latest.tar/gsl-latest/install/include"
LIBS += -L"C:/Users/sm2983/Downloads/gsl-latest.tar/gsl-latest/install/lib/" -llibgsl -llibgslcblas
