#ifndef Q3DWIDGET_H
#define Q3DWIDGET_H

#include <QWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QVBoxLayout>
#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DCore/QEntity>
#include <Qt3DRender/QCamera>
#include <Qt3DRender/QCameraLens>
#include <Qt3DCore/QTransform>
#include <Qt3DCore/QAspectEngine>

#include <Qt3DInput/QInputAspect>
#include <Qt3DCore/QEntity>
#include <Qt3DRender/QCamera>
#include <Qt3DRender/QCameraLens>
#include <Qt3DCore/QTransform>
#include <Qt3DCore/QAspectEngine>
#include <Qt3DRender/QRenderAspect>
#include <Qt3DExtras/QForwardRenderer>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/QCylinderMesh>
#include <Qt3DExtras/QSphereMesh>
#include <Qt3DExtras/QTorusMesh>
#include <Qt3DExtras/QOrbitCameraController>
#include <Qt3DExtras/QFirstPersonCameraController>
#include <Qt3DRender/QPointLight>
#include <QTransform>
#include <Qt3DInput/QInputAspect>
#include <Qt3DRender/QRenderStateSet>
#include <Qt3DRender/QRenderAspect>
#include <Qt3DExtras/QForwardRenderer>
#include <Qt3DExtras/QPerVertexColorMaterial>
#include <Qt3DRender/QGeometryRenderer>
#include <Qt3DRender/QGeometry>
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QBuffer>
#include <Qt3DRender/QRenderSurfaceSelector>
#include <Qt3DRender/QViewport>
#include <Qt3DRender/QCameraSelector>
#include <Qt3DRender/QClearBuffers>
#include <Qt3DRender/QCullFace>
#include <Qt3DRender/QDepthTest>
//#include <Qt3DExtras/M
#include <QApplication>
#include <q3d/plane.h>
#include <Qt3DRender/QParameter>
#include <Qt3DRender/QRenderPass>
#include <Qt3DRender/QShaderProgram>
#include <Qt3DRender/QEffect>
#include <Qt3DRender/QTechnique>
#include <Qt3DRender/QGraphicsApiFilter>
#include <QQmlApplicationEngine>
#include <QFile>
#include <QResource>
#include <QStringLiteral>
#include <Qt3DRender/QBlendEquation>
#include <Qt3DRender/QBlendEquationArguments>
#include <Qt3DRender/QRenderSettings>
#include <QScreen>
#include <QFileInfo>
#include "utils/mthread.h"
#include "Phys/balloon.h"
#include "Phys/membrane.h"
#include "Phys/membranefromobj.h"
#include <math.h>
#include <QDir>

struct IndirectElementDrawBuffer{ // Element Indirect
    uint count;
    uint instancesCount;
    uint firstIndex;
    uint baseVertex;
    uint baseInstance;
};
struct IndirectArrayDrawBuffer{ // Array Indirect
    uint count;
    uint instancesCount;
    uint first;
    uint baseInstance;
};

class Q3dWidget:public QWidget
{
private:
QWidget * _container;
int _snapshotID;
void resetSnapshotID();
QString _templatePath;
QString getSnapshotPath();
public:

    Physics::SurfaceWithPhysics *_mesh;
    Qt3DRender::QGeometry * _triGem,*_triGem2;
    Q_INVOKABLE void tupdate(int);
    Qt3DExtras::Qt3DWindow *_view;
    Q3dWidget(QWidget *parent = nullptr);
    Qt3DCore::QEntity * createScene();
    Qt3DRender::QGeometryRenderer* createPlan(int rows,int columns,double dw,double dh);
    void drawTriangle(int i,Qt3DRender::QGeometry*);
    Qt3DRender::QGeometry* makeTri(Qt3DCore::QTransform *, const QColor& color, Qt3DCore::QEntity *_rootEntity);
    void snapshot();
    void lookAtCenter();
    void init();
};

#endif // Q3DWIDGET_H
