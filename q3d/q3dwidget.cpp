#include "q3dwidget.h"


Qt3DRender::QGeometry * Q3dWidget::makeTri(Qt3DCore::QTransform * trans, const QColor& color, Qt3DCore::QEntity *_rootEntity)
{
    auto *lineEntity = new Qt3DCore::QEntity(_rootEntity);
    auto *geometry = new Qt3DRender::QGeometry(lineEntity);

    // position vertices (start and end)
    QByteArray bufferBytes;
    bufferBytes.resize(3 * 3 * sizeof(float)); // start.x, start.y, start.end + end.x, end.y, end.z
    float *positions = reinterpret_cast<float*>(bufferBytes.data());
    *positions++ = 0;
    *positions++ = 0;
    *positions++ = 0;
    *positions++ = 0;
    *positions++ = 0;
    *positions++ = 1;
    *positions++ = 1.0;
    *positions++ = 1.0;
    *positions++ = 0.5;

    Qt3DRender::QBuffer *triBuff = new Qt3DRender::QBuffer(geometry);
    triBuff->setData(bufferBytes);

    auto *positionAttribute = new Qt3DRender::QAttribute(geometry);
    positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
    positionAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
    positionAttribute->setVertexSize(3);
    positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    positionAttribute->setBuffer(triBuff);
    positionAttribute->setByteStride(3 * sizeof(float));
    positionAttribute->setCount(3);
    geometry->addAttribute(positionAttribute); // We add the vertices in the geometry

    // connectivity between vertices
    QByteArray indexBytes;
    indexBytes.resize(3 * sizeof(unsigned int)); // start to end
    unsigned int *indices = reinterpret_cast<unsigned int*>(indexBytes.data());
    *indices++ = 0;
    *indices++ = 1;
    *indices++ = 2;
    *indices++ = 0;

    auto *indexBuffer = new Qt3DRender::QBuffer(geometry);
    indexBuffer->setData(indexBytes);

    auto *indexAttribute = new Qt3DRender::QAttribute(geometry);
    indexAttribute->setVertexBaseType(Qt3DRender::QAttribute::UnsignedInt);
    indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
    indexAttribute->setBuffer(indexBuffer);
    indexAttribute->setCount(3);
    geometry->addAttribute(indexAttribute); // We add the indices linking the points in the geometry

    // mesh
    auto *line = new Qt3DRender::QGeometryRenderer(lineEntity);
    line->setGeometry(geometry);
    line->setPrimitiveType(Qt3DRender::QGeometryRenderer::LineLoop);
    auto *material = new Qt3DExtras::QPhongMaterial(lineEntity);
    material->setAmbient(color);

    // entity

    lineEntity->addComponent(line);
    lineEntity->addComponent(trans);
    lineEntity->addComponent(material);

    return geometry;
}

Q3dWidget::Q3dWidget(QWidget *parent ):QWidget(parent)
{


}
void Q3dWidget::init(){

    this->_view=new Qt3DExtras::Qt3DWindow();

    Qt3DCore::QEntity *rootEntity = new Qt3DCore::QEntity;

    Qt3DRender::QCamera *camera = _view->camera();
    //_plan=new Geometry::Plane(91,91,10,10);
    //_memb=new Physics::Membrane(0.01,1000,_plan);
    //_balloon=new Physics::Balloon(1e-2);



//    _mesh=new Physics::MembraneFromObj(1e-1);
    _mesh=new Physics::Balloon(1e-1);

    // Camera
    //    camera->lens()->setPerspectiveProjection(45.0f, 16.0f/9.0f, 0.1f, 1000.0f);
    camera->setPosition(QVector3D(0, 0, 40.0f));
    camera->setViewCenter(QVector3D(0, 0, 0));

    /*
    Qt3DCore::QEntity *lightEntity = new Qt3DCore::QEntity(rootEntity);

    Qt3DCore::QEntity *lightEntity2 = new Qt3DCore::QEntity(lightEntity);
    Qt3DCore::QTransform *lightTransform = new Qt3DCore::QTransform();
    Qt3DRender::QPointLight *light = new Qt3DRender::QPointLight(lightEntity2);
    lightTransform->setTranslation(QVector3D(0, -200, 10.0f));
    light->setColor("white");
    light->setIntensity(0.2);
    lightEntity2->addComponent(light);
    lightEntity2->addComponent(lightTransform);

    lightEntity2 = new Qt3DCore::QEntity(lightEntity);
    light = new Qt3DRender::QPointLight(lightEntity2);
    lightTransform = new Qt3DCore::QTransform();
    lightTransform->setTranslation(QVector3D(0, 100, 10.0f));
    light->setColor("white");
    light->setIntensity(0.75);
    lightEntity2->addComponent(light);
    lightEntity2->addComponent(lightTransform);


    lightEntity->addComponent(camera->transform());
*/


    Qt3DCore::QEntity *sceneEntity = createScene();



    Qt3DCore::QTransform *sceneTrans=new Qt3DCore::QTransform();
    sceneEntity->addComponent(sceneTrans);
    sceneEntity->setParent(rootEntity);

    // For camera controls

    Qt3DExtras::QOrbitCameraController *camController = new Qt3DExtras::QOrbitCameraController(sceneTrans);

    camController->setLinearSpeed( 1000.0f );
    camController->setLookSpeed( 100.0f );

    camController->setCamera(camera);






    Qt3DRender::QRenderSurfaceSelector *surfaceSelector = new Qt3DRender::QRenderSurfaceSelector;
    surfaceSelector->setSurface(_view);
    Qt3DRender::QViewport *viewport = new Qt3DRender::QViewport(surfaceSelector);
    viewport->setNormalizedRect(QRectF(0, 0, 1.0, 1.0));
    Qt3DRender::QCameraSelector *cameraSelector = new Qt3DRender::QCameraSelector(viewport);
    cameraSelector->setCamera(camera);
    Qt3DRender::QClearBuffers *clearBuffers = new Qt3DRender::QClearBuffers(cameraSelector);
    clearBuffers->setBuffers(Qt3DRender::QClearBuffers::ColorDepthBuffer);
    clearBuffers->setClearColor(Qt::lightGray);
    Qt3DRender::QRenderStateSet *renderStateSet = new Qt3DRender::QRenderStateSet(clearBuffers);
    Qt3DRender::QCullFace *cullFace = new Qt3DRender::QCullFace(renderStateSet);
    cullFace->setMode(Qt3DRender::QCullFace::NoCulling);
    renderStateSet->addRenderState(cullFace);
    Qt3DRender::QClearBuffers *clB = new Qt3DRender::QClearBuffers;
    clB->setClearColor(QColor(0,255,255));

    Qt3DRender::QDepthTest *depthTest = new Qt3DRender::QDepthTest;
    depthTest->setDepthFunction(Qt3DRender::QDepthTest::Less);
    renderStateSet->addRenderState(depthTest);

    _view->setActiveFrameGraph(surfaceSelector);
    _view->defaultFrameGraph()->setClearColor(QColor::fromRgbF(0.0, 0.5, 1.0, 1.0));

    _view->setRootEntity(rootEntity);

    _view->defaultFrameGraph()->setClearColor(QColor(QRgb(0x4d4d4f)));

    _view->renderSettings()->setRenderPolicy(Qt3DRender::QRenderSettings::OnDemand);


    _container = QWidget::createWindowContainer(this->_view);

    QHBoxLayout *hLayout = new QHBoxLayout(this);
    QVBoxLayout *vLayout = new QVBoxLayout();
    vLayout->setAlignment(Qt::AlignTop);
    hLayout->addWidget(_container, 1);
    hLayout->addLayout(vLayout);
    resetSnapshotID();
}
QString Q3dWidget::getSnapshotPath(){
    return _templatePath.arg(_snapshotID,  6, 'g', -1, '0');
}
void Q3dWidget::resetSnapshotID(){
    QString dir=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_%1_%2)").arg(_mesh->getShape(),_mesh->getTitle());
     QDir dirInfo(dir);
     if(!dirInfo.exists()){
         QDir().mkdir(dir);
     }
    _templatePath=QString(R"(%1\f%2.png)").arg(dir);
    _snapshotID=0;
    while(true){
        if(!QFile::exists(getSnapshotPath())){
            break;
        }
        _snapshotID++;
    }
}
Qt3DRender::QMaterial* createMat(Qt3DCore::QEntity *root ){
    QFile* file=new QFile();
    file->setFileName(":/shaders/robustwireframe.frag");
    file->open(QIODevice::ReadOnly);
    QByteArray  frD=file->readAll();
    file->close();
    file->setFileName(":/shaders/robustwireframe.geom");
    file->open(QIODevice::ReadOnly);

    QByteArray  geomD=file->readAll();
    file->close();
    file->setFileName(":/shaders/robustwireframe.vert");
    file->open(QIODevice::ReadOnly);
    QByteArray  vertD=file->readAll();
    file->close();
    Qt3DRender::QMaterial *ret=new Qt3DRender::QMaterial(root);
    Qt3DRender::QEffect *eff=new Qt3DRender::QEffect(ret);
    eff->addParameter(new Qt3DRender::QParameter("ka", QVariant(QVector4D(0.1, 0.1, 0.1,1))));
    eff->addParameter(new Qt3DRender::QParameter("kd", QVariant(QVector4D(0.7, 0.7, 0.7,1))));
    eff->addParameter(new Qt3DRender::QParameter("ks", QVariant(QVector4D(0.95, 0.95, 0.95,1))));
    eff->addParameter(new Qt3DRender::QParameter("shininess", QVariant(1000.0)));

    Qt3DRender::QTechnique *tec=new Qt3DRender::QTechnique(eff);
    Qt3DRender::QRenderPass *qrpass=new  Qt3DRender::QRenderPass(tec);


    auto rs=new Qt3DRender::QBlendEquation();

    auto m_blendState=new Qt3DRender::QBlendEquationArguments();
    m_blendState->setSourceRgb(Qt3DRender::QBlendEquationArguments::SourceAlpha);
    m_blendState->setDestinationRgb(Qt3DRender::QBlendEquationArguments::OneMinusSourceAlpha);
    rs->setBlendFunction(Qt3DRender::QBlendEquation::Add);

    qrpass->addRenderState(m_blendState);
    qrpass->addRenderState(rs);


    auto depth=new Qt3DRender::QDepthTest();
    depth->setDepthFunction(Qt3DRender::QDepthTest::Less);
    qrpass->addRenderState(depth);

    Qt3DRender::QShaderProgram *qsp=new Qt3DRender::QShaderProgram(qrpass);
    qsp->setVertexShaderCode(vertD);
    qsp->setFragmentShaderCode(frD);
    qsp->setGeometryShaderCode(geomD);

    qrpass->setShaderProgram(qsp);
    tec->addRenderPass(qrpass);
    tec->addParameter(new Qt3DRender::QParameter(QString("light.position"),QVariant(QVector3D(30.0, 0.0, 30.0))));
    tec->addParameter(new Qt3DRender::QParameter(QString("light.intensity"),QVariant(QVector3D(2,2,2))));

    eff->addTechnique(tec);

    Qt3DRender::QFilterKey *m_filterKey = new Qt3DRender::QFilterKey(ret);
    m_filterKey->setName(QStringLiteral("renderingStyle"));
    m_filterKey->setValue(QStringLiteral("forward"));

    // Add the pass to the technique


    // Set the targeted GL version for the technique
    tec->graphicsApiFilter()->setApi(Qt3DRender::QGraphicsApiFilter::OpenGL);
    tec->graphicsApiFilter()->setMajorVersion(3);
    tec->graphicsApiFilter()->setMinorVersion(2);
    tec->graphicsApiFilter()->setProfile(Qt3DRender::QGraphicsApiFilter::CoreProfile);

    // Add filter
    tec->addFilterKey(m_filterKey);

    // Add the technique to the effect


    // Set the effect on the materials
    ret->setEffect(eff);

    /*
    QQmlApplicationEngine eng;
    eng.load();
    auto objs=eng.rootObjects();
    assert(objs.size()!=0);
    Qt3DRender::QMaterial* qe=qobject_cast<Qt3DRender::QMaterial*>(objs[0]);

    //qe=new Qt3DExtras::QPhongMaterial(root);*/
    return  ret;

}
Qt3DCore::QEntity * Q3dWidget::createScene(){
    Qt3DCore::QEntity *rootEntity = new Qt3DCore::QEntity;
    Qt3DRender::QMaterial *material=createMat(rootEntity);
    //Qt3DRender::QMaterial *material = new Qt3DRender::QMaterial(rootEntity);
    material->addParameter(new Qt3DRender::QParameter(QString("ka"),QVariant(QVector4D(0.05,0.05,0.05,1))));
    material->addParameter(new Qt3DRender::QParameter(QString("kd"),QVariant(QVector4D(0.2,0.2,0.2,1 ))));
    material->addParameter(new Qt3DRender::QParameter(QString("ksp"),QVariant(QVector4D(0.8,0.8,0.8,1))));

    material->addParameter(new Qt3DRender::QParameter(QString("line.width"),QVariant(0.5)));
    material->addParameter(new Qt3DRender::QParameter(QString("line.color"),QVariant(QVector4D(0.0,65,0.29,1))));
    material->addParameter(new Qt3DRender::QParameter(QString("shininess"),QVariant(300)));
    material->addParameter(new Qt3DRender::QParameter(QString("alpha"),QVariant(1)));


    /*
    for(auto p:material->parameters()){
        std::cout<<"--"<<p->name().toStdString()<<std::endl;
    }*/
    //material->addParameter()



    Qt3DCore::QEntity *customMeshEntity = new Qt3DCore::QEntity(rootEntity);
    // Transform
    Qt3DCore::QTransform *transform = new Qt3DCore::QTransform;
    transform->setScale(8.0f);
    // Custom Mesh (TetraHedron)
    //Qt3DRender::QGeometryRenderer *customMeshRenderer=_plan->mesh();
    //Qt3DRender::QGeometryRenderer *customMeshRenderer=_balloon->mesh();
    Qt3DRender::QGeometryRenderer *customMeshRenderer=_mesh->mesh();
    customMeshEntity->addComponent(customMeshRenderer);
    customMeshEntity->addComponent(transform);
    customMeshEntity->addComponent(material);
    //fixme
    /*
    _triGem=makeTri(transform,QColor(255,0,0),rootEntity);
    _triGem2=makeTri(transform,QColor(0,0,255),rootEntity);
    Qt3DCore::QEntity * parent=(Qt3DCore::QEntity *)_triGem2->parent();
    parent->setEnabled(false);

    drawTriangle(1,_triGem);
    drawTriangle(1,_triGem2);
*/

    return rootEntity;
}
void Q3dWidget::snapshot(){
    QScreen *screen = _view->screen();
    screen->grabWindow(_view->winId()).save(getSnapshotPath());
    _snapshotID++;
    //_container->grab().save(_templatePath.arg(_snapshotID++));
}
void Q3dWidget::tupdate(int t){

    /*
    if(QThread::currentThread()!=QApplication::instance()->thread()){
        QMetaObject::invokeMethod(this,         // obj
                                  "tupdate",         // member: don't put parameters
                                  Qt::QueuedConnection,     // connection type
                                  Q_ARG(int, t));     // val1
    }else{
        for(int i=0;i<_numVertecies;i++){
            _rawVertexArray[i*9]+=0.01;
        }
    }*/
}
void Q3dWidget::drawTriangle(int ind,Qt3DRender::QGeometry* tri){
    /*
    auto attrs=tri->attributes();

    Qt3DRender::QBuffer *triBuff;
    QString str=Qt3DRender::QAttribute::defaultPositionAttributeName();
    for(auto att:attrs){
        if(att->name()==str){
            triBuff=att->buffer();
        }
    }

    QByteArray bufferBytes=triBuff->data();
    float *positions = reinterpret_cast<float*>(bufferBytes.data());
    float *triPos=_plan->getTri(ind);
    for(int i=0;i<9;i++){
        positions[i]=triPos[i];
    }
    triBuff->updateData(0,bufferBytes);*/

}
void Q3dWidget::lookAtCenter(){
    _view->camera()->setViewCenter(QVector3D(0, 0, 0));
}

