#include "plane.h"
#define BEADGAMMA 1
Geometry::Plane::Plane(int rows,int columns,double w,double h):_rows(rows),_columns(columns),_w(w),_h(h)
{ 
_dw=_w/(double)(_columns-1);
_dh=_h/(double)(_rows-1);
    /*
        for(int j=0;j<len;j++){
            //delete _v[j];
        }
        //delete _v;
    */

    updatable=true;
    _membrn=nullptr;
    _tris=QVector<Triangle*>();
    _beads=new QVector<Physics::BeadInfo*>();
    _edges=new QVector<Geometry::Edge*>();
    _bb=new QVector<Geometry::Edge*>();
    _bt=new QVector<Geometry::Edge*>();
    _bl=new QVector<Geometry::Edge*>();
    _br=new QVector<Geometry::Edge*>();
    _rand=QRandomGenerator::global();
    createPlan();
    rebuild();
#ifdef QT_Debug
    int c=0;
    for(int i=0;i<_edges;i++){
        Geometry::Edge *e=_edges->at(i);
        assert(e->_tris[0]!=nullptr);
        if(e->_tris[1]==nullptr){
            c++;
        }
    }
    assert(c==(rows+columns)*2);
#endif

}
void Geometry::Plane::rebuild(){

    updateVertecies();
    updateIndecies();
    _mesh->setVertexCount(_tris.size()*3);
}
void Geometry::Plane::createTriangleByIndecies(int a,int b,int c,double curv,Geometry::Triangle*& t){
    t=new Geometry::Triangle(this,_tris.size(), _beads->at(a),_beads->at(b),_beads->at(c),curv);
}
Geometry::Edge* Geometry::Plane::createEdgeByIndecies(int a,int b){
    return new Geometry::Edge(this,_beads->at(a),_beads->at(b));
}
void Geometry::Plane::createPlan(){

    _mesh = new Qt3DRender::QGeometryRenderer;
    Qt3DRender::QGeometry *customGeometry = new Qt3DRender::QGeometry(_mesh);
    _verBuff = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::VertexBuffer, customGeometry);
    _indBuff = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::IndexBuffer, customGeometry);

    // vec3 for position
    // vec3 for colors
    // vec3 for normals

    // 4 distinct vertices
    int numVer=_rows*_columns;
    QByteArray verBufferData;
    verBufferData.resize(numVer * (3 + 3 +3) * sizeof(float));
    QVector<QVector3D> vertices = QVector<QVector3D>();
    QVector3D red(1.0,0,0);
    //float dw=_w/(float)(_columns-1),dh=_h/(float)(_rows-1);

    for(int j=0;j<_rows;j++){
        for(int i=0;i<_columns;i++){
            double x=(_rand->generateDouble()-0.5)*0.001;
            double y=(_rand->generateDouble()-0.5)*0.001;
            double z=(_rand->generateDouble())*0.0;
            auto vec=QVector3D((i-(_columns-1)*0.5)*_dw+x,(j-(_rows-1)*0.5)*_dh+y,z);
            vertices<<vec <<QVector3D(0,0,1)<<red;

            _beads->append(new Physics::BeadInfo(this,new Physics::VecD3d(&vec),BEADGAMMA,_beads->size()));
        }
    }
    auto rawVertexArray = reinterpret_cast<float *>(verBufferData.data());
    int idx = 0;
    for(auto v:vertices){
        rawVertexArray[idx++]=v.x();
        rawVertexArray[idx++]=v.y();
        rawVertexArray[idx++]=v.z();
    }

    /*
     * QByteArray indexBufferData;
    QVector<QVector3D> vertices2 = QVector<QVector3D>();
    vertices2<<QVector3D(0,0,0) <<QVector3D(0,0,1)<<red;
    vertices2<<QVector3D(0,0.1,0) <<QVector3D(0,0,1)<<red;
    vertices2<<QVector3D(0.1,0,0) <<QVector3D(0,0,1)<<red;
    indexBufferData.resize(3 * sizeof(ushort));
    ushort *rawIndexArray = reinterpret_cast<ushort *>(indexBufferData.data());

    rawIndexArray[0]=0;
    rawIndexArray[1]=1;
    rawIndexArray[2]=2;
    */

    int numFaces=(_rows-1)*(_columns-1) ;


    // Indices (12)

    QByteArray indexBufferData;

    indexBufferData.resize(numFaces*6 * sizeof(ushort));

    ushort *rawIndexArray = reinterpret_cast<ushort *>(indexBufferData.data());
    for(int i=0;i<numFaces*6;i++){
        rawIndexArray[i]=0;
    }
    int index=0;
    for(int j=0;j<_rows-1;j++){
        Triangle *lt=nullptr;
        for(int i=0;i<_columns-1;i++){

            int col=(j*_columns);

            Geometry::Triangle *t1=nullptr,*t2=nullptr;
            int id=_tris.size();
            if(j==0){
                //t1=new Geometry::Triangle(this,col+ i,col+ i+1,col+i+_columns+1);
                createTriangleByIndecies(col+ i,col+ i+1,col+i+_columns+1,0,t1);
                assert(t1!=nullptr);

                if(i==0){
                    //t2=new Geometry::Triangle(this,new Geometry::Edge(this,_columns,0),t1->_e[2],new Geometry::Edge(this,_columns+1,_columns));
                    t2=new Geometry::Triangle(this,id+1,createEdgeByIndecies(_columns,0),t1->_e[2],createEdgeByIndecies(_columns+1,_columns),0);
                }else{
                    t2=new Geometry::Triangle(this,id+1,lt->_e[1],t1->_e[2],createEdgeByIndecies(i+_columns+1,i+_columns),0);
                }

                _tris<<t1<<t2;
            }else{
                double c=(j==_rows-2||i==0||i==_columns-2)?0:1;

                Triangle *lrt=_tris.at((j-1)*(_columns-1)*2+2*i+1);
                t1=new Triangle(this,id,lrt->_e[2] ,createEdgeByIndecies(col+ i+1,col+i+_columns+1),createEdgeByIndecies(col+i+_columns+1,col+i),c);


                if(i==0){
                    t2=new Triangle(this,id+1,createEdgeByIndecies(col+_columns,col),t1->_e[2],createEdgeByIndecies(col+i+_columns+1,col+i+_columns),c);

                }else{

                    t2=new Triangle(this,id+1,lt->_e[1],t1->_e[2],createEdgeByIndecies(col+i+_columns+1,col+i+_columns),c);
                    //t2=new Triangle(t1->_e[2],new Edge((j+1)*_columns+1,(j+1)*_columns),lt->_e[1]);
                }
                _tris<<t1<<t2;

            }
            lt=t1;
            if(i==0){
                _bl->append(t2->_e[0]);
            }
            if(j==0){
                _bb->append(t1->_e[0]);
            }
            if(i==_columns-2){
                _br->append(t1->_e[1]);
            }
            if(j==_rows-2){
                _bt->append(t2->_e[2]);
            }

            /*
            int i2=(i+j*(_columns));
            assert(i2<numVer);
            assert(t2!=nullptr);*/

            /*real
            rawIndexArray[6*index]=i2;
            rawIndexArray[6*index+1]=i2+1;
            rawIndexArray[6*index+2]=i2+_columns+1;
            //------
            rawIndexArray[6*index+3]=i2+_columns;
            rawIndexArray[6*index+4]=i2;
            rawIndexArray[6*index+5]=i2+_columns+1;
*/

            rawIndexArray[6*index]=t1->_v[0]->ID;

            rawIndexArray[6*index+1]=t1->_v[1]->ID;
            rawIndexArray[6*index+2]=t1->_v[2]->ID;
            //------

            rawIndexArray[6*index+3]=t2->_v[0]->ID;
            rawIndexArray[6*index+4]=t2->_v[1]->ID;
            rawIndexArray[6*index+5]=t2->_v[2]->ID;

            index++;
        }
    }
    //FIXME
    /*
    for(auto bb : *_beads){
        bb->orderConnections();
    }*/
    _verBuff->setData(verBufferData);
    _indBuff->setData(indexBufferData);

    // Attributes
    Qt3DRender::QAttribute *_positionAttribute = new Qt3DRender::QAttribute();
    _positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    _positionAttribute->setBuffer(_verBuff);
    _positionAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
    _positionAttribute->setVertexSize(3);
    _positionAttribute->setByteOffset(0);
    _positionAttribute->setByteStride(9 * sizeof(float));
    _positionAttribute->setCount(numVer);

    _positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());

    Qt3DRender::QAttribute *normalAttribute = new Qt3DRender::QAttribute();
    normalAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    normalAttribute->setBuffer(_verBuff);
    normalAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
    normalAttribute->setVertexSize(3);
    normalAttribute->setByteOffset(3 * sizeof(float));
    normalAttribute->setByteStride(9 * sizeof(float));
    //normalAttribute->setCount(numVer);
    normalAttribute->setCount(3);
    normalAttribute->setName(Qt3DRender::QAttribute::defaultNormalAttributeName());

    Qt3DRender::QAttribute *colorAttribute = new Qt3DRender::QAttribute();
    colorAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    colorAttribute->setBuffer(_verBuff);
    colorAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
    colorAttribute->setVertexSize(3);
    colorAttribute->setByteOffset(6 * sizeof(float));
    colorAttribute->setByteStride(9 * sizeof(float));
    //    colorAttribute->setCount(numVer);
    colorAttribute->setCount(3);
    colorAttribute->setName(Qt3DRender::QAttribute::defaultColorAttributeName());

    Qt3DRender::QAttribute *indexAttribute = new Qt3DRender::QAttribute();
    indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
    indexAttribute->setBuffer(_indBuff);
    indexAttribute->setVertexBaseType(Qt3DRender::QAttribute::UnsignedShort);
    indexAttribute->setVertexSize(1);
    indexAttribute->setByteOffset(0);
    indexAttribute->setByteStride(0);
    //indexAttribute->setCount((_rows-1)*(_columns-1)*6);
    indexAttribute->setCount(3);

    customGeometry->addAttribute(_positionAttribute);
    customGeometry->addAttribute(normalAttribute);
    customGeometry->addAttribute(colorAttribute);
    customGeometry->addAttribute(indexAttribute);

    _mesh->setInstanceCount(1);
    _mesh->setIndexOffset(0);
    _mesh->setFirstInstance(0);
    _mesh->setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);

    _mesh->setGeometry(customGeometry);
    // 4 faces of 3 points
    _mesh->setVertexCount(_tris.size()*3);
    //_mesh->setVertexCount(3);
}

Qt3DRender::QGeometryRenderer * Geometry::Plane::mesh(){
    updateVIN();

    return _mesh;
}
float* Geometry::Plane::getTri(int i){

    auto tri=_tris.at(i);

    auto qba=_verBuff->data();
    auto *rawIndexArray = reinterpret_cast<float *>(qba.data());
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            _triPoses[i*3+j]=rawIndexArray[tri->_v[i]->ID*9+j];
        }
    }
    return _triPoses;
}
void Geometry::Plane::updateIndecies(){

    /*
    auto atts=mesh()->geometry()->attributes();
    Qt3DRender::QBuffer *indBuff=nullptr;
    for(auto att:atts){
        if(att->attributeType()==Qt3DRender::QAttribute::IndexAttribute){
            indBuff=att->buffer();
            break;
        }

    }
    assert(indBuff!=nullptr);*/
    auto bufferBytes=_indBuff->data();
    int size=(_tris.size()) * (3 ) * sizeof(ushort);
    if(bufferBytes.size()!=size){
        //bufferBytes=QByteArray();
        bufferBytes.resize(size);
    }
    ushort *rawIndexArray = reinterpret_cast<ushort *>(bufferBytes.data());
    for(int index=0;index<_tris.count();index++){
        auto tri=_tris.at(index);
        rawIndexArray[3*index]=tri->_v[0]->ID;
        rawIndexArray[3*index+1]=tri->_v[1]->ID;
        rawIndexArray[3*index+2]=tri->_v[2]->ID;
    }
    _indBuff->setData(bufferBytes);

}
void Geometry::Plane::updateNormals(){

    for(auto v:_normals){
        v->zero();
    }

    for(int i=_normals.size();i<_beads->size();i++){
        _normals.append(QSharedPointer<Physics::VecD3d>(new Physics::VecD3d()));
    }
    for(auto tri:_tris){
        auto norm=tri->_norm;
        _normals[tri->_v[0]->ID]->add(norm);
        _normals[tri->_v[1]->ID]->add(norm);
        _normals[tri->_v[2]->ID]->add(norm);
    }
    for(auto v:_normals){
        v->nomilize();
    }
}
void Geometry::Plane::updateVertecies(){
    /*
    auto atts=mesh()->geometry()->attributes();
    Qt3DRender::QBuffer *posBuff=nullptr;
    for(auto att:atts){
        QString str=Qt3DRender::QAttribute::defaultPositionAttributeName();
        if(att->name()==str){
            posBuff=att->buffer();
            break;
        }

    }
QByteArray bufferBytes=posBuff->data();
*/
    updateNormals();

    QByteArray bufferBytes=_verBuff->data();
    int size=_beads->size() * (3 + 3 +3) * sizeof(float);
    if(bufferBytes.size()!=size){
        //bufferBytes=QByteArray();
        //int i=bufferBytes.size();
        bufferBytes.resize(size);

    }
    assert(bufferBytes.size()==size);
    float *positions = reinterpret_cast<float*>(bufferBytes.data());

    for(int i=0;i<_beads->size();i++){
        auto b=_beads->at(i);
        auto v=b->_coords->_coords;
        auto n=_normals.at(b->ID);
        positions[i*9+0]=v[0];
        positions[i*9+1]=v[1];
        positions[i*9+2]=v[2];

        positions[i*9+3]=n->_coords[0];
        positions[i*9+4]=n->_coords[1];
        positions[i*9+5]=n->_coords[2];

        positions[i*9+6]=0;
        positions[i*9+7]=0;
        positions[i*9+8]=1;
    }
    _verBuff->setData(bufferBytes);


}
QVector<Physics::BeadInfo*>* Geometry::Plane::getBeads(){
    return _beads;
}
QVector<Geometry::Edge*>* Geometry::Plane::getEdges(){
    return _edges;
}

Geometry::Plane::~Plane(){
    /*
    for(int i=0;i<_tris.size();i++){
        delete _tris.at(i);
    }*/
}
/*Physics::Membrane* Geometry::Plane::getMembraneInstance(double dt,double k){
    if(_membrn==nullptr){
        _membrn=new Physics::Membrane(_columns,_rows,dt,k,this);
    }
    return  _membrn;
}*/
void Geometry::Plane:: updateTopology(){
    if(!updatable){
        return;
    }
    int n=_tris.size();

    int sel=_rand->generate()%n;
    int edg=_rand->generate()%3;

    auto t1=_tris[sel];
    auto e=t1->_e[edg];
    auto t2=e->_tris[0]==t1?e->_tris[1]:e->_tris[0];
    if(t2!=nullptr&&t2!=t1){
        swapEdge(t1,t2);
    }
    /*
    updateIndecies();
    updateNormals();*/
}
void Geometry::Plane::swapEdge(Geometry::Triangle* t1,Geometry::Triangle* t2){
    auto comEdge=t1->commonEdge(t2);
    /*
#ifdef QT_DEBUG
    assert(comEdge!=nullptr);
    assert(comEdge->triIndex(t1)!=-1);
    assert(comEdge->triIndex(t2)!=-1);
#endif*/
    Geometry::Edge* edge11=t1->edgeWithVertexExclude(comEdge->_vid1,comEdge);
    Geometry::Edge* edge12=t2->edgeWithVertexExclude(comEdge->_vid1,comEdge);

    Geometry::Edge* edge21=t1->edgeWithVertexExclude(comEdge->_vid2,comEdge);
    Geometry::Edge* edge22=t2->edgeWithVertexExclude(comEdge->_vid2,comEdge);
    /*
    Triangle *t11=edge11->_tris[1-edge11->triIndex(t1)];
    Triangle *t12=edge12->_tris[1-edge12->triIndex(t2)];
    Triangle *t21=edge21->_tris[1-edge21->triIndex(t1)];
    Triangle *t22=edge22->_tris[1-edge22->triIndex(t2)];
*/
    Physics::BeadInfo* vid1=edge11->_vid1==comEdge->_vid1?edge11->_vid2:edge11->_vid1;
    Physics::BeadInfo* vid2=edge12->_vid1==comEdge->_vid1?edge12->_vid2:edge12->_vid1;

    //    if(t11!=t12&&t21!=t22)
    //    {

    for(int i=0;i<_edges->size();i++){
        auto ee=_edges->at(i);
        if(ee->equal(vid1,vid2)){
            return;
        }
    }
    /*
#ifdef QT_DEBUG
        assert(edge12->triIndex(t2)!=-1);
        assert(edge21->triIndex(t1)!=-1);

                std::cout<<"DE1\t"<<ee->_tris[0]<<"\t"<<ee->_tris[1]<<"\t"<<i<<"\t"<<_edges->size()<<"\r\n";
                auto tt=ee->_tris[0];
                std::cout<<ee<<"\t"<<_tris.indexOf(tt)<<"\t"<<_tris.indexOf(t1)<<"\t"<<_tris.indexOf(t2)<<std::endl;
                /\*


                std::cout<<tt->commonEdge(t1)<<std::endl;
                std::cout<<tt->commonEdge(t2)<<std::endl;
                std::cout<<"1\r\n";
                std::cout<<tt->commonEdge(t11)<<std::endl;
                std::cout<<tt->commonEdge(t12)<<std::endl;
                std::cout<<t11->commonEdge(t12)<<std::endl;
                std::cout<<"2\r\n";
                std::cout<<tt->commonEdge(t21)<<std::endl;
                std::cout<<tt->commonEdge(t22)<<std::endl;
                std::cout<<t21->commonEdge(t22)<<std::endl;
                std::cout<<"3\r\n";
                tt=ee->_tris[1];
                if(tt!=nullptr){
                    std::cout<<tt->commonEdge(t1)<<std::endl;
                    std::cout<<tt->commonEdge(t2)<<std::endl;

                    std::cout<<"4\r\n";
                    std::cout<<tt->commonEdge(t11)<<std::endl;
                    std::cout<<tt->commonEdge(t12)<<std::endl;
                    std::cout<<t11->commonEdge(t12)<<std::endl;
                    std::cout<<"5\r\n";
                    std::cout<<tt->commonEdge(t21)<<std::endl;
                    std::cout<<tt->commonEdge(t22)<<std::endl;
                    std::cout<<t21->commonEdge(t22)<<std::endl;
                    std::cout<<"6\r\n";
                }
                /\*
                    t1->printEdges();
                    t2->printEdges();
                    if(t11!=nullptr){
                        t11->printEdges();
                    }
                    if(t12!=nullptr){
                        t12->printEdges();
                    }
                    if(t21!=nullptr){
                        t21->printEdges();
                    }
                    if(t22!=nullptr){
                        t22->printEdges();
                    }*\/

                //assert(0);
                updatable=false;
                return;
            }
        }

    auto v1=e->commonVertex(edge11);
    auto v2=e->commonVertex(edge12);
    auto v3=edge11->commonVertex(edge12);
    if(!(v1&&v2&&v3&&v1!=v2&&v1!=v3&&v2!=v3)){
        std::cout<<"T1\r\n";
        assert(0);
    }

    v1=e->commonVertex(edge22);
    v2=e->commonVertex(edge21);
    v3=edge21->commonVertex(edge22);
    if(!(v1&&v2&&v3&&v1!=v2&&v1!=v3&&v2!=v3)){
        std::cout<<"T2\r\n";
        assert(0);
    }
    */
    //#endif

    double dx=DR(vid1,vid2,0);
    double dy=DR(vid1,vid2,1);
    double dz=DR(vid1,vid2,2);
    double k1=1/std::sqrt(dx*dx+dy*dy+dz*dz);

    double dx2=DR(comEdge->_vid1,comEdge->_vid2,0);
    double dy2=DR(comEdge->_vid1,comEdge->_vid2,1);
    double dz2=DR(comEdge->_vid1,comEdge->_vid2,2);
    double k2=1/std::sqrt(dx2*dx2+dy2*dy2+dz2*dz2);

    if(_rand->generateDouble()<k1/(k1+k2)&&t1->_curvater!=0&&t2->_curvater!=0){
        comEdge->_vid1->_connections->removeOne(comEdge);
        comEdge->_vid2->_connections->removeOne(comEdge);
        vid1->_connections->append(comEdge);
        vid2->_connections->append(comEdge);
        comEdge->_vid1=vid1;
        comEdge->_vid2=vid2;

        //


        //comEdge->updateRL();
        t1->_e[0]=edge11;
        t1->_e[1]=edge12;
        t1->_e[2]=comEdge;

        edge12->_tris[edge12->triIndex(t2)]=t1;

        t2->_e[0]=edge21;
        t2->_e[1]=comEdge;
        t2->_e[2]=edge22;

        edge21->_tris[edge21->triIndex(t1)]=t2;
#ifdef QT_DEBUG
        assert(edge12->_tris[0]!=t2&&edge12->_tris[1]!=t2);
        assert(edge21->_tris[0]!=t1&&edge21->_tris[1]!=t1);
#endif

        t1->updateVIndecies();

        t2->updateVIndecies();
//FIXME
        /*
        vid1->orderConnections();
        vid2->orderConnections();*/

    }

}
void Geometry::Plane::addQuadAt(int c,int r){
    int bid=_beads->size();

    Geometry::Triangle *t1,*t2;
    int id=_tris.size();
    if(c==-1){
        Edge* le=_bl->at(r);
        double x=le->_vid2->_coords->_coords[0];
        double y=le->_vid2->_coords->_coords[1];


        if(r==0){
            Physics::BeadInfo *b2=new Physics::BeadInfo(this,new Physics::VecD3d(x-_dw,y+_dh,0),BEADGAMMA,bid+1),
                    *b1=new Physics::BeadInfo(this,new Physics::VecD3d(x-_dw,y,0),BEADGAMMA,bid);
            _beads->append(b1);
            _beads->append(b2);


            //assert(lt->_e[0]==lt->)
            t1=new Geometry::Triangle(this,id,new Edge(this,b1,le->_vid2),le,new Edge(this,le->_vid1,b1),0);
            t2=new Geometry::Triangle(this,id+1,new Edge(this,b2,b1),t1->_e[2],new Edge(this,t1->_e[2]->_vid1,b2),0);

            //_tris.remove(10);

            _bb->insert(0,t1->_e[0]);
            _bt->insert(0,t2->_e[2]);
            _bl->replace(r,t2->_e[0]);
        }else{
            Physics::BeadInfo *b=new Physics::BeadInfo(this,new Physics::VecD3d(x-_dw,y+_dh,0),BEADGAMMA,bid);
            _beads->append(b);

            Geometry::Edge *te=_bt->at(0);

            t1=new Geometry::Triangle(this,id,te,le,new Edge(this,le->_vid1,te->_vid2),0);
            t2=new Geometry::Triangle(this,id+1,new Edge(this,b,te->_vid2),t1->_e[2],new Edge(this,t1->_e[2]->_vid1,b),0);
            _bl->replace(r,t2->_e[0]);
            _bt->replace(0,t2->_e[2]);
            //t2=new Geometry::Triangle(this,new Edge(this,b2,b1),t1->_e[2],new Edge(this,t1->_e[2]->_vid1,b2));
        }
        _tris.append(t1);
        _tris.append(t2);

    }else if(c==_columns){
        Edge* re=_br->at(r);
        double x=re->_vid1->_coords->_coords[0];
        double y=re->_vid1->_coords->_coords[1];



        if(r==0){
            Physics::BeadInfo *b2=new Physics::BeadInfo(this,new Physics::VecD3d(x+_dw,y+_dh,0),BEADGAMMA,bid+1),
                    *b1=new Physics::BeadInfo(this,new Physics::VecD3d(x+_dw,y,0),BEADGAMMA,bid);
            _beads->append(b1);
            _beads->append(b2);
            t1=new Geometry::Triangle(this,id,new Edge(this,re->_vid1,b1),new Edge(this,b1,b2),new Edge(this,b2,re->_vid1),0);
            t2=new Geometry::Triangle(this,id+1,re,t1->_e[2],new Edge(this,b2,re->_vid2),0);

            _bb->append(t1->_e[0]);
            _bt->append(t2->_e[2]);
            _br->replace(r,t1->_e[1]);
        }else{
            Physics::BeadInfo *b=new Physics::BeadInfo(this,new Physics::VecD3d(x+_dw,y+_dh,0),BEADGAMMA,bid);
            _beads->append(b);
            Geometry::Edge *te=_bt->last();
            t1=new Geometry::Triangle(this,id,te,new Edge(this,te->_vid1,b), new Edge(this,b,te->_vid2),0);
            t2=new Geometry::Triangle(this,id+1,re,t1->_e[2],new Edge(this,b,re->_vid2),0);
            _br->replace(r,t1->_e[1]);
            _bt->replace(_bt->size()-1,t2->_e[2]);
            //t2=new Geometry::Triangle(this,new Edge(this,b2,b1),t1->_e[2],new Edge(this,t1->_e[2]->_vid1,b2));
        }
        _tris.append(t1);
        _tris.append(t2);
    }else if(r==-1){
        Edge* be=_bb->at(c);
        double x=be->_vid1->_coords->_coords[0];
        double y=be->_vid1->_coords->_coords[1];
        if(c==0){
            Physics::BeadInfo *b2=new Physics::BeadInfo(this,new Physics::VecD3d(x+_dw,y-_dh,0),BEADGAMMA,bid+1),
                    *b1=new Physics::BeadInfo(this,new Physics::VecD3d(x,y-_dh,0),BEADGAMMA,bid);
            _beads->append(b1);
            _beads->append(b2);
            t1=new Geometry::Triangle(this,id,new Edge(this,b1,b2),new Edge(this,b2,be->_vid2),new Edge(this,be->_vid2,b1),0);
            t2=new Geometry::Triangle(this,id+1,new Edge(this,be->_vid1,b1),t1->_e[2],be,0);
            _tris.append(t1);
            _tris.append(t2);

            _bl->insert(0,t2->_e[0]);
            _br->insert(0,t1->_e[1]);
            _bb->replace(c,t1->_e[0]);
        }else{
            Physics::BeadInfo *b=new Physics::BeadInfo(this,new Physics::VecD3d(x+_dw,y-_dh,0),BEADGAMMA,bid);
            _beads->append(b);
            Edge *re=_br->at(0);


            t1=new Geometry::Triangle(this,id,new Edge(this,re->_vid1,b),new Edge(this,b,be->_vid2),new Edge(this,be->_vid2,re->_vid1),0);
            t2=new Geometry::Triangle(this,id+1,re,t1->_e[2],be,0);
            _tris.append(t1);
            _tris.append(t2);

            //_bl->insert(0,t2->_e[0]);
            _br->replace(0,t1->_e[1]);
            _bb->replace(c,t1->_e[0]);
        }

    }else if(r==_rows){
        Edge* te=_bt->at(c);
        double x=te->_vid2->_coords->_coords[0];
        double y=te->_vid2->_coords->_coords[1];
        if(c==0){
            Physics::BeadInfo *b2=new Physics::BeadInfo(this,new Physics::VecD3d(x+_dw,y+_dh,0),BEADGAMMA,bid+1),
                    *b1=new Physics::BeadInfo(this,new Physics::VecD3d(x,y+_dh,0),BEADGAMMA,bid);
            _beads->append(b1);
            _beads->append(b2);
            t1=new Geometry::Triangle(this,id,te,new Edge(this,te->_vid1,b2),new Edge(this,b2,te->_vid2),0);
            t2=new Geometry::Triangle(this,id+1,new Edge(this,b1,te->_vid2),t1->_e[2],new Edge(this,b2,b1),0);
            _tris.append(t1);
            _tris.append(t2);

            _bl->append(t2->_e[0]);
            _br->append(t1->_e[1]);
            _bt->replace(c,t2->_e[2]);
        }else{
            Physics::BeadInfo *b=new Physics::BeadInfo(this,new Physics::VecD3d(x+_dw,y+_dh,0),BEADGAMMA,bid);
            _beads->append(b);
            Edge *re=_br->last();
            t1=new Geometry::Triangle(this,id,te,new Edge(this,te->_vid1,b),new Edge(this,b,te->_vid2),0);
            t2=new Geometry::Triangle(this,id+1,re,t1->_e[2],new Edge(this,b,re->_vid2),0);
            _tris.append(t1);
            _tris.append(t2);


            _br->replace(_bl->size()-1,t1->_e[1]);
            _bt->replace(c,t2->_e[2]);
        }

    }
}

void Geometry::Plane::expand(){
    for(int i=0;i<_rows-1;i++){
        addQuadAt(-1,i);
    }
    _columns++;
    for(int i=0;i<_rows-1;i++){
        addQuadAt(_columns,i);
    }
    _columns++;
    for(int i=0;i<_columns-1;i++){
        addQuadAt(i,-1);

    }
    _rows++;
    for(int i=0;i<_columns-1;i++){
        addQuadAt(i,_rows);

    }
    _rows++;
}
void Geometry::Plane::addNewEdge(Edge* e){
    _edges->append(e);
}
void Geometry::Plane::updateVIN(){
    updateVertecies();
    updateIndecies();
    updateNormals();

}
