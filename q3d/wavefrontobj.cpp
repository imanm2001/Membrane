#include "wavefrontobj.h"
#define BEADGAMMA 1
Geometry::WaveFrontObj::WaveFrontObj(QString path)
{
    _tris=new QVector<Triangle*>();
    _beads=new QVector<Physics::BeadInfo*>();
    _edges=new QVector<Geometry::Edge*>();

    _mesh = new Qt3DRender::QGeometryRenderer;
    Qt3DRender::QGeometry *customGeometry = new Qt3DRender::QGeometry(_mesh);
    _verBuff = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::VertexBuffer, customGeometry);
    _indBuff = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::IndexBuffer, customGeometry);

    //auto s=QString(R"(C:\Users\sm2983\Documents\Projects\Membrane\Results\Shape\Shape_force_%1.txt)").arg(_title);
    auto f=new QFile(path);

    if(f->exists()&&f->open(QIODevice::ReadOnly | QIODevice::Text)){

        auto fins=new QTextStream(f);
        char c;
        QChar comment=QChar('#');
        auto wordReg=QRegExp("[\\s^/]+");
        while(!fins->atEnd()){
            auto line=fins->readLine();
            if(line.length()>0){
                if(line.at(0)!=comment){
                    auto array=line.split(wordReg, Qt::SkipEmptyParts);
                    auto command=array.at(0);
                    if(command=="v"){
                        double x=array.at(1).toDouble();
                        double y=array.at(2).toDouble();
                        double z=array.at(3).toDouble();
                        auto *b=new Physics::BeadInfo(this,new Physics::VecD3d(x,y,z),BEADGAMMA,_beads->size());
                        _beads->append(b);
                    }else if(command=="f"){
                        int ind1,ind2,ind3;
                        if(array.size()==4){
                            ind1=array[1].toInt()-1;
                            ind2=array[2].toInt()-1;
                            ind3=array[3].toInt()-1;
                        }else if(array.size()==7){
                            ind1=array[1].toInt()-1;
                            ind2=array[3].toInt()-1;
                            ind3=array[5].toInt()-1;
                        }else{
                            std::cout<<line.toStdString()<<std::endl;
                            assert(0);
                        }

                        auto b1=_beads->at(ind1);
                        auto b2=_beads->at(ind2);
                        auto b3=_beads->at(ind3);
                        auto *tri=new Geometry::Triangle(this,_tris->size(),getEdge(b1,b2),getEdge(b2,b3),getEdge(b3,b1),0);
                        tri->getNormal();
                        _tris->append(tri);
                    }
                }

            }
        }



    buildMesh(customGeometry);
    }
}
Geometry::Edge* Geometry::WaveFrontObj::getEdge(Physics::BeadInfo* b1,Physics::BeadInfo* b2){
    Geometry::Edge* ret=nullptr;
    for(int i=0;i<_edges->size();i++){
        auto e=_edges->at(i);
        if(e->equal(b1,b2)){
            ret=e;
            break;
        }
    }
    if(!ret){
        ret=new Geometry::Edge(this,b1,b2);
        _edges->append(ret);
    }
    return ret;
}
void Geometry::WaveFrontObj::updateVertecies(){

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
        positions[i*9+0]=v[0]/3.0;
        positions[i*9+1]=v[1]/3.0;
        positions[i*9+2]=v[2]/3.0;

        positions[i*9+3]=n->_coords[0];
        positions[i*9+4]=n->_coords[1];
        positions[i*9+5]=n->_coords[2];

        positions[i*9+6]=0;
        positions[i*9+7]=0;
        positions[i*9+8]=1;
    }
    _verBuff->setData(bufferBytes);


}
void Geometry::WaveFrontObj::updateNormals(){

    for(auto v:_normals){
        v->zero();
    }

    for(int i=_normals.size();i<_beads->size();i++){
        _normals.append(QSharedPointer<Physics::VecD3d>(new Physics::VecD3d()));
    }
    for(int i=0;i<_tris->size();i++){
        auto tri=_tris->at(i);
        auto norm=tri->_norm;
        _normals[tri->_v[0]->ID]->add(norm);
        _normals[tri->_v[1]->ID]->add(norm);
        _normals[tri->_v[2]->ID]->add(norm);
    }
    for(auto v:_normals){
        v->nomilize();
    }
}
void Geometry::WaveFrontObj::updateIndecies(){

    auto bufferBytes=_indBuff->data();
    int size=(_tris->size()) * (3 ) * sizeof(ushort);
    if(bufferBytes.size()!=size){
        //bufferBytes=QByteArray();
        bufferBytes.resize(size);
    }
    ushort *rawIndexArray = reinterpret_cast<ushort *>(bufferBytes.data());
    for(int index=0;index<_tris->size();index++){
        auto tri=_tris->at(index);
        rawIndexArray[3*index]=tri->_v[0]->ID;
        rawIndexArray[3*index+1]=tri->_v[1]->ID;
        rawIndexArray[3*index+2]=tri->_v[2]->ID;
    }
    _indBuff->setData(bufferBytes);

}
void Geometry::WaveFrontObj::buildMesh(Qt3DRender::QGeometry *customGeometry){
    updateIndecies();
    updateVertecies();
    int numVer=_beads->size();
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
    _mesh->setVertexCount(_tris->size()*3);
}
Qt3DRender::QGeometryRenderer * Geometry::WaveFrontObj::mesh(){
    return _mesh;
}
