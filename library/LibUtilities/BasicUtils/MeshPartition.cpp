///////////////////////////////////////////////////////////////////////////////
//
//  File: MeshPartition.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <LibUtilities/BasicUtils/MeshPartition.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include <tinyxml.h>

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/CompressData.h>
#include <LibUtilities/BasicUtils/FieldIO.h>

#include <LibUtilities/Foundations/Foundations.hpp>


#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/detail/edge.hpp>
#include <boost/format.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        MeshPartitionFactory& GetMeshPartitionFactory()
        {
            typedef Loki::SingletonHolder<MeshPartitionFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy,
                Loki::SingleThreaded> Type;
            return Type::Instance();
        }

        MeshPartition::MeshPartition(const LibUtilities::SessionReaderSharedPtr& pSession) :
                m_isCompressed(false),
                m_numFields(0),
                m_fieldNameToId(),
                m_comm(pSession->GetComm()),
                m_weightingRequired(false),
                m_weightBnd(false),
                m_weightDofs(false)
        {
            ReadConditions(pSession);
            ReadGeometry(pSession);
            ReadExpansions(pSession);
        }

        MeshPartition::~MeshPartition()
        {

        }

        void MeshPartition::PartitionMesh(int nParts, bool shared,
                                          bool overlapping)
        {
            ASSERTL0(m_meshElements.size() >= nParts,
                     "Too few elements for this many processes.");
            m_shared = shared;

            if (m_weightingRequired)
            {
                WeightElements();
            }
            CreateGraph(m_mesh);
            PartitionGraph(m_mesh, nParts, m_localPartition, overlapping);
        }

        void MeshPartition::WriteLocalPartition(LibUtilities::SessionReaderSharedPtr& pSession)
        {
            TiXmlDocument vNew;
            TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
            vNew.LinkEndChild(decl);

            TiXmlElement* vElmtNektar;
            vElmtNektar = new TiXmlElement("NEKTAR");

            int rank = m_comm->GetRowComm()->GetRank();
            OutputPartition(pSession, m_localPartition[rank], vElmtNektar);

            vNew.LinkEndChild(vElmtNektar);

            std::string  dirname = pSession->GetSessionName() + "_xml"; 
            fs::path    pdirname(dirname);
            
            boost::format pad("P%1$07d.xml");
            pad % rank;
            fs::path    pFilename(pad.str());
            
            if(!fs::is_directory(dirname))
            {
                fs::create_directory(dirname);
            }
            
            fs::path fullpath = pdirname / pFilename; 
            vNew.SaveFile(PortablePath(fullpath));
        }

        void MeshPartition::WriteAllPartitions(LibUtilities::SessionReaderSharedPtr& pSession)
        {
            for (int i = 0; i < m_localPartition.size(); ++i)
            {
                TiXmlDocument vNew;
                TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
                vNew.LinkEndChild(decl);

                TiXmlElement* vElmtNektar;
                vElmtNektar = new TiXmlElement("NEKTAR");

                OutputPartition(pSession, m_localPartition[i], vElmtNektar);

                vNew.LinkEndChild(vElmtNektar);

                std::string  dirname = pSession->GetSessionName() + "_xml"; 
                fs::path    pdirname(dirname);
                
                boost::format pad("P%1$07d.xml");
                pad % i;
                fs::path    pFilename(pad.str());
                
                fs::path fullpath = pdirname / pFilename; 
                
                if(!fs::is_directory(dirname))
                {
                    fs::create_directory(dirname);
                }

                vNew.SaveFile(PortablePath(fullpath));
            }
        }

        void MeshPartition::GetCompositeOrdering(CompositeOrdering &composites)
        {
            std::map<int, MeshEntity>::iterator it;
            for (it  = m_meshComposites.begin();
                 it != m_meshComposites.end(); ++it)
            {
                composites[it->first] = it->second.list;
            }
        }

        void MeshPartition::GetBndRegionOrdering(BndRegionOrdering &bndRegs)
        {
            bndRegs = m_bndRegOrder;
        }


        void MeshPartition::ReadExpansions(const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            // Find the Expansions tag
            TiXmlElement *expansionTypes = pSession->GetElement("Nektar/Expansions");

            // Find the Expansion type
            TiXmlElement *expansion = expansionTypes->FirstChildElement();
            std::string   expType   = expansion->Value();

            /// Expansiontypes will contain plenty of data,
            /// where relevant at this stage are composite
            /// ID(s) that this expansion type describes,
            /// nummodes and a list of fields that this
            /// expansion relates to. If this does not exist
            /// the variable is only set to "DefaultVar".

            if(expType == "E")
            {
                while (expansion)
                {
                    std::vector<unsigned int> composite;
                    std::vector<unsigned int> nummodes;
                    std::vector<std::string>  fieldName;

                    const char *nModesStr = expansion->Attribute("NUMMODES");
                    ASSERTL0(nModesStr,"NUMMODES was not defined in EXPANSION section of input");
                    std::string numModesStr = nModesStr;
                    bool valid = ParseUtils::GenerateOrderedVector(numModesStr.c_str(), nummodes);
                    ASSERTL0(valid, "Unable to correctly parse the number of modes.");

                    if (nummodes.size() == 1)
                    {
                        for (int i = 1; i < m_dim; i++)
                        {
                            nummodes.push_back( nummodes[0] );
                        }
                    }
                    ASSERTL0(nummodes.size() == m_dim,"Number of modes should match mesh dimension");


                    const char *fStr = expansion->Attribute("FIELDS");
                    if(fStr)
                    {
                        std::string fieldStr = fStr;
                        bool  valid = ParseUtils::GenerateOrderedStringVector(fieldStr.c_str(),fieldName);
                        ASSERTL0(valid,"Unable to correctly parse the field string in ExpansionTypes.");

                        for (int i = 0; i < fieldName.size(); ++i)
                        {
                            if (m_fieldNameToId.count(fieldName[i]) == 0)
                            {
                                int k = m_fieldNameToId.size();
                                m_fieldNameToId[ fieldName[i] ] = k;
                                m_numFields++;
                            }
                        }
                    }
                    else
                    {
                        fieldName.push_back("DefaultVar");
                        int k = m_fieldNameToId.size();

                        if (m_fieldNameToId.count("DefaultVar") == 0)
                        {
                            ASSERTL0(k == 0,
                                     "Omitting field variables and explicitly listing " \
                                     "them in different ExpansionTypes is wrong practise");

                            m_fieldNameToId[ "DefaultVar" ] = k;
                            m_numFields++;
                        }
                    }

                    std::string compositeStr = expansion->Attribute("COMPOSITE");
                    ASSERTL0(compositeStr.length() > 3, "COMPOSITE must be specified in expansion definition");
                    int beg = compositeStr.find_first_of("[");
                    int end = compositeStr.find_first_of("]");
                    std::string compositeListStr = compositeStr.substr(beg+1,end-beg-1);
                    bool parseGood = ParseUtils::GenerateSeqVector(compositeListStr.c_str(), composite);
                    ASSERTL0(parseGood && !composite.empty(),
                        (std::string("Unable to read composite index range: ") + compositeListStr).c_str());


                    // construct mapping (elmt id, field name) -> nummodes
                    for (int i = 0; i < composite.size(); ++i)
                    {
                        for (int j = 0; j < fieldName.size(); j++)
                        {
                            for (unsigned int k = 0; k < m_meshComposites[composite[i]].list.size(); ++k)
                            {
                                int elid = m_meshComposites[composite[i]].list[k];
                                m_expansions[elid][fieldName[j]] = nummodes;
                                m_shape[elid] = m_meshComposites[composite[i]].type;
                            }
                        }
                    }

                    expansion = expansion->NextSiblingElement("E");
                }
            }
            else if(expType == "F")
            {
                    ASSERTL0(expansion->Attribute("FILE"),
                                "Attribute FILE expected for type F expansion");
                    std::string filenameStr = expansion->Attribute("FILE");
                    ASSERTL0(!filenameStr.empty(),
                                "A filename must be specified for the FILE "
                                "attribute of expansion");

                    // Create fieldIO object to load file
                    //    need a serial communicator to avoid problems with
                    //    shared file system
                    CommSharedPtr comm=
                            GetCommFactory().CreateInstance("Serial", 0, 0);
                    std::string iofmt = FieldIO::GetFileType(
                                filenameStr, comm);
                    FieldIOSharedPtr f = GetFieldIOFactory().CreateInstance(
                                iofmt,
                                comm,
                                pSession->GetSharedFilesystem());
                    // Load field definitions from file
                    std::vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefs;
                    f->Import(filenameStr, fielddefs);

                    // Parse field definitions
                    for (int i = 0; i < fielddefs.size(); ++i)
                    {
                        // Name of fields
                        for (int j = 0; j < fielddefs[i]->m_fields.size(); ++j)
                        {
                            std::string fieldName = fielddefs[i]->m_fields[j];
                            if (m_fieldNameToId.count(fieldName) == 0)
                            {
                                int k = m_fieldNameToId.size();
                                m_fieldNameToId[ fieldName ] = k;
                                m_numFields++;
                            }
                        }
                        // Number of modes and shape for each element
                        int numHomoDir = fielddefs[i]->m_numHomogeneousDir;
                        int cnt = 0;
                        for (int j = 0; j < fielddefs[i]->m_elementIDs.size(); ++j)
                        {
                            int elid = fielddefs[i]->m_elementIDs[j];
                            std::vector<unsigned int> nummodes;
                            for (int k = 0; k < m_dim; k++)
                            {
                                nummodes.push_back(fielddefs[i]->m_numModes[cnt++]);
                            }
                            if (fielddefs[i]->m_uniOrder)
                            {
                                cnt = 0;
                            }
                            else
                            {
                                cnt += numHomoDir;
                            }
                            for (int k = 0; k < fielddefs[i]->m_fields.size(); k++)
                            {
                                std::string fieldName = fielddefs[i]->m_fields[k];
                                m_expansions[elid][fieldName] = nummodes;
                            }
                            switch (fielddefs[i]->m_shapeType)
                            {
                                case eSegment:
                                {
                                    m_shape[elid] = 'S';
                                    break;
                                }
                                case eTriangle:
                                {
                                    m_shape[elid] = 'T';
                                    break;
                                }
                                case eQuadrilateral:
                                {
                                    m_shape[elid] = 'Q';
                                    break;
                                }
                                case eTetrahedron:
                                {
                                    m_shape[elid] = 'A';
                                    break;
                                }
                                case ePyramid:
                                {
                                    m_shape[elid] = 'R';
                                    break;
                                }
                                case ePrism:
                                {
                                    m_shape[elid] = 'P';
                                    break;
                                }
                                case eHexahedron:
                                {
                                    m_shape[elid] = 'H';
                                    break;
                                }
                                default:
                                    ASSERTL0 (false, "Shape not recognized.");
                                    break;
                            }
                        }
                    }
            }
            else
            {
                ASSERTL0(false,"Expansion type not defined or not supported at the moment");
            }
        }




        void MeshPartition::ReadGeometry(const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            TiXmlElement* x;
            TiXmlElement *vGeometry, *vSubElement;

            vGeometry = pSession->GetElement("Nektar/Geometry");
            m_dim = atoi(vGeometry->Attribute("DIM"));

            // Read mesh vertices
            vSubElement = pSession->GetElement("Nektar/Geometry/Vertex");

            // Retrieve any VERTEX attributes specifying mesh transforms
            std::string attr[] = {"XSCALE", "YSCALE", "ZSCALE",
                                  "XMOVE",  "YMOVE",  "ZMOVE" };
            for (int i = 0; i < 6; ++i)
            {
                const char *val =  vSubElement->Attribute(attr[i].c_str());
                if (val)
                {
                    m_vertexAttributes[attr[i]] = std::string(val);
                }
            }

            // check to see if compressed
            std::string IsCompressed;
            vSubElement->QueryStringAttribute("COMPRESSED",&IsCompressed); 

            if(IsCompressed.size()) 
            {
                ASSERTL0(boost::iequals(IsCompressed,
                                        CompressData::GetCompressString()),
                        "Compressed formats do not match. Expected :"
                             + CompressData::GetCompressString()
                             + "but got "+ std::string(IsCompressed));

                m_isCompressed = true;

                // Extract the vertex body
                TiXmlNode* vertexChild = vSubElement->FirstChild();
                ASSERTL0(vertexChild, "Unable to extract the data "
                         "from the compressed vertex tag.");

                std::string vertexStr;
                if (vertexChild->Type() == TiXmlNode::TINYXML_TEXT)
                {
                    vertexStr += vertexChild->ToText()->ValueStr();
                }

                std::vector<MeshVertex> vertData;
                CompressData::ZlibDecodeFromBase64Str(vertexStr,vertData);

                for(int i = 0; i < vertData.size(); ++i)
                {
                    m_meshVertices[vertData[i].id] = vertData[i];
                }
            }
            else
            {
                x = vSubElement->FirstChildElement();

                while(x)
                {
                    TiXmlAttribute* y = x->FirstAttribute();
                    ASSERTL0(y, "Failed to get attribute.");
                    MeshVertex v;
                    v.id = y->IntValue();
                    std::vector<std::string> vCoords;
                    std::string vCoordStr = x->FirstChild()->ToText()->Value();
                    boost::split(vCoords, vCoordStr, boost::is_any_of("\t "));
                    v.x = atof(vCoords[0].c_str());
                    v.y = atof(vCoords[1].c_str());
                    v.z = atof(vCoords[2].c_str());
                    m_meshVertices[v.id] = v;
                    x = x->NextSiblingElement();
                }
            }

            // Read mesh edges
            if (m_dim >= 2)
            {
                vSubElement = pSession->GetElement("Nektar/Geometry/Edge");
                ASSERTL0(vSubElement, "Cannot read edges");

                // check to see if compressed
                std::string IsCompressed;
                vSubElement->QueryStringAttribute("COMPRESSED",&IsCompressed); 

                if(IsCompressed.size()) 
                {
                    ASSERTL0(boost::iequals(IsCompressed,
                                        CompressData::GetCompressString()),
                            "Compressed formats do not match. Expected :"
                            + CompressData::GetCompressString()
                            + " but got "
                            + boost::lexical_cast<std::string>(IsCompressed));

                    m_isCompressed = true;

                    // Extract the edge body
                    TiXmlNode* edgeChild = vSubElement->FirstChild();
                    ASSERTL0(edgeChild,
                             "Unable to extract the data from the compressed "
                             "edge tag.");

                    std::string edgeStr;

                    if (edgeChild->Type() == TiXmlNode::TINYXML_TEXT)
                    {
                        edgeStr += edgeChild->ToText()->ValueStr();
                    }

                    std::vector<MeshEdge> edgeData;
                    CompressData::ZlibDecodeFromBase64Str(edgeStr,edgeData);

                    for(int i = 0; i < edgeData.size(); ++i)
                    {
                        MeshEntity e;
                        e.id = edgeData[i].id;
                        e.list.push_back(edgeData[i].v0);
                        e.list.push_back(edgeData[i].v1);
                        m_meshEdges[e.id] = e;
                    }
                }
                else
                {
                    x = vSubElement->FirstChildElement();

                    while(x)
                    {
                        TiXmlAttribute* y = x->FirstAttribute();
                        ASSERTL0(y, "Failed to get attribute.");
                        MeshEntity e;
                        e.id = y->IntValue();
                        e.type = 'E';
                        std::vector<std::string> vVertices;
                        std::string vVerticesString = x->FirstChild()->ToText()->Value();
                        boost::split(vVertices, vVerticesString, boost::is_any_of("\t "));
                        e.list.push_back(atoi(vVertices[0].c_str()));
                        e.list.push_back(atoi(vVertices[1].c_str()));
                        m_meshEdges[e.id] = e;
                        x = x->NextSiblingElement();
                    }
                }
            }
            
            // Read mesh faces
            if (m_dim == 3)
            {
                vSubElement = pSession->GetElement("Nektar/Geometry/Face");
                ASSERTL0(vSubElement, "Cannot read faces.");
                x = vSubElement->FirstChildElement();

                while(x)
                {
                    // check to see if compressed
                    std::string IsCompressed;
                    x->QueryStringAttribute("COMPRESSED",&IsCompressed); 

                    if(IsCompressed.size()) 
                    {
                        ASSERTL0(boost::iequals(IsCompressed,
                                        CompressData::GetCompressString()),
                                 "Compressed formats do not match. Expected :"
                                 + CompressData::GetCompressString()
                                 + " but got "
                                 + boost::lexical_cast<std::string>(
                                                                IsCompressed));
                        m_isCompressed = true;

                        // Extract the edge body
                        TiXmlNode* faceChild = x->FirstChild();
                        ASSERTL0(faceChild, "Unable to extract the data from "
                                            "the compressed edge tag.");

                        std::string faceStr;
                        if (faceChild->Type() == TiXmlNode::TINYXML_TEXT)
                        {
                            faceStr += faceChild->ToText()->ValueStr();
                        }

                        // uncompress and fill in values.
                        const std::string val= x->Value();

                        if(boost::iequals(val,"T")) // triangle
                        {
                            std::vector<MeshTri> faceData;
                            CompressData::ZlibDecodeFromBase64Str(faceStr,
                                                                  faceData);

                            for(int i= 0; i < faceData.size(); ++i)
                            {
                                MeshEntity f;

                                f.id = faceData[i].id;
                                f.type = 'T';
                                for (int j = 0; j < 3; ++j)
                                {
                                    f.list.push_back(faceData[i].e[j]);
                                }
                                m_meshFaces[f.id] = f;
                            }
                        }
                        else if (boost::iequals(val,"Q"))
                        {
                            std::vector<MeshQuad> faceData;
                            CompressData::ZlibDecodeFromBase64Str(faceStr,
                                                                  faceData);

                            for(int i= 0; i < faceData.size(); ++i)
                            {
                                MeshEntity f;

                                f.id = faceData[i].id;
                                f.type = 'Q';
                                for (int j = 0; j < 4; ++j)
                                {
                                    f.list.push_back(faceData[i].e[j]);
                                }
                                m_meshFaces[f.id] = f;
                            }
                        }
                        else
                        {
                            ASSERTL0(false,"Unrecognised face tag");
                        }
                    }
                    else
                    {
                        TiXmlAttribute* y = x->FirstAttribute();
                        ASSERTL0(y, "Failed to get attribute.");
                        MeshEntity f;
                        f.id = y->IntValue();
                        f.type = x->Value()[0];
                        std::vector<std::string> vEdges;
                        std::string vEdgeStr = x->FirstChild()->ToText()->Value();
                        boost::split(vEdges, vEdgeStr, boost::is_any_of("\t "));
                        for (int i = 0; i < vEdges.size(); ++i)
                        {
                            f.list.push_back(atoi(vEdges[i].c_str()));
                        }
                        m_meshFaces[f.id] = f;
                    }
                    x = x->NextSiblingElement();
                }
            }

            // Read mesh elements
            vSubElement = pSession->GetElement("Nektar/Geometry/Element");
            ASSERTL0(vSubElement, "Cannot read elements.");
            x = vSubElement->FirstChildElement();
            while(x)
            {
                // check to see if compressed
                std::string IsCompressed;
                x->QueryStringAttribute("COMPRESSED",&IsCompressed); 

                if(IsCompressed.size()) 
                {
                    ASSERTL0(boost::iequals(IsCompressed,
                                        CompressData::GetCompressString()),
                             "Compressed formats do not match. Expected :"
                             + CompressData::GetCompressString()
                             + " but got "
                             + boost::lexical_cast<std::string>(IsCompressed));
                    m_isCompressed = true;

                    // Extract the body
                    TiXmlNode* child = x->FirstChild();
                    ASSERTL0(child, "Unable to extract the data from the "
                                    "compressed element tag.");

                    std::string elmtStr;
                    if (child->Type() == TiXmlNode::TINYXML_TEXT)
                    {
                        elmtStr += child->ToText()->ValueStr();
                    }

                    // uncompress and fill in values.
                    const std::string val= x->Value();

                    switch(val[0])
                    {
                        case 'S': // segment
                        {
                            std::vector<MeshEdge> data;
                            CompressData::ZlibDecodeFromBase64Str(elmtStr,data);

                            for(int i= 0; i < data.size(); ++i)
                            {
                                MeshEntity e;
                                e.id = data[i].id;
                                e.type = 'S';
                                e.list.push_back(data[i].v0);
                                e.list.push_back(data[i].v1);
                                m_meshElements[e.id] = e;
                            }
                        }
                        break;
                        case 'T': // triangle
                        {
                            std::vector<MeshTri> data;
                            CompressData::ZlibDecodeFromBase64Str(elmtStr,data);

                            for(int i= 0; i < data.size(); ++i)
                            {
                                MeshEntity f;
                                f.id = data[i].id;
                                f.type = 'T';
                                for (int j = 0; j < 3; ++j)
                                {
                                    f.list.push_back(data[i].e[j]);
                                }
                                m_meshElements[f.id] = f;
                            }
                        }
                        break;
                        case 'Q': // quadrilateral
                        {
                            std::vector<MeshQuad> data;
                            CompressData::ZlibDecodeFromBase64Str(elmtStr,data);

                            for(int i= 0; i < data.size(); ++i)
                            {
                                MeshEntity f;
                                f.id = data[i].id;
                                f.type = 'Q';
                                for (int j = 0; j < 4; ++j)
                                {
                                    f.list.push_back(data[i].e[j]);
                                }
                                m_meshElements[f.id] = f;
                            }
                        }
                        break;
                        case 'A': // tetrahedron
                        {
                            std::vector<MeshTet> data;
                            CompressData::ZlibDecodeFromBase64Str(elmtStr,data);

                            for(int i= 0; i < data.size(); ++i)
                            {
                                MeshEntity f;
                                f.id = data[i].id;
                                f.type = 'A';
                                for (int j = 0; j < 4; ++j)
                                {
                                    f.list.push_back(data[i].f[j]);
                                }
                                m_meshElements[f.id] = f;
                            }
                        }
                        break;
                        case 'P': // prism
                        {
                            std::vector<MeshPyr> data;
                            CompressData::ZlibDecodeFromBase64Str(elmtStr,data);

                            for(int i= 0; i < data.size(); ++i)
                            {
                                MeshEntity f;
                                f.id = data[i].id;
                                f.type = 'P';
                                for (int j = 0; j < 5; ++j)
                                {
                                    f.list.push_back(data[i].f[j]);
                                }
                                m_meshElements[f.id] = f;
                            }
                        }
                        break;
                        case 'R': // pyramid
                        {
                            std::vector<MeshPrism> data;
                            CompressData::ZlibDecodeFromBase64Str(elmtStr,data);

                            for(int i= 0; i < data.size(); ++i)
                            {
                                MeshEntity f;
                                f.id = data[i].id;
                                f.type = 'R';
                                for (int j = 0; j < 5; ++j)
                                {
                                    f.list.push_back(data[i].f[j]);
                                }
                                m_meshElements[f.id] = f;
                            }
                        }
                        break;
                        case 'H': // hexahedron
                        {
                            std::vector<MeshHex> data;
                            CompressData::ZlibDecodeFromBase64Str(elmtStr,data);

                            for(int i= 0; i < data.size(); ++i)
                            {
                                MeshEntity f;
                                f.id = data[i].id;
                                f.type = 'H';
                                for (int j = 0; j < 6; ++j)
                                {
                                    f.list.push_back(data[i].f[j]);
                                }
                                m_meshElements[f.id] = f;
                            }
                        }
                        break;
                        default:
                            ASSERTL0(false,"Unrecognised element tag");
                    }
                }
                else
                {
                    TiXmlAttribute* y = x->FirstAttribute();
                    ASSERTL0(y, "Failed to get attribute.");
                    MeshEntity e;
                    e.id = y->IntValue();
                    std::vector<std::string> vItems;
                    std::string vItemStr = x->FirstChild()->ToText()->Value();
                    boost::split(vItems, vItemStr, boost::is_any_of("\t "));
                    for (int i = 0; i < vItems.size(); ++i)
                    {
                        e.list.push_back(atoi(vItems[i].c_str()));
                    }
                    e.type = x->Value()[0];
                    m_meshElements[e.id] = e;
                }
                x = x->NextSiblingElement();
            }

            // Read mesh curves
            if (pSession->DefinesElement("Nektar/Geometry/Curved"))
            {
                vSubElement = pSession->GetElement("Nektar/Geometry/Curved");

                // check to see if compressed
                std::string IsCompressed;
                vSubElement->QueryStringAttribute("COMPRESSED",&IsCompressed); 

                x = vSubElement->FirstChildElement();
                while(x)
                {
                    if(IsCompressed.size()) 
                    {
                        ASSERTL0(boost::iequals(IsCompressed,
                                            CompressData::GetCompressString()),
                                "Compressed formats do not match. Expected :"
                                + CompressData::GetCompressString()
                                + " but got "
                                + boost::lexical_cast<std::string>(
                                                            IsCompressed));

                        m_isCompressed = true;

                        const char *entitytype = x->Value();
                        // The compressed curved information is stored
                        // in two parts: MeshCurvedInfo and
                        // MeshCurvedPts.  MeshCurvedPts is just a
                        // list of MeshVertex values of unique vertex
                        // values from which we can make edges and
                        // faces.
                        //
                        // Then there is a list of NekInt64 pieces of
                        // information which make a MeshCurvedInfo
                        // struct. This contains information such as
                        // the curve id, the entity id, the number of
                        // curved points and the offset of where these
                        // points are stored in the pts vector of
                        // MeshVertex values. Finally the point type
                        // is also stored but in NekInt64 format
                        // rather than an enum for binary stride
                        // compatibility.
                        if(boost::iequals(entitytype,"E")||
                           boost::iequals(entitytype,"F"))
                        {
                            // read in data
                            std::string elmtStr;
                            TiXmlNode* child = x->FirstChild();

                            if (child->Type() == TiXmlNode::TINYXML_TEXT)
                            {
                                elmtStr = child->ToText()->ValueStr();
                            }

                            std::vector<MeshCurvedInfo> cinfo;
                            CompressData::ZlibDecodeFromBase64Str(elmtStr,
                                                                  cinfo);

                            // unpack list of curved edge or faces
                            for(int i = 0; i < cinfo.size(); ++i)
                            {
                                MeshCurved c;
                                c.id         = cinfo[i].id;
                                c.entitytype = entitytype[0];
                                c.entityid   = cinfo[i].entityid;
                                c.npoints    = cinfo[i].npoints;
                                c.type       = kPointsTypeStr[cinfo[i].ptype];
                                c.ptid       = cinfo[i].ptid;
                                c.ptoffset   = cinfo[i].ptoffset;
                                m_meshCurved[std::make_pair(c.entitytype,
                                                            c.id)] = c;
                            }
                        }
                        else if(boost::iequals(entitytype,"DATAPOINTS"))
                        {
                            MeshCurvedPts cpts;
                            NekInt id;

                            ASSERTL0(x->Attribute("ID", &id),
                                     "Failed to get ID from PTS section");
                            cpts.id = id;

                            // read in data
                            std::string elmtStr;

                            TiXmlElement* DataIdx =
                                x->FirstChildElement("INDEX");
                            ASSERTL0(DataIdx,
                                     "Cannot read data index tag in compressed "
                                     "curved section");

                            TiXmlNode* child = DataIdx->FirstChild();
                            if (child->Type() == TiXmlNode::TINYXML_TEXT)
                            {
                                elmtStr = child->ToText()->ValueStr();
                            }

                            CompressData::ZlibDecodeFromBase64Str(elmtStr,
                                                                  cpts.index);

                            TiXmlElement* DataPts
                                    = x->FirstChildElement("POINTS");
                            ASSERTL0(DataPts,
                                     "Cannot read data pts tag in compressed "
                                     "curved section");

                            child = DataPts->FirstChild();
                            if (child->Type() == TiXmlNode::TINYXML_TEXT)
                            {
                                elmtStr = child->ToText()->ValueStr();
                            }

                            CompressData::ZlibDecodeFromBase64Str(elmtStr,
                                                                  cpts.pts);

                            m_meshCurvedPts[cpts.id] = cpts;
                        }
                        else
                        {
                            ASSERTL0(false, "Unknown tag in curved section");
                        }

                        x = x->NextSiblingElement();
                    }
                    else
                    {
                        MeshCurved c;
                        ASSERTL0(x->Attribute("ID", &c.id),
                                 "Failed to get attribute ID");
                        c.type = std::string(x->Attribute("TYPE"));
                        ASSERTL0(!c.type.empty(),
                                 "Failed to get attribute TYPE");
                        ASSERTL0(x->Attribute("NUMPOINTS", &c.npoints),
                                 "Failed to get attribute NUMPOINTS");
                        c.data = x->FirstChild()->ToText()->Value();
                        c.entitytype = x->Value()[0];

                        if (c.entitytype == "E")
                        {
                            ASSERTL0(x->Attribute("EDGEID", &c.entityid),
                                     "Failed to get attribute EDGEID");
                        }
                        else if (c.entitytype == "F")
                        {
                            ASSERTL0(x->Attribute("FACEID", &c.entityid),
                                     "Failed to get attribute FACEID");
                        }
                        else
                        {
                            ASSERTL0(false, "Unknown curve type.");
                        }

                        m_meshCurved[std::make_pair(c.entitytype, c.id)] = c;
                        x = x->NextSiblingElement();
                    }
                }
            }

            // Read composites
            vSubElement = pSession->GetElement("Nektar/Geometry/Composite");
            ASSERTL0(vSubElement, "Cannot read composites.");
            x = vSubElement->FirstChildElement();
            while(x)
            {
                TiXmlAttribute* y = x->FirstAttribute();
                ASSERTL0(y, "Failed to get attribute.");
                MeshEntity c;
                c.id = y->IntValue();
                std::string vSeqStr = x->FirstChild()->ToText()->Value();
                c.type = vSeqStr[0];
                std::string::size_type indxBeg = vSeqStr.find_first_of('[') + 1;
                std::string::size_type indxEnd = vSeqStr.find_last_of(']') - 1;
                vSeqStr = vSeqStr.substr(indxBeg, indxEnd - indxBeg + 1);

                std::vector<unsigned int> vSeq;
                ParseUtils::GenerateSeqVector(vSeqStr.c_str(), vSeq);

                for (int i = 0; i < vSeq.size(); ++i)
                {
                    c.list.push_back(vSeq[i]);
                }
                m_meshComposites[c.id] = c;
                x = x->NextSiblingElement();
            }

            // Read Domain
            vSubElement = pSession->GetElement("Nektar/Geometry/Domain");
            ASSERTL0(vSubElement, "Cannot read domain");
            std::string vSeqStr = vSubElement->FirstChild()->ToText()->Value();
            std::string::size_type indxBeg = vSeqStr.find_first_of('[') + 1;
            std::string::size_type indxEnd = vSeqStr.find_last_of(']') - 1;
            vSeqStr = vSeqStr.substr(indxBeg, indxEnd - indxBeg + 1);
            ParseUtils::GenerateSeqVector(vSeqStr.c_str(), m_domain);
        }

        void MeshPartition::PrintPartInfo(std::ostream &out)
        {
            int nElmt = boost::num_vertices(m_mesh);
            int nPart = m_localPartition.size();

            out << "# Partition information:" << std::endl;
            out << "# No. elements  : " << nElmt << std::endl;
            out << "# No. partitions: " << nPart << std::endl;
            out << "# ID  nElmt  nLocDof  nBndDof" << std::endl;

            BoostVertexIterator vertit, vertit_end;
            std::vector<int> partElmtCount(nPart, 0);
            std::vector<int> partLocCount (nPart, 0);
            std::vector<int> partBndCount (nPart, 0);

            std::map<int, int> elmtSizes;
            std::map<int, int> elmtBndSizes;

            for (std::map<int, NummodesPerField>::iterator expIt =
                    m_expansions.begin(); expIt != m_expansions.end(); ++expIt)
            {
                int elid = expIt->first;
                NummodesPerField npf = expIt->second;

                for (NummodesPerField::iterator it = npf.begin(); it != npf.end(); ++it)
                {
                    ASSERTL0(it->second.size() == m_dim,
                        " Number of directional" \
                        " modes in expansion spec for element id = " +
                        boost::lexical_cast<std::string>(elid) +
                        " and field " +
                        boost::lexical_cast<std::string>(it->first) +
                        " does not correspond to mesh dimension");

                    int na = it->second[0];
                    int nb = 0;
                    int nc = 0;
                    if (m_dim >= 2)
                    {
                        nb = it->second[1];
                    }
                    if (m_dim == 3)
                    {
                        nc = it->second[2];
                    }

                    elmtSizes[elid]    = CalculateElementWeight(
                        m_shape[elid], false, na, nb, nc);
                    elmtBndSizes[elid] = CalculateElementWeight(
                        m_shape[elid], true,  na, nb, nc);
                }
            }

            for (boost::tie(vertit, vertit_end) = boost::vertices(m_mesh);
                 vertit != vertit_end; ++vertit)
            {
                int partId = m_mesh[*vertit].partition;
                partElmtCount[partId]++;
                partLocCount [partId] += elmtSizes[m_mesh[*vertit].id];
                partBndCount [partId] += elmtBndSizes[m_mesh[*vertit].id];
            }

            for (int i = 0; i < nPart; ++i)
            {
                out << i << " " << partElmtCount[i] << " " << partLocCount[i] << " " << partBndCount[i] << std::endl;
            }
        }

        void MeshPartition::ReadConditions(const SessionReaderSharedPtr& pSession)
        {
            if (!pSession->DefinesElement("Nektar/Conditions/SolverInfo"))
            {
                // No SolverInfo = no change of default action to weight
                // mesh graph.
                return;
            }

            TiXmlElement* solverInfoElement = 
                    pSession->GetElement("Nektar/Conditions/SolverInfo");

            TiXmlElement* solverInfo = 
                    solverInfoElement->FirstChildElement("I");
            ASSERTL0(solverInfo, "Cannot read SolverInfo tags");

            while (solverInfo)
            {
                // read the property name
                ASSERTL0(solverInfo->Attribute("PROPERTY"),
                         "Missing PROPERTY attribute in solver info "
                         "section. ");
                std::string solverProperty = 
                    solverInfo->Attribute("PROPERTY");
                ASSERTL0(!solverProperty.empty(),
                         "Solver info properties must have a non-empty "
                         "name. ");
                // make sure that solver property is capitalised
                std::string solverPropertyUpper =
                    boost::to_upper_copy(solverProperty);


                // read the value
                ASSERTL0(solverInfo->Attribute("VALUE"),
                        "Missing VALUE attribute in solver info section. ");
                std::string solverValue    = solverInfo->Attribute("VALUE");
                ASSERTL0(!solverValue.empty(),
                         "Solver info properties must have a non-empty value");
                // make sure that property value is capitalised
                std::string propertyValueUpper =
                    boost::to_upper_copy(solverValue);

                if (solverPropertyUpper == "WEIGHTPARTITIONS") 
                {
                    if (propertyValueUpper == "DOF")
                    {
                        m_weightingRequired = true;
                        m_weightDofs        = true;
                    }
                    else if (propertyValueUpper == "BOUNDARY")
                    {
                        m_weightingRequired = true;
                        m_weightBnd        = true;
                    }
                    else if (propertyValueUpper == "BOTH")
                    {
                        m_weightingRequired = true;
                        m_weightDofs        = true;
                        m_weightBnd        = true;
                    }
                    return;
                }
                solverInfo = solverInfo->NextSiblingElement("I");
            }
        }


        /*
         * Calculate element weights based on
         *   - element type (Q,T,H,P,R,A)
         *   - nummodes in expansion which this element belongs to via composite.
         *
         * For each element we prepare two vertex weightings, one associated
         * with the number of matrix elements associated with it (to balance
         * matrix multiplication work) and another associated
         * with all work which scales linearly with the number of its 
         * coefficients: communication, vector updates etc.
         *
         * \todo Refactor this code to explicitly represent performance model
         * and flexibly generate graph vertex weights depending on perf data.
         */
        void MeshPartition::WeightElements()
        {
            std::vector<unsigned int> weight(m_numFields, 1);
            std::map<int, MeshEntity>::iterator eIt;
            for (eIt = m_meshElements.begin(); eIt != m_meshElements.end(); ++eIt)
            {
                m_vertWeights[eIt->first] = weight;
                m_vertBndWeights[eIt->first] = weight;
                m_edgeWeights[eIt->first] = weight;
            }

            for (std::map<int, NummodesPerField>::iterator expIt =
                    m_expansions.begin(); expIt != m_expansions.end(); ++expIt)
            {
                int elid = expIt->first;
                NummodesPerField npf = expIt->second;

                for (NummodesPerField::iterator it = npf.begin(); it != npf.end(); ++it)
                {
                    ASSERTL0(it->second.size() == m_dim,
                        " Number of directional" \
                        " modes in expansion spec for element id = " +
                        boost::lexical_cast<std::string>(elid) +
                        " and field " +
                        boost::lexical_cast<std::string>(it->first) +
                        " does not correspond to mesh dimension");

                    int na = it->second[0];
                    int nb = 0;
                    int nc = 0;
                    if (m_dim >= 2)
                    {
                        nb = it->second[1];
                    }
                    if (m_dim == 3)
                    {
                        nc = it->second[2];
                    }

                    m_vertWeights[elid][m_fieldNameToId[it->first]] =
                            CalculateElementWeight(m_shape[elid], false,
                                                   na, nb, nc);
                    m_vertBndWeights[elid][m_fieldNameToId[it->first]] =
                            CalculateElementWeight(m_shape[elid], true,
                                                   na, nb, nc);
                    m_edgeWeights[elid][m_fieldNameToId[it->first]] =
                            CalculateEdgeWeight(m_shape[elid],
                                                   na, nb, nc);
                }
            } // for i
        }

        void MeshPartition::CreateGraph(BoostSubGraph& pGraph)
        {
            // Maps edge/face to first mesh element id.
            // On locating second mesh element id, graph edge is created instead.
            std::map<int, int> vGraphEdges;
            std::map<int, MeshEntity>::iterator eIt;
            int vcnt = 0;

            for (eIt = m_meshElements.begin(); eIt != m_meshElements.end();
                 ++eIt, ++vcnt)
            {
                BoostVertex v = boost::add_vertex(pGraph);
                pGraph[v].id = eIt->first;
                pGraph[v].partition = 0;

                if (m_weightingRequired)
                {
                    pGraph[v].weight     = m_vertWeights[eIt->first];
                    pGraph[v].bndWeight  = m_vertBndWeights[eIt->first];
                    pGraph[v].edgeWeight = m_edgeWeights[eIt->first];
                }

                // Process element entries and add graph edges
                for (unsigned j = 0; j < eIt->second.list.size(); ++j)
                {
                    int eId = eIt->second.list[j];

                    // Look to see if we've examined this edge/face before
                    // If so, we've got both graph vertices so add edge
                    if (vGraphEdges.find(eId) != vGraphEdges.end())
                    {
                        BoostEdge e = boost::add_edge(vcnt, vGraphEdges[eId], pGraph).first;
                        pGraph[e].id = vcnt;
                    }
                    else
                    {
                        vGraphEdges[eId] = vcnt;
                    }
                }
            }
        }

        /**
         * @brief Partition the graph.
         *
         * This routine partitions the graph @p pGraph into @p nParts, producing
         * subgraphs that are populated in @p pLocalPartition. If the @p
         * overlapping option is set (which is used for post-processing
         * purposes), the resulting partitions are extended to cover
         * neighbouring elements by additional vertex on the dual graph, which
         * produces overlapping partitions (i.e. the intersection of two
         * connected partitions is non-empty).
         *
         * @param pGraph           Graph to be partitioned.
         * @param nParts           Number of partitions.
         * @param pLocalPartition  Vector of sub-graphs representing each
         *                         partition.
         * @param overlapping      True if resulting partitions should overlap.
         */
        void MeshPartition::PartitionGraph(BoostSubGraph& pGraph,
                                           int nParts,
                                           std::vector<BoostSubGraph>& pLocalPartition,
                                           bool overlapping)
        {
            int i;
            int nGraphVerts = boost::num_vertices(pGraph);
            int nGraphEdges = boost::num_edges(pGraph);

            int ncon = 1;
            if (m_weightDofs && m_weightBnd)
            {
                ncon = 2;
            }
            // Convert boost graph into CSR format
            BoostVertexIterator    vertit, vertit_end;
            BoostAdjacencyIterator adjvertit, adjvertit_end;
            Array<OneD, int> part(nGraphVerts,0);

            if (m_comm->GetRowComm()->TreatAsRankZero())
            {
                int acnt = 0;
                int vcnt = 0;
                int nWeight = ncon*nGraphVerts;
                Array<OneD, int> xadj(nGraphVerts+1,0);
                Array<OneD, int> adjncy(2*nGraphEdges);
                Array<OneD, int> adjwgt(2*nGraphEdges, 1);
                Array<OneD, int> vwgt(nWeight, 1);
                Array<OneD, int> vsize(nGraphVerts, 1);

                for ( boost::tie(vertit, vertit_end) = boost::vertices(pGraph);
                      vertit != vertit_end;
                      ++vertit)
                {
                    for ( boost::tie(adjvertit, adjvertit_end) = boost::adjacent_vertices(*vertit,pGraph);
                          adjvertit != adjvertit_end;
                          ++adjvertit)
                    {
                        adjncy[acnt++] = *adjvertit;
                        if (m_weightingRequired)
                        {
                            adjwgt[acnt-1] = pGraph[*vertit].edgeWeight[0];
                        }
                    }

                    xadj[++vcnt] = acnt;

                    if (m_weightingRequired)
                    {
                        int ccnt = 0;
                        if (m_weightDofs)
                        {
                            vwgt[ncon*(vcnt-1)+ccnt] = pGraph[*vertit].weight[0];
                            ccnt++;
                        }
                        if (m_weightBnd)
                        {
                            vwgt[ncon*(vcnt-1)+ccnt] = pGraph[*vertit].bndWeight[0];
                        }
                    }
                }

                // Call Metis and partition graph
                int vol = 0;

                try
                {
                    //////////////////////////////////////////////////////
                    // On a cartesian communicator do mesh partiotion just on the first column
                    // so there is no doubt the partitions are all the same in all the columns
                    if(m_comm->GetColumnComm()->GetRank() == 0)
                    {
                        // Attempt partitioning using METIS.
                        PartitionGraphImpl(nGraphVerts, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nParts, vol, part);

                        // Check METIS produced a valid partition and fix if not.
                        CheckPartitions(nParts, part);
                        if (!m_shared)
                        {
                            // distribute among columns
                            for (i = 1; i < m_comm->GetColumnComm()->GetSize(); ++i)
                            {
                                m_comm->GetColumnComm()->Send(i, part);
                            }
                        }
                    }
                    else
                    {
                        m_comm->GetColumnComm()->Recv(0, part);
                    }
                    if (!m_shared)
                    {
                        m_comm->GetColumnComm()->Block();

                        //////////////////////////////////
                        // distribute among rows
                        for (i = 1; i < m_comm->GetRowComm()->GetSize(); ++i)
                        {
                            m_comm->GetRowComm()->Send(i, part);
                        }
                    }
                }
                catch (...)
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Error in calling metis to partition graph.");
                }
            }
            else
            {
                m_comm->GetRowComm()->Recv(0, part);
            }

            // Create boost subgraph for this process's partitions
            int nCols = nParts;
            pLocalPartition.resize(nCols);
            for (i = 0; i < nCols; ++i)
            {
                pLocalPartition[i] = pGraph.create_subgraph();
            }

            // Populate subgraph
            i = 0;
            for ( boost::tie(vertit, vertit_end) = boost::vertices(pGraph);
                  vertit != vertit_end;
                  ++vertit, ++i)
            {
                pGraph[*vertit].partition = part[i];
                boost::add_vertex(i, pLocalPartition[part[i]]);
            }

            // If the overlapping option is set (for post-processing purposes),
            // add vertices that correspond to the neighbouring elements.
            if (overlapping)
            {
                for ( boost::tie(vertit, vertit_end) = boost::vertices(pGraph);
                      vertit != vertit_end;
                      ++vertit)
                {
                    for (boost::tie(adjvertit, adjvertit_end) = boost::adjacent_vertices(*vertit,pGraph);
                         adjvertit != adjvertit_end; ++adjvertit)
                    {
                        if(part[*adjvertit] != part[*vertit])
                        {
                            boost::add_vertex(*adjvertit, pLocalPartition[part[*vertit]]);
                        }
                    }
                }
            }
        }


        void MeshPartition::CheckPartitions(int nParts, Array<OneD, int> &pPart)
        {
            unsigned int       i     = 0;
            unsigned int       cnt   = 0;
            bool               valid = true;

            // Check that every process has at least one element assigned
            for (i = 0; i < nParts; ++i)
            {
                cnt = std::count(pPart.begin(), pPart.end(), i);
                if (cnt == 0)
                {
                    valid = false;
                }
            }

            // If METIS produced an invalid partition, repartition naively.
            // Elements are assigned to processes in a round-robin fashion.
            // It is assumed that METIS failure only occurs when the number of
            // elements is approx. the number of processes, so this approach
            // should not be too inefficient communication-wise.
            if (!valid)
            {
                for (i = 0; i < pPart.num_elements(); ++i)
                {
                    pPart[i] = i % nParts;
                }
            }
        }


        void MeshPartition::OutputPartition(
                LibUtilities::SessionReaderSharedPtr& pSession,
                BoostSubGraph& pGraph,
                TiXmlElement* pNektar)
        {
            // Write Geometry data
            std::string vDim   = pSession->GetElement("Nektar/Geometry")->Attribute("DIM");
            std::string vSpace = pSession->GetElement("Nektar/Geometry")->Attribute("SPACE");
            std::string vPart  = boost::lexical_cast<std::string>(pGraph[*boost::vertices(pGraph).first].partition);
            TiXmlElement* vElmtGeometry = new TiXmlElement("GEOMETRY");
            vElmtGeometry->SetAttribute("DIM", vDim);
            vElmtGeometry->SetAttribute("SPACE", vSpace);
            vElmtGeometry->SetAttribute("PARTITION", vPart);

            TiXmlElement *vVertex  = new TiXmlElement("VERTEX");
            TiXmlElement *vEdge    = new TiXmlElement("EDGE");
            TiXmlElement *vFace    = new TiXmlElement("FACE");
            TiXmlElement *vElement = new TiXmlElement("ELEMENT");
            TiXmlElement *vCurved  = new TiXmlElement("CURVED");
            TiXmlElement *vComposite = new TiXmlElement("COMPOSITE");
            TiXmlElement *vDomain  = new TiXmlElement("DOMAIN");

            TiXmlElement *x;
            TiXmlText    *y;

            BoostVertexIterator    vertit, vertit_end;
            int id;

            std::map<int, MeshEntity> vComposites;
            std::map<int, MeshEntity> vElements;
            std::map<int, MeshEntity> vEdges;
            std::map<int, MeshEntity> vFaces;
            std::map<int, MeshVertex> vVertices;
            std::map<int, MeshEntity>::iterator vIt;
            std::map<int, MeshVertex>::iterator vVertIt;
            std::map<std::string, std::string>::iterator vAttrIt;

            std::vector<unsigned int> idxList;

            // Populate lists of elements, edges and vertices required.
            for (boost::tie(vertit, vertit_end) = boost::vertices(pGraph);
                 vertit != vertit_end;
                 ++vertit)
            {
                id = pGraph[*vertit].id;
                vElements[id] = m_meshElements[pGraph[*vertit].id];
            }

            std::map<int, MeshEntity> * vNext = &vElements;
            switch (m_dim)
            {
                case 3:
                {
                    // Compile list of faces
                    for (vIt = vNext->begin(); vIt != vNext->end(); vIt++)
                    {
                        for (unsigned int j = 0; j < vIt->second.list.size(); ++j)
                        {
                            id = vIt->second.list[j];
                            vFaces[id] = m_meshFaces[id];
                        }
                    }
                    vNext = &vFaces;
                }
                case 2:
                {
                    // Compile list of edges
                    for (vIt = vNext->begin(); vIt != vNext->end(); vIt++)
                    {
                        for (unsigned int j = 0; j < vIt->second.list.size(); ++j)
                        {
                            id = vIt->second.list[j];
                            vEdges[id] = m_meshEdges[id];
                        }
                    }
                    vNext = &vEdges;
                }
                case 1:
                {
                    // Compile list of vertices
                    for (vIt = vNext->begin(); vIt != vNext->end(); vIt++)
                    {
                        for (unsigned int j = 0; j < vIt->second.list.size(); ++j)
                        {
                            id = vIt->second.list[j];
                            vVertices[id] = m_meshVertices[id];
                        }
                    }
                }
            }

            // Generate XML data for these mesh entities
            if(m_isCompressed)
            {
                std::vector<MeshVertex> vertInfo;
                for (vVertIt = vVertices.begin();
                     vVertIt != vVertices.end(); vVertIt++)
                {
                    MeshVertex v;
                    v.id = vVertIt->first;
                    v.x = vVertIt->second.x;
                    v.y = vVertIt->second.y;
                    v.z = vVertIt->second.z;
                    vertInfo.push_back(v);
                }
                std::string vertStr;
                CompressData::ZlibEncodeToBase64Str(vertInfo,vertStr);
                vVertex->SetAttribute("COMPRESSED",
                                      CompressData::GetCompressString());
                vVertex->SetAttribute("BITSIZE",
                                      CompressData::GetBitSizeStr());
                vVertex->LinkEndChild(new TiXmlText(vertStr));
            }
            else
            {
                for (vVertIt = vVertices.begin(); vVertIt != vVertices.end(); vVertIt++)
                {
                    x = new TiXmlElement("V");
                    x->SetAttribute("ID", vVertIt->first);
                    std::stringstream vCoords;
                    vCoords.precision(12);
                    vCoords << std::setw(15) << vVertIt->second.x << " "
                            << std::setw(15) << vVertIt->second.y << " "
                            << std::setw(15) << vVertIt->second.z << " ";
                    y = new TiXmlText(vCoords.str());
                    x->LinkEndChild(y);
                    vVertex->LinkEndChild(x);
                }
            }

            // Apply transformation attributes to VERTEX section
            for (vAttrIt  = m_vertexAttributes.begin();
                 vAttrIt != m_vertexAttributes.end();
                 ++ vAttrIt)
            {
                vVertex->SetAttribute(vAttrIt->first, vAttrIt->second);
            }

            if (m_dim >= 2)
            {
                if(m_isCompressed)
                {
                    std::vector<MeshEdge> edgeInfo;
                    for (vIt = vEdges.begin(); vIt != vEdges.end(); vIt++)
                    {
                        MeshEdge e;
                        e.id = vIt->first;
                        e.v0 = vIt->second.list[0];
                        e.v1 = vIt->second.list[1];
                        edgeInfo.push_back(e);
                    }
                    std::string edgeStr;
                    CompressData::ZlibEncodeToBase64Str(edgeInfo,edgeStr);
                    vEdge->SetAttribute("COMPRESSED",
                                        CompressData::GetCompressString());
                    vEdge->SetAttribute("BITSIZE",
                                        CompressData::GetBitSizeStr());
                    vEdge->LinkEndChild(new TiXmlText(edgeStr));
                }
                else
                {
                    for (vIt = vEdges.begin(); vIt != vEdges.end(); vIt++)
                    {
                        x = new TiXmlElement("E");
                        x->SetAttribute("ID", vIt->first);
                        std::stringstream vVertices;
                        vVertices << std::setw(10) << vIt->second.list[0]
                                  << std::setw(10) << vIt->second.list[1]
                                  << " ";
                        y = new TiXmlText(vVertices.str());
                        x->LinkEndChild(y);
                        vEdge->LinkEndChild(x);
                    }
                }
            }

            if (m_dim >= 3)
            {

                if(m_isCompressed)
                {
                    std::vector<MeshTri>  TriFaceInfo;
                    std::vector<MeshQuad> QuadFaceInfo;

                    for (vIt = vFaces.begin(); vIt != vFaces.end(); vIt++)
                    {
                        switch(vIt->second.list.size())
                        {
                            case 3:
                            {
                                MeshTri f;
                                f.id = vIt->first;
                                for(int i = 0; i < 3; ++i)
                                {
                                    f.e[i] = vIt->second.list[i];
                                }
                                TriFaceInfo.push_back(f);
                            }
                            break;
                            case 4:
                            {
                                MeshQuad f;
                                f.id = vIt->first;
                                for(int i = 0; i < 4; ++i)
                                {
                                    f.e[i] = vIt->second.list[i];
                                }
                                QuadFaceInfo.push_back(f);
                            }
                            break;
                            default:
                                ASSERTL0(false,"Unknown face type.");
                        }
                    }

                    if(TriFaceInfo.size())
                    {
                        std::string vType("T");
                        x = new TiXmlElement(vType);

                        std::string faceStr;
                        CompressData::ZlibEncodeToBase64Str(TriFaceInfo,
                                                            faceStr);
                        x->SetAttribute("COMPRESSED",
                                        CompressData::GetCompressString());
                        x->SetAttribute("BITSIZE",
                                        CompressData::GetBitSizeStr());
                        x->LinkEndChild(new TiXmlText(faceStr));
                        vFace->LinkEndChild(x);
                    }

                    if(QuadFaceInfo.size())
                    {
                        std::string vType("Q");
                        x = new TiXmlElement(vType);
                        std::string faceStr;
                        CompressData::ZlibEncodeToBase64Str(QuadFaceInfo,
                                                            faceStr);
                        x->SetAttribute("COMPRESSED",
                                        CompressData::GetCompressString());
                        x->SetAttribute("BITSIZE",
                                        CompressData::GetBitSizeStr());
                        x->LinkEndChild(new TiXmlText(faceStr));
                        vFace->LinkEndChild(x);
                    }
                }
                else
                {
                    for (vIt = vFaces.begin(); vIt != vFaces.end(); vIt++)
                    {
                        std::string vType("F");
                        vType[0] = vIt->second.type;
                        x = new TiXmlElement(vType);
                        x->SetAttribute("ID", vIt->first);
                        std::stringstream vListStr;
                        for (unsigned int i = 0; i < vIt->second.list.size(); ++i)
                        {
                            vListStr << std::setw(10) << vIt->second.list[i];
                        }
                        vListStr << " ";
                        y = new TiXmlText(vListStr.str());
                        x->LinkEndChild(y);
                        vFace->LinkEndChild(x);
                    }
                }
            }


            if(m_isCompressed)
            {
                std::vector<MeshEdge>  SegInfo;
                std::vector<MeshTri>   TriInfo;
                std::vector<MeshQuad>  QuadInfo;
                std::vector<MeshTet>   TetInfo;
                std::vector<MeshPyr>   PyrInfo;
                std::vector<MeshPrism> PrismInfo;
                std::vector<MeshHex>   HexInfo;

                //gather elemements in groups.
                for (vIt = vElements.begin(); vIt != vElements.end(); vIt++)
                {
                    switch(vIt->second.type)
                    {
                        case 'S':
                        {
                            MeshEdge e;
                            e.id = vIt->first;
                            e.v0 = vIt->second.list[0];
                            e.v1 = vIt->second.list[1];
                            SegInfo.push_back(e);
                        }
                        break;
                        case 'T':
                        {
                            MeshTri f;
                            f.id = vIt->first;
                            for(int i = 0; i < 3; ++i)
                            {
                                f.e[i] = vIt->second.list[i];
                            }
                            TriInfo.push_back(f);
                        }
                        break;
                        case 'Q':
                        {
                            MeshQuad f;
                            f.id = vIt->first;
                            for(int i = 0; i < 4; ++i)
                            {
                                f.e[i] = vIt->second.list[i];
                            }
                            QuadInfo.push_back(f);
                        }
                        break;
                        case 'A':
                        {
                            MeshTet vol;
                            vol.id = vIt->first;
                            for(int i = 0; i < 4; ++i)
                            {
                                vol.f[i] = vIt->second.list[i];
                            }
                            TetInfo.push_back(vol);
                        }
                        break;
                        case 'P':
                        {
                            MeshPyr vol;
                            vol.id = vIt->first;
                            for(int i = 0; i < 5; ++i)
                            {
                                vol.f[i] = vIt->second.list[i];
                            }
                            PyrInfo.push_back(vol);
                        }
                        break;
                        case 'R':
                        {
                            MeshPrism vol;
                            vol.id = vIt->first;
                            for(int i = 0; i < 5; ++i)
                            {
                                vol.f[i] = vIt->second.list[i];
                            }
                            PrismInfo.push_back(vol);
                        }
                        break;
                        case 'H':
                        {
                            MeshHex vol;
                            vol.id = vIt->first;
                            for(int i = 0; i < 6; ++i)
                            {
                                vol.f[i] = vIt->second.list[i];
                            }
                            HexInfo.push_back(vol);
                        }
                        break;
                        default:
                            ASSERTL0(false,"Unknown element type");
                    }
                }

                if(SegInfo.size())
                {
                    std::string vType("S");
                    x = new TiXmlElement(vType);

                    std::string segStr;
                    CompressData::ZlibEncodeToBase64Str(SegInfo,segStr);
                    x->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                    x->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());
                    x->LinkEndChild(new TiXmlText(segStr));
                    vElement->LinkEndChild(x);
                }

                if(TriInfo.size())
                {
                    std::string vType("T");
                    x = new TiXmlElement(vType);

                    std::string faceStr;
                    CompressData::ZlibEncodeToBase64Str(TriInfo,faceStr);
                    x->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                    x->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());
                    x->LinkEndChild(new TiXmlText(faceStr));
                    vElement->LinkEndChild(x);
                }

                if(QuadInfo.size())
                {
                    std::string vType("Q");
                    x = new TiXmlElement(vType);
                    std::string faceStr;
                    CompressData::ZlibEncodeToBase64Str(QuadInfo,faceStr);
                    x->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                    x->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());
                    x->LinkEndChild(new TiXmlText(faceStr));
                    vElement->LinkEndChild(x);
                }

                if(TetInfo.size())
                {
                    std::string vType("A");
                    x = new TiXmlElement(vType);
                    std::string volStr;
                    CompressData::ZlibEncodeToBase64Str(TetInfo,volStr);
                    x->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                    x->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());
                    x->LinkEndChild(new TiXmlText(volStr));
                    vElement->LinkEndChild(x);
                }

                if(PyrInfo.size())
                {
                    std::string vType("P");
                    x = new TiXmlElement(vType);
                    std::string volStr;
                    CompressData::ZlibEncodeToBase64Str(PyrInfo,volStr);
                    x->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                    x->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());
                    x->LinkEndChild(new TiXmlText(volStr));
                    vElement->LinkEndChild(x);
                }

                if(PrismInfo.size())
                {
                    std::string vType("R");
                    x = new TiXmlElement(vType);
                    std::string volStr;
                    CompressData::ZlibEncodeToBase64Str(PrismInfo,volStr);
                    x->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                    x->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());
                    x->LinkEndChild(new TiXmlText(volStr));
                    vElement->LinkEndChild(x);
                }

                if(HexInfo.size())
                {
                    std::string vType("H");
                    x = new TiXmlElement(vType);
                    std::string volStr;
                    CompressData::ZlibEncodeToBase64Str(HexInfo,volStr);
                    x->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                    x->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());
                    x->LinkEndChild(new TiXmlText(volStr));
                    vElement->LinkEndChild(x);
                }

            }
            else
            {
                for (vIt = vElements.begin(); vIt != vElements.end(); vIt++)
                {
                    std::string vType("T");
                    vType[0] = vIt->second.type;
                    x = new TiXmlElement(vType.c_str());
                    x->SetAttribute("ID", vIt->first);
                    std::stringstream vEdges;
                    for (unsigned i = 0; i < vIt->second.list.size(); ++i)
                    {
                        vEdges << std::setw(10) << vIt->second.list[i];
                    }
                    vEdges << " ";
                    y = new TiXmlText(vEdges.str());
                    x->LinkEndChild(y);
                    vElement->LinkEndChild(x);
                }
            }

            if (m_dim >= 2)
            {
                std::map<MeshCurvedKey, MeshCurved>::const_iterator vItCurve;

                if(m_isCompressed)
                {
                    std::vector<MeshCurvedInfo> edgeinfo;
                    std::vector<MeshCurvedInfo> faceinfo;
                    MeshCurvedPts  curvedpts;
                    curvedpts.id = 0; // assume all points are going in here
                    int ptoffset = 0;
                    int newidx   = 0;
                    std::map<int,int> idxmap;

                    for (vItCurve  = m_meshCurved.begin();
                         vItCurve != m_meshCurved.end();
                         ++vItCurve)
                    {
                        MeshCurved c = vItCurve->second;

                        bool IsEdge = boost::iequals(c.entitytype,"E");
                        bool IsFace = boost::iequals(c.entitytype,"F");

                        if((IsEdge&&vEdges.find(c.entityid) != vEdges.end())||
                           (IsFace&&vFaces.find(c.entityid) != vFaces.end()))
                        {
                            MeshCurvedInfo cinfo;
                            // add in
                            cinfo.id       = c.id;
                            cinfo.entityid = c.entityid;
                            cinfo.npoints  = c.npoints;
                            for(int i = 0; i < SIZE_PointsType; ++i)
                            {
                                if(c.type.compare(kPointsTypeStr[i]) == 0)
                                {
                                    cinfo.ptype = (PointsType) i;
                                    break;
                                }
                            }
                            cinfo.ptid   = 0; // set to just one point set
                            cinfo.ptoffset = ptoffset;
                            ptoffset += c.npoints;

                            if (IsEdge)
                            {
                                edgeinfo.push_back(cinfo);
                            }
                            else
                            {
                                faceinfo.push_back(cinfo);
                            }

                            // fill in points to list.
                            for(int i =0; i < c.npoints; ++i)
                            {
                                // get index from full list;
                                int idx = m_meshCurvedPts[c.ptid]
                                    .index[c.ptoffset+i];

                                // if index is not already in curved
                                // points add it or set index to location
                                if(idxmap.count(idx) == 0)
                                {
                                    idxmap[idx] = newidx;
                                    curvedpts.index.push_back(newidx);
                                    curvedpts.pts.push_back(
                                            m_meshCurvedPts[c.ptid].pts[idx]);
                                    newidx++;
                                }
                                else
                                {
                                    curvedpts.index.push_back(idxmap[idx]);
                                }
                            }
                        }
                    }

                    // add xml information
                    if(edgeinfo.size())
                    {
                        vCurved->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                        vCurved->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());

                        x = new TiXmlElement("E");
                        std::string dataStr;
                        CompressData::ZlibEncodeToBase64Str(edgeinfo,dataStr);
                        x->LinkEndChild(new TiXmlText(dataStr));
                        vCurved->LinkEndChild(x);
                    }

                    if(faceinfo.size())
                    {
                        vCurved->SetAttribute("COMPRESSED",
                                    CompressData::GetCompressString());
                        vCurved->SetAttribute("BITSIZE",
                                    CompressData::GetBitSizeStr());

                        x = new TiXmlElement("F");
                        std::string dataStr;
                        CompressData::ZlibEncodeToBase64Str(faceinfo,dataStr);
                        x->LinkEndChild(new TiXmlText(dataStr));
                        vCurved->LinkEndChild(x);
                    }

                    if(edgeinfo.size()||faceinfo.size())
                    {
                        x = new TiXmlElement("DATAPOINTS");
                        x->SetAttribute("ID", curvedpts.id);

                        TiXmlElement *subx = new TiXmlElement("INDEX");
                        std::string dataStr;
                        CompressData::ZlibEncodeToBase64Str(curvedpts.index,
                                                            dataStr);
                        subx->LinkEndChild(new TiXmlText(dataStr));
                        x->LinkEndChild(subx);

                        subx = new TiXmlElement("POINTS");
                        CompressData::ZlibEncodeToBase64Str(curvedpts.pts,
                                                            dataStr);
                        subx->LinkEndChild(new TiXmlText(dataStr));
                        x->LinkEndChild(subx);

                        vCurved->LinkEndChild(x);
                    }
                }
                else
                {
                    for (vItCurve  = m_meshCurved.begin();
                         vItCurve != m_meshCurved.end();
                         ++vItCurve)
                    {
                        MeshCurved c = vItCurve->second;

                        bool IsEdge = boost::iequals(c.entitytype,"E");
                        bool IsFace = boost::iequals(c.entitytype,"F");

                        if((IsEdge&&vEdges.find(c.entityid) != vEdges.end())||
                           (IsFace&&vFaces.find(c.entityid) != vFaces.end()))
                        {
                            x = new TiXmlElement(c.entitytype);
                            x->SetAttribute("ID", c.id);
                            if (IsEdge)
                            {
                                x->SetAttribute("EDGEID", c.entityid);
                            }
                            else
                            {
                                x->SetAttribute("FACEID", c.entityid);
                            }
                            x->SetAttribute("TYPE", c.type);
                            x->SetAttribute("NUMPOINTS", c.npoints);
                            y = new TiXmlText(c.data);
                            x->LinkEndChild(y);
                            vCurved->LinkEndChild(x);
                        }
                    }
                }
            }

            // Generate composites section comprising only those mesh entities
            // which belong to this partition.
            for (vIt = m_meshComposites.begin(); vIt != m_meshComposites.end(); ++vIt)
            {
                idxList.clear();

                for (unsigned int j = 0; j < vIt->second.list.size(); ++j)
                {
                    // Based on entity type, check if in this partition
                    switch (vIt->second.type)
                    {
                    case 'V':
                        if (vVertices.find(vIt->second.list[j]) == vVertices.end())
                        {
                            continue;
                        }
                        break;
                    case 'E':
                        if (vEdges.find(vIt->second.list[j]) == vEdges.end())
                        {
                            continue;
                        }
                        break;
                    case 'F':
                        if (vFaces.find(vIt->second.list[j]) == vFaces.end())
                        {
                            continue;
                        }
                        break;
                    default:
                        if (vElements.find(vIt->second.list[j]) == vElements.end())
                        {
                            continue;
                        }
                        break;
                    }

                    idxList.push_back(vIt->second.list[j]);
                }

                std::string vCompositeStr = ParseUtils::GenerateSeqString(idxList);

                if (vCompositeStr.length() > 0)
                {
                    vComposites[vIt->first] = vIt->second;
                    x = new TiXmlElement("C");
                    x->SetAttribute("ID", vIt->first);
                    vCompositeStr = "X[" + vCompositeStr + "]";
                    vCompositeStr[0] = vIt->second.type;
                    y = new TiXmlText(vCompositeStr.c_str());
                    x->LinkEndChild(y);
                    vComposite->LinkEndChild(x);
                }
            }

            idxList.clear();
            std::string vDomainListStr;
            for (unsigned int i = 0; i < m_domain.size(); ++i)
            {
                if (vComposites.find(m_domain[i]) != vComposites.end())
                {
                    idxList.push_back(m_domain[i]);
                }
            }
            vDomainListStr = "C[" + ParseUtils::GenerateSeqString(idxList) + "]";
            TiXmlText* vDomainList = new TiXmlText(vDomainListStr);
            vDomain->LinkEndChild(vDomainList);

            vElmtGeometry->LinkEndChild(vVertex);
            if (m_dim >= 2)
            {
                vElmtGeometry->LinkEndChild(vEdge);
            }
            if (m_dim >= 3)
            {
                vElmtGeometry->LinkEndChild(vFace);
            }
            vElmtGeometry->LinkEndChild(vElement);
            if (m_dim >= 2)
            {
                vElmtGeometry->LinkEndChild(vCurved);
            }
            vElmtGeometry->LinkEndChild(vComposite);
            vElmtGeometry->LinkEndChild(vDomain);

            pNektar->LinkEndChild(vElmtGeometry);

            if (pSession->DefinesElement("Nektar/Conditions"))
            {
                std::set<int> vBndRegionIdList;
                TiXmlElement* vConditions    = new TiXmlElement(*pSession->GetElement("Nektar/Conditions"));
                TiXmlElement* vBndRegions    = vConditions->FirstChildElement("BOUNDARYREGIONS");
                TiXmlElement* vBndConditions = vConditions->FirstChildElement("BOUNDARYCONDITIONS");
                TiXmlElement* vItem;

                if (vBndRegions)
                {
                    TiXmlElement* vNewBndRegions = new TiXmlElement("BOUNDARYREGIONS");
                    vItem = vBndRegions->FirstChildElement();
                    while (vItem)
                    {
                        std::string vSeqStr = vItem->FirstChild()->ToText()->Value();
                        std::string::size_type indxBeg = vSeqStr.find_first_of('[') + 1;
                        std::string::size_type indxEnd = vSeqStr.find_last_of(']') - 1;
                        vSeqStr = vSeqStr.substr(indxBeg, indxEnd - indxBeg + 1);
                        std::vector<unsigned int> vSeq;
                        ParseUtils::GenerateSeqVector(vSeqStr.c_str(), vSeq);

                        idxList.clear();

                        for (unsigned int i = 0; i < vSeq.size(); ++i)
                        {
                            if (vComposites.find(vSeq[i]) != vComposites.end())
                            {
                                idxList.push_back(vSeq[i]);
                            }
                        }
                        int p = atoi(vItem->Attribute("ID"));

                        std::string vListStr = ParseUtils::GenerateSeqString(idxList);

                        if (vListStr.length() == 0)
                        {
                            TiXmlElement* tmp = vItem;
                            vItem = vItem->NextSiblingElement();
                            vBndRegions->RemoveChild(tmp);
                        }
                        else
                        {
                            vListStr = "C[" + vListStr + "]";
                            TiXmlText* vList = new TiXmlText(vListStr);
                            TiXmlElement* vNewElement = new TiXmlElement("B");
                            vNewElement->SetAttribute("ID", p);
                            vNewElement->LinkEndChild(vList);
                            vNewBndRegions->LinkEndChild(vNewElement);
                            vBndRegionIdList.insert(p);
                            vItem = vItem->NextSiblingElement();
                        }

                        // Store original order of boundary region.
                        m_bndRegOrder[p] = vSeq;
                        
                    }
                    vConditions->ReplaceChild(vBndRegions, *vNewBndRegions);
                }

                if (vBndConditions)
                {
                    vItem = vBndConditions->FirstChildElement();
                    while (vItem)
                    {
                        std::set<int>::iterator x;
                        if ((x = vBndRegionIdList.find(atoi(vItem->Attribute("REF")))) != vBndRegionIdList.end())
                        {
                            vItem->SetAttribute("REF", *x);
                            vItem = vItem->NextSiblingElement();
                        }
                        else
                        {
                            TiXmlElement* tmp = vItem;
                            vItem = vItem->NextSiblingElement();
                            vBndConditions->RemoveChild(tmp);
                        }
                    }
                }
                pNektar->LinkEndChild(vConditions);
            }

            // Distribute other sections of the XML to each process as is.
            TiXmlElement* vSrc = pSession->GetElement("Nektar")
                                                    ->FirstChildElement();
            while (vSrc)
            {
                std::string vName = boost::to_upper_copy(vSrc->ValueStr());
                if (vName != "GEOMETRY" && vName != "CONDITIONS")
                {
                    pNektar->LinkEndChild(new TiXmlElement(*vSrc));
                }
                vSrc = vSrc->NextSiblingElement();
            }
        }

        void MeshPartition::GetElementIDs(const int procid, std::vector<unsigned int> &elmtid)
        {
            BoostVertexIterator    vertit, vertit_end;

            ASSERTL0(procid < m_localPartition.size(),"procid is less than the number of partitions");
            
            // Populate lists of elements, edges and vertices required.
            for ( boost::tie(vertit, vertit_end) = boost::vertices(m_localPartition[procid]);
                  vertit != vertit_end;
                  ++vertit)
            {
                elmtid.push_back(m_meshElements[m_localPartition[procid][*vertit].id].id);
            }
        }

        int MeshPartition::CalculateElementWeight(
            char elmtType,
            bool bndWeight,
            int  na,
            int  nb,
            int  nc)
        {
            int weight = 0;

            switch (elmtType)
            {
                case 'A':
                    weight = bndWeight ?
                        StdTetData  ::getNumberOfBndCoefficients(na, nb, nc) :
                        StdTetData  ::getNumberOfCoefficients   (na, nb, nc);
                    break;
                case 'R':
                    weight = bndWeight ?
                        StdPrismData::getNumberOfBndCoefficients(na, nb, nc) :
                        StdPrismData::getNumberOfCoefficients   (na, nb, nc);
                    break;
                case 'H':
                    weight = bndWeight ?
                        StdHexData  ::getNumberOfBndCoefficients(na, nb, nc) :
                        StdHexData  ::getNumberOfCoefficients   (na, nb, nc);
                    break;
                case 'P':
                    weight = bndWeight ?
                        StdPyrData  ::getNumberOfBndCoefficients(na, nb, nc) :
                        StdPyrData  ::getNumberOfCoefficients   (na, nb, nc);
                    break;
                case 'Q':
                    weight = bndWeight ?
                        StdQuadData ::getNumberOfBndCoefficients(na, nb) :
                        StdQuadData ::getNumberOfCoefficients   (na, nb);
                    break;
                case 'T':
                    weight = bndWeight ?
                        StdTriData  ::getNumberOfBndCoefficients(na, nb) :
                        StdTriData  ::getNumberOfCoefficients   (na, nb);
                    break;
                case 'S':
                    weight = bndWeight ?
                        StdSegData  ::getNumberOfBndCoefficients(na) :
                        StdSegData  ::getNumberOfCoefficients   (na);
                    break;
                case 'V':
                    weight = 1;
                    break;
                default:
                    break;
            }

            return weight;
        }

        /**
         *     Calculate the number of modes needed for communication when
         *        in partition boundary, to be used as weighting for edges.
         *     Since we do not know exactly which face this refers to, assume
         *        the max order and quad face (for prisms) as arbitrary choices
         */
        int MeshPartition::CalculateEdgeWeight(
            char elmtType,
            int  na,
            int  nb,
            int  nc)
        {
            int weight = 0;
            int n = std::max ( na, std::max(nb, nc));
            switch (elmtType)
            {
                case 'A':
                    weight =
                        StdTriData ::getNumberOfCoefficients   (n, n);
                    break;
                case 'R':
                    weight =
                        StdQuadData ::getNumberOfCoefficients   (n, n);
                    break;
                case 'H':
                    weight =
                        StdQuadData ::getNumberOfCoefficients   (n, n);
                    break;
                case 'P':
                    weight =
                        StdQuadData ::getNumberOfCoefficients   (n, n);
                    break;
                case 'Q':
                case 'T':
                    weight = n;
                    break;
                case 'S':
                    weight = 1;
                    break;
                default:
                    break;
            }

            return weight;
        }
    }
}
