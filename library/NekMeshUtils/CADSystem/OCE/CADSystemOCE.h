////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_CADSYSTEM_OCE_CADSYSTEMOCE
#define NekMeshUtils_CADSYSTEM_OCE_CADSYSTEMOCE

#include <NekMeshUtils/CADSystem/CADSystem.h>
#include <NekMeshUtils/CADSystem/OCE/OpenCascade.h>

namespace Nektar
{
namespace NekMeshUtils
{

class CADSystemOCE : public CADSystem
{
public:

    static CADSystemSharedPtr create(std::string name)
    {
        return MemoryManager<CADSystemOCE>::AllocateSharedPtr(name);
    }

    static std::string key;

    /**
     * @brief Default constructor.
     */
    CADSystemOCE(std::string name) : CADSystem(name) {}
    ~CADSystemOCE(){};

    bool LoadCAD();

    Array<OneD, NekDouble> GetBoundingBox();

private:
    /// Function to add curve to CADSystem::m_verts.
    void AddVert(int i, TopoDS_Shape in);
    /// Function to add curve to CADSystem::m_curves.
    void AddCurve(int i, TopoDS_Shape in);
    /// Function to add surface to CADSystem::m_surfs.
    void AddSurf(int i, TopoDS_Shape in);

    TopoDS_Shape BuildNACA(std::string naca);
    /// OCC master object
    TopoDS_Shape shape;
    TopTools_IndexedMapOfShape mapOfVerts, mapOfEdges, mapOfFaces;
};


}
}

#endif
