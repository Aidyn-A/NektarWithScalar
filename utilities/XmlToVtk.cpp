//////////////////////////////////////////////////////////////////////////////
//
// File: XmlToVtk.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Output VTK file of XML mesh, optionally with Jacobian.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iomanip>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>

using namespace std;
using namespace Nektar;

Array<OneD, NekDouble> GetQ(
        SpatialDomains::DerivStorage &deriv,
        LocalRegions::ExpansionSharedPtr exp);

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        cerr << "Usage: XmlToVtk  meshfile" << endl;
        exit(1);
    }

    LibUtilities::SessionReader::RegisterCmdLineFlag(
        "jacobian", "j", "Output Jacobian as scalar field");
    LibUtilities::SessionReader::RegisterCmdLineFlag(
        "quality", "q", "Output distribution of scaled Jacobians");

    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);

    bool jac = vSession->DefinesCmdLineArgument("jacobian");
    bool quality = vSession->DefinesCmdLineArgument("quality");

    jac = quality ? true : jac;

    // Read in mesh from input file
    string meshfile(argv[argc-1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt =
        SpatialDomains::MeshGraph::Read(vSession); //meshfile);

    // Set up Expansion information
    SpatialDomains::ExpansionMap emap = graphShPt->GetExpansions();
    SpatialDomains::ExpansionMapIter it;

    for (it = emap.begin(); it != emap.end(); ++it)
    {
        for (int i = 0; i < it->second->m_basisKeyVector.size(); ++i)
        {
            LibUtilities::BasisKey  tmp1 = it->second->m_basisKeyVector[i];
            LibUtilities::PointsKey tmp2 = tmp1.GetPointsKey();
            it->second->m_basisKeyVector[i] = LibUtilities::BasisKey(
                tmp1.GetBasisType(), tmp1.GetNumModes(),
                LibUtilities::PointsKey(tmp1.GetNumModes(),
                                        LibUtilities::ePolyEvenlySpaced));
        }
    }

    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(1);

    switch(expdim)
    {
        case 1:
        {
            if(vSession->DefinesSolverInfo("HOMOGENEOUS"))
            {
                std::string HomoStr = vSession->GetSolverInfo("HOMOGENEOUS");
                MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;

                ASSERTL0(
                    HomoStr == "HOMOGENEOUS1D" || HomoStr == "Homogeneous1D" ||
                    HomoStr == "1D"            || HomoStr == "Homo1D",
                    "Only 3DH1D supported for XML output currently.");

                int nplanes;
                vSession->LoadParameter("HomModesZ", nplanes);

                // choose points to be at evenly spaced points at nplanes + 1
                // points
                const LibUtilities::PointsKey Pkey(
                    nplanes + 1, LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(
                    LibUtilities::eFourier, nplanes, Pkey);
                NekDouble lz = vSession->GetParameter("LZ");

                Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>
                    ::AllocateSharedPtr(
                        vSession, Bkey, lz, false, false, graphShPt);
                Exp[0] = Exp2DH1;
            }
            else
            {
                MultiRegions::ExpList1DSharedPtr Exp1D;
                Exp1D = MemoryManager<MultiRegions::ExpList1D>
                    ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] = Exp1D;
            }

            break;
        }
        case 2:
        {
            if(vSession->DefinesSolverInfo("HOMOGENEOUS"))
            {
                std::string HomoStr = vSession->GetSolverInfo("HOMOGENEOUS");
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                ASSERTL0(
                    HomoStr == "HOMOGENEOUS1D" || HomoStr == "Homogeneous1D" ||
                    HomoStr == "1D"            || HomoStr == "Homo1D",
                    "Only 3DH1D supported for XML output currently.");

                int nplanes;
                vSession->LoadParameter("HomModesZ", nplanes);

                // choose points to be at evenly spaced points at nplanes + 1
                // points
                const LibUtilities::PointsKey Pkey(
                    nplanes + 1, LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(
                    LibUtilities::eFourier, nplanes, Pkey);
                NekDouble lz = vSession->GetParameter("LZ");

                Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                    ::AllocateSharedPtr(
                        vSession, Bkey, lz, false, false, graphShPt);
                Exp[0] = Exp3DH1;
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                    ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] =  Exp2D;
            }
            break;
        }
        case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp3D;
            break;
        }
        default:
        {
            ASSERTL0(false,"Expansion dimension not recognised");
            break;
        }
    }

    // Write out VTK file.
    string   outname(strtok(argv[argc-1],"."));
    outname += ".vtu";
    ofstream outfile(outname.c_str());

    Exp[0]->WriteVtkHeader(outfile);

    if (jac)
    {
        // Find minimum Jacobian.
        Array<OneD, NekDouble> tmp;
        Array<OneD, NekDouble> x0 (Exp[0]->GetNpoints());
        Array<OneD, NekDouble> x1 (Exp[0]->GetNpoints());
        Array<OneD, NekDouble> x2 (Exp[0]->GetNpoints());
        Exp[0]->GetCoords(x0, x1, x2);

        // Write out field containing Jacobian.
        for(int i = 0; i < Exp[0]->GetExpSize(); ++i)
        {
            LocalRegions::ExpansionSharedPtr e = Exp[0]->GetExp(i);
            SpatialDomains::GeometrySharedPtr geom = e->GetGeom();
            StdRegions::StdExpansionSharedPtr chi  = geom->GetXmap();
            LibUtilities::PointsKeyVector     p    = chi->GetPointsKeys();
            SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

            // Grab Jacobian entries (on collapsed quadrature points)
            SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);
            const int npts  = chi->GetTotPoints();


            Array<OneD, const NekDouble> q = GetQ(deriv, e);
            Vmath::Vcopy(npts, q, 1,
                         tmp = Exp[0]->UpdatePhys()
                                 + Exp[0]->GetPhys_Offset(i), 1);



            Exp[0]->WriteVtkPieceHeader(outfile, i);
            Exp[0]->WriteVtkPieceData  (outfile, i, "Q");
            Exp[0]->WriteVtkPieceFooter(outfile, i);
        }

        unsigned int n
            = Vmath::Imin(Exp[0]->GetNpoints(), Exp[0]->GetPhys(), 1);
        cout << "- Minimum Jacobian: "
             << Vmath::Vmin(Exp[0]->GetNpoints(), Exp[0]->GetPhys(), 1)
             << " at coords (" << x0[n] << ", " << x1[n] << ", " << x2[n] << ")"
             << endl;

    }
    else
    {
        // For each field write header and footer, since there is no field data.
        for(int i = 0; i < Exp[0]->GetExpSize(); ++i)
        {
            Exp[0]->WriteVtkPieceHeader(outfile, i);
            Exp[0]->WriteVtkPieceFooter(outfile, i);
        }
    }

    Exp[0]->WriteVtkFooter(outfile);

    return 0;
}

Array<OneD, NekDouble> GetQ(
        SpatialDomains::DerivStorage &deriv,
        LocalRegions::ExpansionSharedPtr exp)
{
    StdRegions::StdExpansionSharedPtr chi = exp->GetGeom()->GetXmap();
    const int pts = deriv[0][0].num_elements();
    const int nq  = chi->GetTotPoints();
    const int expDim = chi->GetNumBases();

    Array<OneD, NekDouble> jac(pts, 0.0);
    Array<OneD, NekDouble> eta(nq);

    switch (expDim)
    {
        case 2:
        {
            Vmath::Vvtvvtm(pts, &deriv[0][0][0], 1, &deriv[1][1][0], 1,
                           &deriv[1][0][0], 1, &deriv[0][1][0], 1,
                           &jac[0],         1);
            break;
        }
        case 3:
        {
            Array<OneD, NekDouble> tmp(pts, 0.0);
            Vmath::Vvtvvtm(pts, &deriv[1][1][0], 1, &deriv[2][2][0], 1,
                           &deriv[2][1][0], 1, &deriv[1][2][0], 1,
                           &tmp[0],         1);
            Vmath::Vvtvp  (pts, &deriv[0][0][0], 1, &tmp[0],         1,
                           &jac[0],         1, &jac[0],         1);
            Vmath::Vvtvvtm(pts, &deriv[2][1][0], 1, &deriv[0][2][0], 1,
                           &deriv[0][1][0], 1, &deriv[2][2][0], 1,
                           &tmp[0],         1);
            Vmath::Vvtvp  (pts, &deriv[1][0][0], 1, &tmp[0],         1,
                           &jac[0],         1, &jac[0],         1);
            Vmath::Vvtvvtm(pts, &deriv[0][1][0], 1, &deriv[1][2][0], 1,
                           &deriv[1][1][0], 1, &deriv[0][2][0], 1,
                           &tmp[0],         1);
            Vmath::Vvtvp  (pts, &deriv[2][0][0], 1, &tmp[0],         1,
                           &jac[0],         1, &jac[0],         1);
            break;
        }
    }

    for (int k = 0; k < pts; ++k)
    {
        NekDouble frob = 0.0;

        for (int i = 0; i < deriv.num_elements(); ++i)
        {
            for (int j = 0; j < deriv[i].num_elements(); ++j)
            {
                frob += deriv[i][j][k]*deriv[i][j][k];
            }
        }

        NekDouble delta = 0.0; // fixme
        NekDouble sigma = 0.5*(jac[k] + sqrt(jac[k]*jac[k] + 4*delta*delta));

        eta[k] = expDim * pow(sigma,2.0/expDim) / frob;
        if(eta[k] > 1 || eta[k] < 0)
            cout << eta[k] << endl;
    }

    if (pts == 1)
    {
        Vmath::Fill(nq-1, eta[0], &eta[1], 1);
    }

    return eta;
}
