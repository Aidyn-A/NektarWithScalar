////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessInterpPoints.cpp
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
//  Description: Interpolate  field to a series of specified points.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/geometry.hpp>
#include "ProcessInterpPoints.h"

#include <FieldUtils/Interpolator.h>

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessInterpPoints::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "interppoints"),
        ProcessInterpPoints::create,
        "Interpolates a set of points to another, requires fromfld and "
        "fromxml to be defined, a line, plane or block of points can be "
        "defined");

ProcessInterpPoints::ProcessInterpPoints(FieldSharedPtr f) : ProcessModule(f)
{

    m_config["fromxml"] = ConfigOption(
        false, "NotSet", "Xml file from which to interpolate field");

    ASSERTL0(m_config["fromxml"].as<string>().compare("NotSet") != 0,
             "Need to specify fromxml=file.xml");

    m_config["fromfld"] = ConfigOption(
        false, "NotSet", "Fld file from which to interpolate field");

    ASSERTL0(m_config["fromfld"].as<string>().compare("NotSet") != 0,
             "Need to specify fromfld=file.fld ");

    m_config["clamptolowervalue"] =
        ConfigOption(false, "-10000000", "Lower bound for interpolation value");
    m_config["clamptouppervalue"] =
        ConfigOption(false, "10000000", "Upper bound for interpolation value");
    m_config["defaultvalue"] =
        ConfigOption(false, "0", "Default value if point is outside domain");
    m_config["line"] =
        ConfigOption(false, "NotSet", "Specify a line of N points using "
                                      "line=N,x0,y0,z0,z1,y1,z1");
    m_config["plane"] = ConfigOption(
        false, "NotSet", "Specify a plane of N1 x N2 points using "
                         "plane=N1,N2,x0,y0,z0,z1,y1,z1,x2,y2,z2,x3,"
                         "y3,z3");
    m_config["box"] = ConfigOption(
        false, "NotSet", "Specify a rectangular box of N1 x N2 x N3 points "
                         "using a box of points limited by box="
                         "N1,N2,N3,xmin,xmax,ymin,ymax,zmin,zmax");

    m_config["cp"] =
        ConfigOption(false, "NotSet",
                     "Parameters p0 and q to determine pressure coefficients "
                     "(box only currently)");
}

ProcessInterpPoints::~ProcessInterpPoints()
{
}

void ProcessInterpPoints::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "ProcessInterpPoints: interpolating to points..." << endl;
        }
    }

    int rank   = m_f->m_comm->GetRank();
    int nprocs = m_f->m_comm->GetSize();

    // Check for command line point specification if no .pts file
    // specified
    if (m_f->m_fieldPts == LibUtilities::NullPtsField)
    {
        if (m_config["line"].as<string>().compare("NotSet") != 0)
        {
            string help = m_config["line"].as<string>();
            vector<NekDouble> values;
            ASSERTL0(ParseUtils::GenerateUnOrderedVector(
                         m_config["line"].as<string>().c_str(), values),
                     "Failed to interpret line string");

            ASSERTL0(values.size() > 2,
                     "line string should contain 2 Dim+1 values "
                     "N,x0,y0,z0,x1,y1,z1");

            double tmp;
            ASSERTL0(std::modf(values[0], &tmp) == 0.0, "N is not an integer");
            ASSERTL0(values[0] > 1, "N is not a valid number");

            int dim  = (values.size() - 1) / 2;
            int npts = values[0];
            Array<OneD, Array<OneD, NekDouble> > pts(dim);

            for (int i = 0; i < dim; ++i)
            {
                pts[i] = Array<OneD, NekDouble>(npts);
            }

            for (int i = 0; i < npts; ++i)
            {
                pts[0][i] =
                    values[1] +
                    i / ((NekDouble)(npts - 1)) * (values[dim + 1] - values[1]);
                if (dim > 1)
                {
                    pts[1][i] = values[2] +
                                i / ((NekDouble)(npts - 1)) *
                                    (values[dim + 2] - values[2]);

                    if (dim > 2)
                    {
                        pts[2][i] = values[3] +
                                    i / ((NekDouble)(npts - 1)) *
                                        (values[dim + 3] - values[3]);
                    }
                }
            }

            vector<int> ppe;
            ppe.push_back(npts);
            m_f->m_fieldPts =
                MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(dim,
                                                                         pts);
            m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsLine);
            m_f->m_fieldPts->SetPointsPerEdge(ppe);
        }
        else if (m_config["plane"].as<string>().compare("NotSet") != 0)
        {
            string help = m_config["plane"].as<string>();
            vector<NekDouble> values;
            ASSERTL0(ParseUtils::GenerateUnOrderedVector(
                         m_config["plane"].as<string>().c_str(), values),
                     "Failed to interpret plane string");

            ASSERTL0(values.size() > 9,
                     "plane string should contain 4 Dim+2 values "
                     "N1,N2,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3");

            double tmp;
            ASSERTL0(std::modf(values[0], &tmp) == 0.0, "N1 is not an integer");
            ASSERTL0(std::modf(values[1], &tmp) == 0.0, "N2 is not an integer");

            ASSERTL0(values[0] > 1, "N1 is not a valid number");
            ASSERTL0(values[1] > 1, "N2 is not a valid number");

            int dim = (values.size() - 2) / 4;

            int npts1 = values[0];
            int npts2 = values[1];

            Array<OneD, Array<OneD, NekDouble> > pts(dim);

            int totpts  = npts1 * npts2;
            int nlocpts = totpts / nprocs;

            if (rank < nprocs - 1)
            {
                for (int i = 0; i < dim; ++i)
                {
                    pts[i] = Array<OneD, NekDouble>(nlocpts);
                }

                int cnt    = 0;
                int cntloc = 0;

                for (int j = 0; j < npts2; ++j)
                {
                    for (int i = 0; i < npts1; ++i)
                    {

                        if ((cnt >= rank * nlocpts) &&
                            (cnt < (rank + 1) * nlocpts))
                        {
                            pts[0][cntloc] =
                                (values[2] +
                                 i / ((NekDouble)(npts1 - 1)) *
                                     (values[dim + 2] - values[2])) *
                                    (1.0 - j / ((NekDouble)(npts2 - 1))) +
                                (values[3 * dim + 2] +
                                 i / ((NekDouble)(npts1 - 1)) *
                                     (values[2 * dim + 2] -
                                      values[3 * dim + 2])) *
                                    (j / ((NekDouble)(npts2 - 1)));

                            pts[1][cntloc] =
                                (values[3] +
                                 i / ((NekDouble)(npts1 - 1)) *
                                     (values[dim + 3] - values[3])) *
                                    (1.0 - j / ((NekDouble)(npts2 - 1))) +
                                (values[3 * dim + 3] +
                                 i / ((NekDouble)(npts1 - 1)) *
                                     (values[2 * dim + 3] -
                                      values[3 * dim + 3])) *
                                    (j / ((NekDouble)(npts2 - 1)));

                            if (dim > 2)
                            {
                                pts[2][cntloc] =
                                    (values[4] +
                                     i / ((NekDouble)(npts1 - 1)) *
                                         (values[dim + 4] - values[4])) *
                                        (1.0 - j / ((NekDouble)(npts2 - 1))) +
                                    (values[3 * dim + 4] +
                                     i / ((NekDouble)(npts1 - 1)) *
                                         (values[2 * dim + 4] -
                                          values[3 * dim + 4])) *
                                        (j / ((NekDouble)(npts2 - 1)));
                            }
                            cntloc++;
                        }
                        cnt++;
                    }
                }
            }
            else
            {
                totpts = totpts - rank * nlocpts;

                for (int i = 0; i < dim; ++i)
                {
                    pts[i] = Array<OneD, NekDouble>(totpts);
                }

                int cnt    = 0;
                int cntloc = 0;

                for (int j = 0; j < npts2; ++j)
                {
                    for (int i = 0; i < npts1; ++i)
                    {

                        if (cnt >= rank * nlocpts)
                        {
                            pts[0][cntloc] =
                                (values[2] +
                                 i / ((NekDouble)(npts1 - 1)) *
                                     (values[dim + 2] - values[2])) *
                                    (1.0 - j / ((NekDouble)(npts2 - 1))) +
                                (values[3 * dim + 2] +
                                 i / ((NekDouble)(npts1 - 1)) *
                                     (values[2 * dim + 2] -
                                      values[3 * dim + 2])) *
                                    (j / ((NekDouble)(npts2 - 1)));

                            pts[1][cntloc] =
                                (values[3] +
                                 i / ((NekDouble)(npts1 - 1)) *
                                     (values[dim + 3] - values[3])) *
                                    (1.0 - j / ((NekDouble)(npts2 - 1))) +
                                (values[3 * dim + 3] +
                                 i / ((NekDouble)(npts1 - 1)) *
                                     (values[2 * dim + 3] -
                                      values[3 * dim + 3])) *
                                    (j / ((NekDouble)(npts2 - 1)));

                            if (dim > 2)
                            {
                                pts[2][cntloc] =
                                    (values[4] +
                                     i / ((NekDouble)(npts1 - 1)) *
                                         (values[dim + 4] - values[4])) *
                                        (1.0 - j / ((NekDouble)(npts2 - 1))) +
                                    (values[3 * dim + 4] +
                                     i / ((NekDouble)(npts1 - 1)) *
                                         (values[2 * dim + 4] -
                                          values[3 * dim + 4])) *
                                        (j / ((NekDouble)(npts2 - 1)));
                            }
                            cntloc++;
                        }
                        cnt++;
                    }
                }
            }

            vector<int> ppe;
            ppe.push_back(npts1);
            ppe.push_back(npts2);
            m_f->m_fieldPts =
                MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(dim,
                                                                         pts);
            m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsPlane);
            m_f->m_fieldPts->SetPointsPerEdge(ppe);
        }
        else if (m_config["box"].as<string>().compare("NotSet") != 0)
        {
            string help = m_config["box"].as<string>();
            vector<NekDouble> values;
            ASSERTL0(ParseUtils::GenerateUnOrderedVector(
                         m_config["box"].as<string>().c_str(), values),
                     "Failed to interpret box string");

            ASSERTL0(values.size() == 9,
                     "box string should contain 9 values "
                     "N1,N2,N3,xmin,xmax,ymin,ymax,zmin,zmax");

            int dim = 3;

            int npts1 = values[0];
            int npts2 = values[1];
            int npts3 = values[2];

            Array<OneD, Array<OneD, NekDouble> > pts(dim);

            int totpts  = npts1 * npts2 * npts3;
            int nlocpts = totpts / nprocs;

            if (rank < nprocs - 1) // for rank 0 to nproc-1
            {
                totpts = nlocpts;

                for (int i = 0; i < dim; ++i)
                {
                    pts[i] = Array<OneD, NekDouble>(totpts);
                }

                int cnt    = 0;
                int cntloc = 0;

                for (int k = 0; k < npts3; ++k)
                {
                    for (int j = 0; j < npts2; ++j)
                    {
                        for (int i = 0; i < npts1; ++i)
                        {
                            if ((cnt >= rank * nlocpts) &&
                                (cnt < (rank + 1) * nlocpts))
                            {
                                pts[0][cntloc] = values[3] +
                                                 i / ((NekDouble)(npts1 - 1)) *
                                                     (values[4] - values[3]);
                                pts[1][cntloc] = values[5] +
                                                 j / ((NekDouble)(npts2 - 1)) *
                                                     (values[6] - values[5]);
                                pts[2][cntloc] = values[7] +
                                                 k / ((NekDouble)(npts3 - 1)) *
                                                     (values[8] - values[7]);
                                cntloc++;
                            }
                            cnt++;
                        }
                    }
                }
            }
            else // give last rank all remaining points
            {
                totpts = totpts - rank * nlocpts;

                for (int i = 0; i < dim; ++i)
                {
                    pts[i] = Array<OneD, NekDouble>(totpts);
                }

                int cnt    = 0;
                int cntloc = 0;

                for (int k = 0; k < npts3; ++k)
                {
                    for (int j = 0; j < npts2; ++j)
                    {
                        for (int i = 0; i < npts1; ++i)
                        {
                            if (cnt >= rank * nlocpts)
                            {
                                pts[0][cntloc] = values[3] +
                                                 i / ((NekDouble)(npts1 - 1)) *
                                                     (values[4] - values[3]);
                                pts[1][cntloc] = values[5] +
                                                 j / ((NekDouble)(npts2 - 1)) *
                                                     (values[6] - values[5]);
                                pts[2][cntloc] = values[7] +
                                                 k / ((NekDouble)(npts3 - 1)) *
                                                     (values[8] - values[7]);
                                cntloc++;
                            }
                            cnt++;
                        }
                    }
                }
            }

            vector<int> ppe;
            ppe.push_back(npts1);
            ppe.push_back(npts2);
            ppe.push_back(npts3);
            m_f->m_fieldPts =
                MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(dim,
                                                                         pts);
            m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsBox);
            m_f->m_fieldPts->SetPointsPerEdge(ppe);
            vector<NekDouble> boxdim;
            boxdim.assign(&values[3], &values[3] + 6);
            m_f->m_fieldPts->SetBoxSize(boxdim);
        }
    }

    FieldSharedPtr fromField = boost::shared_ptr<Field>(new Field());

    std::vector<std::string> files;
    ParseUtils::GenerateOrderedStringVector(m_config["fromxml"].as<string>().c_str(), files);
    // set up session file for from field
    fromField->m_session =
        LibUtilities::SessionReader::CreateInstance(0, 0, files);

    // Set up range based on min and max of local parallel partition
    SpatialDomains::DomainRangeShPtr rng =
        MemoryManager<SpatialDomains::DomainRange>::AllocateSharedPtr();

    int coordim = m_f->m_fieldPts->GetDim();
    int npts    = m_f->m_fieldPts->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > pts;
    m_f->m_fieldPts->GetPts(pts);

    rng->m_checkShape = false;
    rng->m_zmin       = -1;
    rng->m_zmax       = 1;
    rng->m_ymin       = -1;
    rng->m_ymax       = 1;
    switch (coordim)
    {
        case 3:
            rng->m_doZrange = true;
            rng->m_zmin     = Vmath::Vmin(npts, pts[2], 1);
            rng->m_zmax     = Vmath::Vmax(npts, pts[2], 1);
            if (rng->m_zmax == rng->m_zmin)
            {
                rng->m_zmin -= 1;
                rng->m_zmax += 1;
            }
        case 2:
            rng->m_doYrange = true;
            rng->m_ymin     = Vmath::Vmin(npts, pts[1], 1);
            rng->m_ymax     = Vmath::Vmax(npts, pts[1], 1);
        case 1:
            rng->m_doXrange = true;
            rng->m_xmin     = Vmath::Vmin(npts, pts[0], 1);
            rng->m_xmax     = Vmath::Vmax(npts, pts[0], 1);
            break;
        default:
            ASSERTL0(false, "too many values specfied in range");
    }

    // setup rng parameters.
    fromField->m_graph =
        SpatialDomains::MeshGraph::Read(fromField->m_session, rng);

    // Read in local from field partitions
    const SpatialDomains::ExpansionMap &expansions =
        fromField->m_graph->GetExpansions();

    Array<OneD, int> ElementGIDs(expansions.size());
    SpatialDomains::ExpansionMap::const_iterator expIt;

    int i = 0;
    for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
    {
        ElementGIDs[i++] = expIt->second->m_geomShPtr->GetGlobalID();
    }

    // check to see that we do have some elmement in teh domain since
    // possibly all points could be outside of the domain
    ASSERTL0(i > 0, "No elements are set. Are the interpolated points "
                    "wihtin the domain given by the xml files?");

    string fromfld = m_config["fromfld"].as<string>();
    m_f->FieldIOForFile(fromfld)->Import(
        fromfld, fromField->m_fielddef, fromField->m_data,
        LibUtilities::NullFieldMetaDataMap, ElementGIDs);

    int NumHomogeneousDir = fromField->m_fielddef[0]->m_numHomogeneousDir;

    //----------------------------------------------
    // Set up Expansion information to use mode order from field
    fromField->m_graph->SetExpansions(fromField->m_fielddef);

    int nfields = fromField->m_fielddef[0]->m_fields.size();

    fromField->m_exp.resize(nfields);
    fromField->m_exp[0] = fromField->SetUpFirstExpList(NumHomogeneousDir, true);

    m_f->m_exp.resize(nfields);

    // declare auxiliary fields.
    for (i = 1; i < nfields; ++i)
    {
        fromField->m_exp[i] = fromField->AppendExpList(NumHomogeneousDir);
    }

    // load field into expansion in fromfield.
    for (int j = 0; j < nfields; ++j)
    {
        for (i = 0; i < fromField->m_fielddef.size(); i++)
        {
            fromField->m_exp[j]->ExtractDataToCoeffs(
                fromField->m_fielddef[i], fromField->m_data[i],
                fromField->m_fielddef[0]->m_fields[j],
                fromField->m_exp[j]->UpdateCoeffs());
        }
        fromField->m_exp[j]->BwdTrans(fromField->m_exp[j]->GetCoeffs(),
                                      fromField->m_exp[j]->UpdatePhys());

        Array<OneD, NekDouble> newPts(m_f->m_fieldPts->GetNpoints());
        m_f->m_fieldPts->AddField(newPts,
                                  fromField->m_fielddef[0]->m_fields[j]);
    }

    NekDouble clamp_low = m_config["clamptolowervalue"].as<NekDouble>();
    NekDouble clamp_up  = m_config["clamptouppervalue"].as<NekDouble>();
    NekDouble def_value = m_config["defaultvalue"].as<NekDouble>();

    InterpolateFieldToPts(fromField->m_exp, m_f->m_fieldPts, clamp_low,
                          clamp_up, def_value);

    if (!boost::iequals(m_config["cp"].as<string>(), "NotSet"))
    {
        calcCp0();
    }
}

void ProcessInterpPoints::InterpolateFieldToPts(
    vector<MultiRegions::ExpListSharedPtr> &field0,
    LibUtilities::PtsFieldSharedPtr &pts,
    NekDouble clamp_low,
    NekDouble clamp_up,
    NekDouble def_value)
{
    ASSERTL0(pts->GetNFields() >= field0.size(), "ptField has too few fields");

    int nfields = field0.size();

    Interpolator interp;
    if (m_f->m_comm->GetRank() == 0)
    {
        interp.SetProgressCallback(&ProcessInterpPoints::PrintProgressbar,
                                   this);
    }
    interp.Interpolate(field0, pts);
    if (m_f->m_comm->GetRank() == 0)
    {
        cout << endl;
    }

    for (int f = 0; f < nfields; ++f)
    {
        for (int i = 0; i < pts->GetNpoints(); ++i)
        {
            if (pts->GetPointVal(f, i) > clamp_up)
            {
                pts->SetPointVal(f, i, clamp_up);
            }
            else if (pts->GetPointVal(f, i) < clamp_low)
            {
                pts->SetPointVal(f, i, clamp_low);
            }
        }
    }
}

void ProcessInterpPoints::calcCp0()
{
    LibUtilities::PtsFieldSharedPtr pts = m_f->m_fieldPts;
    int dim = pts->GetDim();
    int nq1 = pts->GetNpoints();
    int r, f;
    int pfield = -1;
    NekDouble p0,qinv;
    vector<int> velid;

    vector<NekDouble> values;
    ASSERTL0(ParseUtils::GenerateUnOrderedVector(
                    m_config["cp"].as<string>().c_str(),values),
                "Failed to interpret cp string");

    ASSERTL0(values.size() == 2,
                "cp string should contain 2 values "
                "p0 and q (=1/2 rho u^2)");

    p0  =  values[0];
    qinv = 1.0/values[1];

    for(int i = 0; i < pts->GetNFields(); ++i)
    {
        if(boost::iequals(pts->GetFieldName(i),"p"))
        {
            pfield = i;
        }

        if(boost::iequals(pts->GetFieldName(i),"u")||
            boost::iequals(pts->GetFieldName(i),"v")||
            boost::iequals(pts->GetFieldName(i),"w"))
        {
            velid.push_back(i);
        }
    }

    if(pfield != -1)
    {
        if(!velid.size())
        {
            WARNINGL0(false,"Did not find velocity components for Cp0");
        }
    }
    else
    {
        WARNINGL0(false,"Failed to find 'p' field to determine cp0");
    }

    // Allocate data storage
    Array<OneD, Array< OneD, NekDouble> > data(2);

    for (f = 0; f < 2; ++f)
    {
        data[f] = Array< OneD, NekDouble>(nq1, 0.0);
    }

    for (r = 0; r < nq1; r++)
    {
        if(pfield != -1) // calculate cp
        {
            data[0][r] = qinv*(pts->GetPointVal(dim + pfield, r) - p0);

            if(velid.size()) // calculate cp0
            {
                NekDouble q = 0;
                for(int i = 0; i < velid.size(); ++i)
                {
                    q += 0.5*pts->GetPointVal(dim + velid[i], r)*
                             pts->GetPointVal(dim + velid[i], r);
                }
                data[1][r] = qinv*(pts->GetPointVal(dim + pfield, r)+q - p0);
            }
        }
    }

    if(pfield != -1)
    {
        pts->AddField(data[0], "Cp");
        if(velid.size())
        {
            pts->AddField(data[1], "Cp0");
        }
    }
}

void ProcessInterpPoints::PrintProgressbar(const int position,
                                           const int goal) const
{
    LibUtilities::PrintProgressbar(position, goal, "Interpolating");
}
}
}
