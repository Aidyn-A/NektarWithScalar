///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFE.cpp
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
// Description: Navier Stokes equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/NavierStokesCFE.h>

using namespace std;

namespace Nektar
{
    string NavierStokesCFE::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "NavierStokesCFE", NavierStokesCFE::create,
            "NavierStokes equations in conservative variables.");

    NavierStokesCFE::NavierStokesCFE(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleFlowSystem(pSession)
    {
    }

    NavierStokesCFE::~NavierStokesCFE()
    {

    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void NavierStokesCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

        m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance(diffName, diffName);

        m_diffusion->SetFluxVectorNS(&NavierStokesCFE::
                                      GetViscousFluxVector, this);

        // Concluding initialisation of diffusion operator
        m_diffusion->InitObject         (m_session, m_fields);
    }

    void NavierStokesCFE::v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();
        int nScalars  = nvariables-(2+m_spacedim);

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

        Array<OneD, Array<OneD, NekDouble> > scalars(nScalars);
        Array<OneD, Array<OneD, NekDouble> > scalarsFwd(nScalars);
        Array<OneD, Array<OneD, NekDouble> > scalarsBwd(nScalars);

        Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables-1);
        Array<OneD, Array<OneD, NekDouble> > inFwd(nvariables-1);
        Array<OneD, Array<OneD, NekDouble> > inBwd(nvariables-1);

        for (i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints);
        }

        for (i = 0; i < nvariables-1; ++i)
        {
            inarrayDiff[i] = Array<OneD, NekDouble>(npoints);
            inFwd[i]       = Array<OneD, NekDouble>(nTracePts);
            inBwd[i]       = Array<OneD, NekDouble>(nTracePts);
        }

        for (i = 0; i < nScalars; ++i)
        {
            scalars[i] = Array<OneD, NekDouble>(npoints);
            scalarsFwd[i] = Array<OneD, NekDouble>(nTracePts);
            scalarsBwd[i] = Array<OneD, NekDouble>(nTracePts);
        }
        // Extract pressure
        //    (use inarrayDiff[0] as a temporary storage for the pressure)
        m_varConv->GetPressure(inarray, inarrayDiff[0]);

        // Extract temperature
        m_varConv->GetTemperature(inarray, inarrayDiff[0],
                inarrayDiff[m_spacedim]);

        // Extract velocities
        m_varConv->GetVelocityVector(inarray, inarrayDiff);

        // Extract scalars
        m_varConv->GetScalars(inarray, scalars);
        for (i = 0; i < nScalars; ++i)
        {
            Vmath::Vcopy(npoints, scalars[i], 1,
                                  inarrayDiff[m_spacedim+1+i], 1);
        }       

        // Repeat calculation for trace space
        if (pFwd == NullNekDoubleArrayofArray || 
            pBwd == NullNekDoubleArrayofArray)
        {
            inFwd = NullNekDoubleArrayofArray;
            inBwd = NullNekDoubleArrayofArray;
        }
        else
        {
            m_varConv->GetPressure(pFwd,    inFwd[0]);
            m_varConv->GetPressure(pBwd,    inBwd[0]);

            m_varConv->GetTemperature(pFwd,    inFwd[0],
                inFwd[m_spacedim]);
            m_varConv->GetTemperature(pBwd,    inBwd[0],
                inBwd[m_spacedim]);

            m_varConv->GetVelocityVector(pFwd, inFwd);
            m_varConv->GetVelocityVector(pBwd, inBwd);

            m_varConv->GetScalars(pFwd, scalarsFwd);
            m_varConv->GetScalars(pBwd, scalarsBwd);

            for (i = 0; i < nScalars; ++i)
            {
                Vmath::Vcopy(nTracePts, scalarsFwd[i], 1,
                                      inFwd[m_spacedim+1+i], 1);
                Vmath::Vcopy(nTracePts, scalarsBwd[i], 1,
                                      inBwd[m_spacedim+1+i], 1);
            }
        }

        // Diffusion term in physical rhs form
        m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff,
                             inFwd, inBwd);

        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vadd(npoints,
                        outarrayDiff[i], 1,
                        outarray[i], 1,
                        outarray[i], 1);
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > diffusivity        (nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        thermalConductivity, 1);
        }

        Vmath::Smul(nPts, 1.0/m_Schmidt, mu, 1, diffusivity, 1);

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, derivativesO1[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, derivativesO1[i][j], 1,
                                  derivativesO1[j][i], 1,
                                  viscousTensor[i][j+1], 1);

                Vmath::Vmul(nPts, mu, 1,
                                  viscousTensor[i][j+1], 1,
                                  viscousTensor[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, viscousTensor[i][j+1], 1,
                                  divVel, 1,
                                  viscousTensor[i][j+1], 1);
                }
                else
                {
                    // Copy to make symmetric
                    Vmath::Vcopy(nPts, viscousTensor[i][j+1], 1,
                                       viscousTensor[j][i+1], 1);
                }
            }
        }

        // Terms for the energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, physfield[j], 1,
                               viscousTensor[i][j+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, thermalConductivity, 1,
                               derivativesO1[i][m_spacedim], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
        }

        // Terms for the scalars
        for (i = 0; i < m_spacedim; ++i)
        {
            for(j=m_spacedim+2; j<nVariables; ++j)
            {
                Vmath::Zero(nPts, viscousTensor[i][j], 1);

                // Add D*Yj_i
                Vmath::Vvtvp(nPts, diffusivity, 1,
                                   derivativesO1[i][j-1], 1,
                                   viscousTensor[i][j], 1,
                                   viscousTensor[i][j], 1);

            }
        }
    }
}
