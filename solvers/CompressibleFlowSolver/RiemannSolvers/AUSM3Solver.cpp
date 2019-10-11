///////////////////////////////////////////////////////////////////////////////
//
// File: AUSM3Solver.cpp
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
// Description: AUSM3 Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/AUSM3Solver.h>

namespace Nektar
{
    std::string AUSM3Solver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "AUSM3",
            AUSM3Solver::create,
            "AUSM3 Riemann solver");

    AUSM3Solver::AUSM3Solver() : CompressibleSolver()
    {

    }

    /**
     * @brief AUSM3 Riemann solver
     *
     * @param rhoL      Density left state.
     * @param rhoR      Density right state.  
     * @param rhouL     x-momentum component left state.  
     * @param rhouR     x-momentum component right state.  
     * @param rhovL     y-momentum component left state.  
     * @param rhovR     y-momentum component right state.  
     * @param rhowL     z-momentum component left state.  
     * @param rhowR     z-momentum component right state.
     * @param EL        Energy left state.  
     * @param ER        Energy right state. 
     * @param rhof      Computed Riemann flux for density.
     * @param rhouf     Computed Riemann flux for x-momentum component 
     * @param rhovf     Computed Riemann flux for y-momentum component 
     * @param rhowf     Computed Riemann flux for z-momentum component 
     * @param Ef        Computed Riemann flux for energy.
     */
    void AUSM3Solver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {
        static NekDouble gamma = m_params["gamma"]();

        // Left and Right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;

        // Pressure
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));

        // Speeds of sound
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        
        // Average speeds of sound
        NekDouble cA = 0.5 * (cL + cR);
        
        // Local Mach numbers
        NekDouble ML = uL / cA;
        NekDouble MR = uR / cA;
        
        // Parameters for specify the upwinding
        // Note: if fa = 1 then AUSM3 = AUSM3
        NekDouble Mco    = 0.01;
        NekDouble Mtilde = 0.5 * (ML * ML + MR * MR);
        NekDouble Mo     = std::min(1.0, std::max(Mtilde, Mco*Mco));
        NekDouble fa     = Mo * (2.0 - Mo);
        NekDouble beta   = 0.125;
        NekDouble alpha  = 0.1875;
        NekDouble sigma  = 1.0;
        NekDouble Kp     = 0.25;
        NekDouble Ku     = 0.75;
        NekDouble rhoA   = 0.5 * (rhoL + rhoR);
        NekDouble Mp     = -(Kp / fa) * ((pR - pL) / (rhoA * cA * cA)) * 
                            std::max(1.0 - sigma * Mtilde, 0.0);
        
        NekDouble Mbar   = M4Function(0, beta, ML) + 
                           M4Function(1, beta, MR) + Mp;
        
        NekDouble pu     = -2.0 * Ku * rhoA * cA * cA * (MR - ML) * 
                           P5Function(0, alpha, ML) * P5Function(1, alpha, MR);
        
        NekDouble pbar   = pL * P5Function(0, alpha, ML) + 
                           pR * P5Function(1, alpha, MR) + pu;
        
        if (Mbar >= 0.0)
        {
            rhof  = cA * Mbar * rhoL;
            rhouf = cA * Mbar * rhoL * uL + pbar;
            rhovf = cA * Mbar * rhoL * vL;
            rhowf = cA * Mbar * rhoL * wL;
            Ef    = cA * Mbar * (EL + pL);
        }
        else
        {
            rhof  = cA * Mbar * rhoR;
            rhouf = cA * Mbar * rhoR * uR + pbar;
            rhovf = cA * Mbar * rhoR * vR;
            rhowf = cA * Mbar * rhoR * wR;
            Ef    = cA * Mbar * (ER + pR);
        }
    }


    void AUSM3Solver::v_ArraySolve(
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux,
        const int nDim)
    {
        static NekDouble gamma = m_params["gamma"]();
        for (int i = 0; i < Fwd[0].num_elements(); i++)
        {
            double  rhoL = 0; 
            double  rhouL = 0; 
            double  rhovL = 0; 
            double  rhowL = 0; 
            double  EL = 0;
        
            double  rhoR = 0; 
            double  rhouR = 0; 
            double  rhovR = 0; 
            double  rhowR = 0; 
            double  ER = 0;

            if(nDim == 1)
            {
                rhoL = Fwd[0][i]; 
                rhouL = Fwd[1][i]; 
                EL = Fwd[2][i];

                        
                rhoR = Bwd[0][i]; 
                rhouR = Bwd[1][i]; 
                ER = Bwd[2][i];
            }
            else if(nDim == 2)
            {
                rhoL = Fwd[0][i]; 
                rhouL = Fwd[1][i]; 
                rhovL = Fwd[2][i]; 
                EL = Fwd[3][i];

                        
                rhoR = Bwd[0][i]; 
                rhouR = Bwd[1][i]; 
                rhovR = Bwd[2][i]; 
                ER = Bwd[3][i];
            }
            else if(nDim == 3)
            {
                rhoL = Fwd[0][i]; 
                rhouL = Fwd[1][i]; 
                rhovL = Fwd[2][i]; 
                rhowL = Fwd[3][i]; 
                EL = Fwd[4][i];

                        
                rhoR = Bwd[0][i]; 
                rhouR = Bwd[1][i]; 
                rhovR = Bwd[2][i]; 
                rhowR = Bwd[3][i]; 
                ER = Bwd[4][i];
            }

            // Left and Right velocities
            NekDouble uL = rhouL / rhoL;
            NekDouble vL = rhovL / rhoL;
            NekDouble wL = rhowL / rhoL;
            NekDouble uR = rhouR / rhoR;
            NekDouble vR = rhovR / rhoR;
            NekDouble wR = rhowR / rhoR;

            // Pressure
            NekDouble pL = (gamma - 1.0) *
                (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
            NekDouble pR = (gamma - 1.0) *
                (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));

            // Speeds of sound
            NekDouble cL = sqrt(gamma * pL / rhoL);
            NekDouble cR = sqrt(gamma * pR / rhoR);
            
            // Average speeds of sound
            NekDouble cA = 0.5 * (cL + cR);
            
            // Local Mach numbers
            NekDouble ML = uL / cA;
            NekDouble MR = uR / cA;
            
            // Parameters for specify the upwinding
            // Note: if fa = 1 then AUSM3 = AUSM2
            NekDouble Mco    = 0.01;
            NekDouble Mtilde = 0.5 * (ML * ML + MR * MR);
            NekDouble Mo     = sqrt(std::min(1.0, std::max(Mtilde, Mco*Mco)));
            NekDouble fa     = Mo * (2.0 - Mo);
            NekDouble beta   = 0.125;
            NekDouble alpha  = 0.1875 * (5.0 * fa * fa - 4.0);
            NekDouble sigma  = 1.0;
            NekDouble Kp     = 0.25;
            NekDouble Ku     = 0.75;
            NekDouble rhoA   = 0.5 * (rhoL + rhoR);
            NekDouble Mp     = -(Kp / fa) * ((pR - pL) / (rhoA * cA * cA)) * 
                                std::max(1.0 - sigma * Mtilde, 0.0);
            
            NekDouble Mbar   = M4Function(0, beta, ML) + 
                               M4Function(1, beta, MR) + Mp;
            
            NekDouble pu     = -2.0 * Ku * rhoA * cA * fa * (uR - uL) * 
                               P5Function(0, alpha, ML) * P5Function(1, alpha, MR);

            NekDouble pbar   = pL * P5Function(0, alpha, ML) + 
                               pR * P5Function(1, alpha, MR) + pu;

            if (Mbar >= 0.0)
            {
                NekDouble mdot = cA * Mbar * rhoL;
                flux[0][i] = mdot;

                for(int k=1; k<Fwd.num_elements(); k++)
                {
                    flux[k][i] = mdot * Fwd[k][i]/Fwd[0][i];
                }

                flux[1][i] = mdot * Fwd[1][i]/Fwd[0][i] + pbar;

                flux[nDim+1][i] = mdot * (Fwd[nDim+1][i] + pL)/Fwd[0][i];
            }
            else
            {
                NekDouble mdot = cA * Mbar * rhoR;
                flux[0][i] = mdot;

                for(int k=1; k<Bwd.num_elements(); k++)
                {
                    flux[k][i] = mdot * Bwd[k][i]/Bwd[0][i];
                }

                flux[1][i] = mdot * Bwd[1][i]/Bwd[0][i] + pbar;

                flux[nDim+1][i] = mdot * (Bwd[nDim+1][i] + pR)/Bwd[0][i];
            }
        }
    }
    
    double AUSM3Solver::M1Function(int A, double M)
    {
        double out;
        
        if (A == 0)
        {
            out = 0.5 * (M + fabs(M));
        }
        else 
        {
            out = 0.5 * (M - fabs(M));
        }
        
        return out; 
    }
    
    double AUSM3Solver::M2Function(int A, double M)
    {
        double out;
        
        if (A == 0)
        {
            out = 0.25 * (M + 1.0) * (M + 1.0);
        }
        else
        {
            out = -0.25 * (M - 1.0) * (M - 1.0);
        }
        
        return out;
    }
    
    double AUSM3Solver::M4Function(int A, double beta, double M)
    {
        double out;
        
        if (fabs(M) >= 1.0)
        {
            out = M1Function(A, M);
        }
        else
        {
            out = M2Function(A, M);
            
            if (A == 0)
            {
                out *= 1.0 - 16.0 * beta * M2Function(1, M);
            }
            else
            {
                out *= 1.0 + 16.0 * beta * M2Function(0, M);
            }
        }
        
        return out;
    }
    
    double AUSM3Solver::P5Function(int A, double alpha, double M)
    {
        double out;
        
        if (fabs(M) >= 1.0)
        {
            out = (1.0 / M) * M1Function(A, M);
        }
        else
        {
            out = M2Function(A, M);
            
            if (A == 0)
            {
                out *= (2.0 - M) - 16.0 * alpha * M * M2Function(1, M);
            }
            else
            {
                out *= (-2.0 - M) + 16.0 * alpha * M * M2Function(0, M);
            }
        }
        
        return out;
    }
}
