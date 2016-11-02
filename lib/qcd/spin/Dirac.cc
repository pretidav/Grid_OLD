    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/spin/Dirac.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid.h>

namespace Grid {

  namespace QCD {

    Gamma::GammaMatrix  Gamma::GammaMatrices [] = {
      Gamma::Identity,
      Gamma::GammaX,
      Gamma::GammaY,
      Gamma::GammaZ,
      Gamma::GammaT,
      Gamma::Gamma5,
      Gamma::MinusIdentity,
      Gamma::MinusGammaX,
      Gamma::MinusGammaY,
      Gamma::MinusGammaZ,
      Gamma::MinusGammaT,
      Gamma::MinusGamma5,
      Gamma::GammaXGamma5,
      Gamma::GammaYGamma5,
      Gamma::GammaZGamma5,
      Gamma::GammaTGamma5,
      Gamma::SigmaXY,
      Gamma::SigmaXZ,
      Gamma::SigmaYZ,
      Gamma::SigmaXT,
      Gamma::SigmaYT,
      Gamma::SigmaZT
    };
    const char *Gamma::GammaMatrixNames[] = { 
      "Identity ",
      "GammaX   ",
      "GammaY   ",
      "GammaZ   ",
      "GammaT   ",
      "Gamma5   ",
      "-Identity",
      "-GammaX  ",
      "-GammaY  ",
      "-GammaZ  ",
      "-GammaT  ",
      "-Gamma5  ",
      "GammaXGamma5",
      "GammaYGamma5",
      "GammaZGamma5",
      "GammaTGamma5",
      "SigmaXY",
      "SigmaXZ",
      "SigmaYZ",
      "SigmaXT",
      "SigmaYT",
      "SigmaZT",
      "         "
    };
    
    SpinMatrix makeGammaProd(const unsigned int i)
    {
      SpinMatrix g;
      
      g = 1.;
      if (i & 0x1)
      {
        g = g*Gamma(Gamma::GammaMatrix::GammaX);
      }
      if (i & 0x2)
      {
        g = g*Gamma(Gamma::GammaMatrix::GammaY);
      }
      if (i & 0x4)
      {
        g = g*Gamma(Gamma::GammaMatrix::GammaZ);
      }
      if (i & 0x8)
      {
        g = g*Gamma(Gamma::GammaMatrix::GammaT);
      }
      
      return g;
    }

    //    void sprojMul( vHalfSpinColourVector &out,vColourMatrix &u, vSpinColourVector &in){
    //      vHalfSpinColourVector hspin;
    //      spProjXp(hspin,in);
    //      mult(&out,&u,&hspin);
    //    }
  }
}
