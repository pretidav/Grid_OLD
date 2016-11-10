/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonAdjointFermionGauge.cc

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include "Grid/Grid.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid {
namespace QCD {

// Here change the allowed (higher) representations
typedef Representations< FundamentalRepresentation, TwoIndexAntiSymmetricRepresentation > TheRepresentations;


class HmcRunner : public NerscHmcRunnerHirep< TheRepresentations > {
 public:
  void BuildTheAction(int argc, char **argv)

  {
    typedef WilsonTwoIndexAntiSymmetricImplR ImplPolicy; // gauge field implemetation for the pseudofermions
    typedef WilsonTwoIndexAntiSymmetricFermionR FermionAction; // type of lattice fermions (Wilson, DW, ...)
    typedef typename FermionAction::FermionField FermionField;

    UGrid = SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
        GridDefaultMpi());
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

    FGrid = UGrid;
    FrbGrid = UrbGrid;

    // temporarily need a gauge field
    TwoIndexAntiSymmetricRepresentation::LatticeField U(UGrid);

    // Gauge action
    WilsonGaugeActionR Waction(10.2);

    Real mass = -0.04743083;//
    FermionAction FermOp(U, *FGrid, *FrbGrid, mass);

    ConjugateGradient<FermionField> CG(1.0e-10, 10000, false);

    // Pass two solvers: one for the force computation and one for the action
    TwoFlavourPseudoFermionAction<ImplPolicy> Nf2(FermOp, CG, CG);

    // Set smearing (true/false), default: false
    Nf2.is_smeared = false;

    // Collect actions
    ActionLevel<LatticeGaugeField, TheRepresentations > Level1(1);
    Level1.push_back(&Nf2);

    ActionLevel<LatticeGaugeField, TheRepresentations > Level2(4);
    Level2.push_back(&Waction);

    TheAction.push_back(Level1);
    TheAction.push_back(Level2);

    Run(argc, argv);
  };
};
}
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads
            << " threads" << std::endl;

  HmcRunner TheHMC;

  TheHMC.BuildTheAction(argc, argv);
}
