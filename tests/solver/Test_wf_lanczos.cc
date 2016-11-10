    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_lanczos.cc

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
#include <Grid/Grid.h>


using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  LatticeGaugeField Umu(UGrid);
  

  
  //NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(std::string("ckpoint_SU3new_lat"),
  //                                                std::string("ckpoint_SU3new_rng"), 1);
// int CNFGSTART=29;
//  int CNFGEND=29;
//  int CNFGSTEP=1;
  
  //  for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){
  std::vector<int>         seeds5({0,1,2,3});
  GridParallelRNG          RNG5(UGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG5rb(UrbGrid);  RNG5.SeedFixedIntegers(seeds5);
  SU<Nc>::TepidConfiguration(RNG5, Umu);  

    RealD mass=-2.0; // am 
    WilsonFermionR D(Umu, *UGrid, *UrbGrid, mass);
    
    //        MdagMLinearOperator<WilsonFermionR,LatticeFermion> g5D(D);
    //        SchurDiagTwoOperator<WilsonFermionR,LatticeFermion> g5D(D);
    Gamma5R5HermitianLinearOperator<WilsonFermionR,LatticeFermion> g5D(D); //FIXME: I DO NOT CONVERGE 
   
    const int Nstop = 20;
    const int Nk = 60;
    const int Np = 60;
    const int Nm = Nk+Np;
    const int MaxIt= 10000;
    RealD resid = 1.0e-0;
    
    std::vector<double> Coeffs { 0.,1. };
    Polynomial<LatticeFermion> PolyX(Coeffs);
    Chebyshev<LatticeFermion> Cheb(0.01,15.,20);//0.2,5.0,11
    
    //  ChebyshevLanczos<LatticeFermion> Cheb(9.,1.,0.,20);
    //  Cheb.csv(std::cout);
    //  exit(-24);
    ImplicitlyRestartedLanczos<LatticeFermion> IRL(g5D,Cheb,Nstop,Nk,Nm,resid,MaxIt);
  
    std::vector<RealD> eval(Nm);
    LatticeFermion src(UGrid); 
    gaussian(RNG5,src);
    std::vector<LatticeFermion> evec(Nm,UGrid);

    for(int i=0;i<1;i++){
      std::cout << i<<" / "<< Nm<< " grid pointer "<<evec[i]._grid<<std::endl;
    }
    
    int Nconv;
    IRL.calc(eval,evec,src,Nconv); 
    
    std::cout << " -- Smallest eval = " << eval[Nconv-1] <<std::endl;
  
    //  }
  Grid_finalize();
}
