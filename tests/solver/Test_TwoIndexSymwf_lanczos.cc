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

RealD AllZero(RealD x){ return 0.;}



typedef WilsonTwoIndexSymmetricFermionR FermionAction; // type of lattice fermions (Wilson, DW, ...)                                                                             
typedef typename FermionAction::FermionField FermionField;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid;
  GridRedBlackCartesian * FrbGrid;

  FGrid=UGrid;
  FrbGrid=UrbGrid;

  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n",UGrid,UrbGrid,FGrid,FrbGrid);

  LatticeGaugeField Umu(UGrid);
  TwoIndexSymmetricRepresentation HiRep(UGrid);
  
  NerscField header;                                                          
  
  int CNFGSTART,CNFGEND,CNFGSTEP;
  
  CNFGSTART=2999;
  CNFGEND=3000;
  CNFGSTEP=1;
  
  for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){
    
    std::vector<int>         seeds5({cnfg,cnfg+1,cnfg+2,cnfg+3});
    GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
    GridParallelRNG          RNG5rb(FrbGrid);  RNG5.SeedFixedIntegers(seeds5);

    std::string file =  "./ckpoint_lat." +  std::to_string(cnfg);
    std::cout << file << std::endl;

    NerscIO::readConfiguration(Umu,header,file);
    
    HiRep.update_representation(Umu);  
    typename TwoIndexRep<Nc, Symmetric>::LatticeField Umu_R = HiRep.U;
    
    RealD mass=-1.32; // am 
 
    FermionAction DW(Umu_R, *FGrid, *FrbGrid, mass);
    
    MdagMLinearOperator<FermionAction,FermionField> HermOp(DW);
    

    const int Nstop = 30;
    const int Nk = 60;
    const int Np = 60;
    const int Nm = Nk+Np;
    const int MaxIt= 10000;
    RealD resid = 1.0e-6;
    
    std::vector<double> Coeffs { 0.,-1. };
    Polynomial<FermionField> PolyX(Coeffs);
    Chebyshev<FermionField> Cheb(0.2,5.,11);//0.2,5.0,11
    
    //  ChebyshevLanczos<LatticeFermion> Cheb(9.,1.,0.,20);
    //  Cheb.csv(std::cout);
    //  exit(-24);
    ImplicitlyRestartedLanczos<FermionField> IRL(HermOp,Cheb,Nstop,Nk,Nm,resid,MaxIt);
  
  
    std::vector<RealD> eval(Nm);
    FermionField src(FGrid); 
    gaussian(RNG5rb,src);
    std::vector<FermionField> evec(Nm,FGrid);

    for(int i=0;i<1;i++){
      std::cout << i<<" / "<< Nm<< " grid pointer "<<evec[i]._grid<<std::endl;
    }
    
    int Nconv;
    
    IRL.calc(eval,evec,src,Nconv); 
    
    std::cout << " -- Smallest eval = " << eval[Nconv-1] <<std::endl;
  
  }
  Grid_finalize();
}
