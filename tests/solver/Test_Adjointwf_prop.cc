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


typedef WilsonAdjFermionR FermionOp; // type of lattice fermions (Wilson, DW, ...)
typedef typename FermionOp::FermionField FermionField;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = UGrid;
  GridRedBlackCartesian * FrbGrid = UrbGrid;
 

  GridSerialRNG sRNG;
  GridParallelRNG pRNG(UGrid);


  LatticeGaugeField Umu(UGrid);
  AdjointRep<Nc> HiRep(UGrid);  
  
  int CNFGSTART,CNFGEND,CNFGSTEP;
  
  CNFGSTART=1;
  CNFGEND=1;
  CNFGSTEP=1;
  


  RealD mass=0.0; // am 


  
  NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(std::string("ckpoint_lat_SU3Fund"),
						  std::string("ckpoint_rng_SU3Fund"), 1);

  for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){   //loop over saved cnfg 
    Checkpoint.CheckpointRestore(cnfg,Umu, sRNG, pRNG);
    HiRep.update_representation(Umu);
    typename AdjointRep<Nc>::LatticeField Umu_R = HiRep.U;

    FermionOp Dw(Umu_R, *FGrid, *FrbGrid, mass);

    std::vector<int>         seeds_src({cnfg,cnfg+1,cnfg+2,cnfg+3}); //do something bettere with seeds
    GridParallelRNG          RNGsrc(UGrid);  RNGsrc.SeedFixedIntegers(seeds_src);

    
    //POINT SOURCE
    typename AdjointRep<Nc>::LatticePropagator src(UGrid);
    typename AdjointRep<Nc>::LatticePropagator Prop(UGrid);
    FermionField src_vec(UGrid);
    src = zero;
    src_vec = zero;
    std::vector<int> position(4);
    for (int i=0;i<3;i++) position[i]=0; position[3]=0;
    typename AdjointRep<Nc>::SpinColourMatrix one;
    one=1.0;
    pokeSite(one,src,position);

    std::cout<< "POINT SOURCE in "  << position  <<std::endl;

    MdagMLinearOperator<WilsonAdjFermionR, FermionField> MdagMOp(Dw);
    for (unsigned int s = 0; s < Ns; s++){
      for (unsigned int c = 0; c < Nc*Nc-1; c++){
	PropToFerm< AdjointRep<Nc>::LatticePropagator,FermionField, AdjointRep<Nc>::parms>(src_vec,src,s,c);
	std::cout<< "s = " << s << "c = " << c << std::endl;
	ConjugateGradient<FermionField> CG(1.0e-8, 10000);
	FermionField Psi(UGrid);
	Psi = zero;
	CG(MdagMOp, src_vec, Psi); // one-to-all                                                                                                                                
	FermToProp< AdjointRep<Nc>::LatticePropagator,FermionField, AdjointRep<Nc>::parms>(Prop,Psi,s,c);
      }
    }

    //PP
    LatticeComplex C(UGrid);
    std::vector <TComplex> Ct;
    C = trace(Prop*adj(Prop));
    sliceSum(C,Ct,3);
    for (int t=0;t<20;t++){
      std::cout << "t     <PP>" << std::endl;
      std::cout << " " << t << "     " << TensorRemove(Ct[t]) << std::endl;
    }


  }
  Grid_finalize();
}
