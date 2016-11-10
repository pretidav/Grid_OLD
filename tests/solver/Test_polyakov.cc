   /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_unprec.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG  pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);
  GridSerialRNG    sRNG; 

  LatticeGaugeField Umu(&Grid); 
  //  SU3::TepidConfiguration(pRNG, Umu);  


  NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(std::string("ckpoint_SU3new_lat"),
                                                  std::string("ckpoint_SU3new_rng"), 1);

  int CNFGSTART=29;
  int CNFGEND=29;
  int CNFGSTEP=1;


  for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){   //loop over saved cnfg                                                                                                
    Checkpoint.CheckpointRestore(cnfg,Umu, sRNG, pRNG);
    
#define gaugetransf    
    //Gauge Transformation
#ifdef gaugetransf
      std::vector<int> seeds2({5,6,7,8});      
      GridParallelRNG pRNG2(&Grid);  pRNG2.SeedFixedIntegers(seeds2);
      LatticeColourMatrix Omega(&Grid);
      LatticeColourMatrix ShiftedOmega(&Grid);
      LatticeGaugeField Up(&Grid); Up=zero;
      LatticeColourMatrix U(&Grid); U=zero;
      LatticeColourMatrix Up_mu(&Grid); Up_mu=zero;
      SU<Nc>::LieRandomize(pRNG2, Omega, 1.0);
      for (int mu=0;mu<Nd;mu++){
	U=peekLorentz(Umu,mu);
	ShiftedOmega=Cshift(Omega,mu,1);
	Up_mu=Omega*U*adj(ShiftedOmega);
	pokeLorentz(Up,Up_mu,mu);
	Up_mu=zero;
	U=zero;
      }
      Umu=zero;
      Umu=Up;
#endif

      int T = Grid.GlobalDimensions()[3];
      int X = Grid.GlobalDimensions()[0];
      int Y = Grid.GlobalDimensions()[1];
      int Z = Grid.GlobalDimensions()[2];
      std::cout << "lattice dim = " << X << "x"<< Y <<"x"<< Z <<"x"<< T << std::endl;
      
      std::vector<LatticeColourMatrix> U_peaked(4, &Grid);
      LatticeColourMatrix Ut(&Grid); Ut=zero;
      LatticeColourMatrix P(&Grid); P=zero;
      
      Ut = peekLorentz(Umu,3);
      P = Ut;
      for (int t=1;t<T;t++){ 
	P = PeriodicGimplR::CovShiftForward(Ut,3,P);
      }
      Complex Psum;
      RealD V3=X*Y*Z;
      Psum=sum(trace(P));
      Psum/=V3; 
      std::cout << "tr(Polyakov Loop) = " << TensorRemove(Psum).real() << "     " << TensorRemove(Psum).imag() << std::endl;
      
  } 
  Grid_finalize();
}

