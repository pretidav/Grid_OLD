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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
 
  NerscField header;                                                          
  LatticeGaugeField Umu(UGrid);
  AdjointRep<Nc> AdjRep(UGrid);  


  int CNFGSTART,CNFGEND,CNFGSTEP;
  
  CNFGSTART=2999;
  CNFGEND=3000;
  CNFGSTEP=1;
  
  GridSerialRNG sRNG;
  GridParallelRNG pRNG(UGrid);
  int cnfg=CNFGSTART;
  
  NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(std::string("ckpoint_lat"),
						  std::string("ckpoint_rng"), 1);

  //for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){
    
    /*
    std::string file =  "./ckpoint_lat_SU3." +  std::to_string(cnfg);
    
    std::cout << file << std::endl;
    
    NerscIO::readConfiguration(Umu,header,file);
  

    std::cout << "header.checksum = " << header.checksum << std::endl;   // checksum from cnfg is wrong, not header one.
    */
  
    
     
      Checkpoint.CheckpointRestore(2999,Umu, sRNG, pRNG);
    
     
      Checkpoint.CheckpointRestore(3000,Umu, sRNG, pRNG);

    
    //    AdjRep.update_representation(Umu);  
    //  typename AdjointRep<Nc>::LatticeField Umu_Adj = AdjRep.U;
   

      //}
  Grid_finalize();
}
