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

template<class d>
struct scal {
  d internal;
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4});


  GridParallelRNG  pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);
  GridSerialRNG    sRNG; 

  LatticeGaugeField Umu(&Grid); 
  //    SU3::TepidConfiguration(pRNG, Umu);  
                                                                                                                                                                                    

  RealD mass=0.166666666666667;


#define POINT_SOURCES

  
  NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(std::string("ckpoint_WG_lat"),
                                                  std::string("ckpoint_WG_rng"), 1);
  
  int CNFGSTART=20;
  int CNFGEND=20;
  int CNFGSTEP=1;

  
  for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){   //loop over saved cnfg                                                                                                
    Checkpoint.CheckpointRestore(cnfg,Umu, sRNG, pRNG);
      
WilsonFermionR Dw(Umu,Grid,RBGrid,mass);
    
#ifdef POINT_SOURCES
    LatticePropagator src(&Grid);
    LatticePropagator Prop(&Grid); Prop=zero;
    LatticeFermion src_vec(&Grid);
    LatticeFermion Mdagsrc_vec(&Grid); Mdagsrc_vec=zero;
    src = zero;
    src_vec = zero;
    std::vector<int> position(4);
    for (int i=0;i<3;i++) position[i]=0; position[3]=0;
    SpinColourMatrix one;
    one=1.0;
    pokeSite(one,src,position);

    std::cout<< "POINT SOURCE in "  << position  <<std::endl;
    for (unsigned int s = 0; s < Ns; s++){
      for (unsigned int c = 0; c < Nc; c++){
	PropToFerm<LatticePropagator,LatticeFermion,qcd>(src_vec,src,s,c);
	std::cout<< "s = " << s << "c = " << c << std::endl;
	ConjugateGradient<LatticeFermion> CG(1.0e-6, 10000);
	LatticeFermion Psi(&Grid);
	MdagMLinearOperator<WilsonFermionR, LatticeFermion> MdagMOp(Dw);
	Psi = zero;
	MdagMOp.AdjOp(src_vec,Mdagsrc_vec);

	CG(MdagMOp, Mdagsrc_vec, Psi); // all-to-one                                                                                                   
	FermToProp<LatticePropagator,LatticeFermion,qcd>(Prop,Psi,s,c);
      }
    }  
#endif
        
#ifdef WALL_SOURCES
    LatticePropagator src(&Grid);
    LatticePropagator Prop(&Grid); Prop=zero;
    LatticeFermion src_vec(&Grid);
    Lattice<iScalar<vInteger> > t(&Grid);
    int Ta,Tb;
    
    Ta=0;Tb=0;
    
    LatticeComplex eta(&Grid);
    LatticeCoordinate(t, 3);
  
    Complex shift(1.,1.);
    Complex I(0,1.0);
   
    //U(1)
    /*
    random(pRNG,eta);
    eta=cos(eta) + I * sin(eta);
    */  
    //Z2 * Z2 
     bernoulli(pRNG, eta);
     eta = (2.*eta - shift)*(1./::sqrt(2.));
    //Gauss
    //    gaussian(pRNG,eta);
    eta = where((t >= Ta) and (t <= Tb), eta, 0.*eta);
    src = 1.;
    src = src*eta;
    if (Ta==Tb){
      std::cout<< "WALL SOURCE in "  << Ta <<std::endl;
    } else {
      std::cout<< "VOLUME SOURCE in "  << Ta << "<->"<< Tb <<std::endl;
    }
    for (unsigned int s = 0; s < Ns; s++){                                                                                                                                       
      for (unsigned int c = 0; c < Nc; c++){                                                                                                                                      
	PropToFerm<LatticePropagator,LatticeFermion,qcd>(src_vec,src,s,c);  
	std::cout<< "s = " << s << "c = " << c << std::endl;  
	ConjugateGradient<LatticeFermion> CG(1.0e-8, 10000);
	LatticeFermion Psi(&Grid);
	MdagMLinearOperator<WilsonFermionR, LatticeFermion> MdagMOp(Dw);
	FermToProp<LatticePropagator,LatticeFermion,qcd>(Prop,Psi,s,c);
	Psi = zero;
	CG(MdagMOp, src_vec, Psi); // all-to-one                                                                                                    
	FermToProp<LatticePropagator,LatticeFermion,qcd>(Prop,Psi,s,c);  
      }
    }      
#endif
    
#ifdef ALL2ALL_SOURCES
    
    int Nnoise = 12; 
    LatticePropagator src(&Grid);
    LatticePropagator Prop(&Grid); Prop = zero;
    LatticePropagator Prop2(&Grid); Prop2 = zero;
    LatticeFermion src_vec(&Grid);
    Lattice<iScalar<vInteger> > t(&Grid);
    int Ta,Tb;
    
    Ta=0;Tb=19;
    
    LatticeComplex eta(&Grid); 
    LatticeCoordinate(t, 3);
    
    Complex shift(1.,1.);
    Complex I(0,1.0);
    for (int nsrc=0;nsrc<Nnoise;nsrc++){
      //U(1)
      /*
	random(pRNG,eta);
	eta=cos(eta) + I * sin(eta);
      */
      //Z2 * Z2
      bernoulli(pRNG, eta);
      eta = (2.*eta - shift)*(1./::sqrt(2.));
      //Gauss
      //    gaussian(pRNG,eta);
      eta = where((t >= Ta) and (t <= Tb), eta, 0.*eta);
      src = 1.;
      src = src*eta; // * only defined if eta complex
      if (Ta==Tb){
	std::cout<< "WALL SOURCE in "  << Ta <<std::endl; // you want all-to-all
      } else { 
	std::cout<< "VOLUME SOURCE n"<<nsrc<<" in "  << Ta << "<->"<< Tb <<std::endl;
      }  

      for (unsigned int s = 0; s < Ns; s++){
	for (unsigned int c = 0; c < Nc; c++){
	  PropToFerm<LatticePropagator,LatticeFermion,qcd>(src_vec,src,s,c);
	  std::cout<< "s = " << s << "c = " << c << std::endl;
	  ConjugateGradient<LatticeFermion> CG(1.0e-8, 10000);
	  LatticeFermion Psi(&Grid);
	  MdagMLinearOperator<WilsonFermionR, LatticeFermion> MdagMOp(Dw);	  
	  Psi = zero;
	  CG(MdagMOp, src_vec, Psi);   
	  FermToProp<LatticePropagator,LatticeFermion,qcd>(Prop2,Psi,s,c);
	}
      }
      Prop=Prop + Prop2*adj(src);
    }
    RealD norm=1./Nnoise;
    Prop=Prop*norm;
#endif
    
        
    const Gamma::GammaMatrix *g = Gamma::GammaMatrices;
    const char **list           = Gamma::GammaMatrixNames;
    
    //Gamma(g[mu]); how to do product of gamma matrices ??? 
    
    int T = Grid.GlobalDimensions()[3];
    int X = Grid.GlobalDimensions()[0];
    int Y = Grid.GlobalDimensions()[1];
    int Z = Grid.GlobalDimensions()[2];

    std::cout << "lattice dim = " << X << "x"<< Y <<"x"<< Z <<"x"<< T << std::endl;
    
    int mu;
    int nu;

    LatticeComplex C(&Grid);
    std::vector <TComplex> Ct;

    for (mu=0;mu<6;mu++){
      for (nu=0;nu<6;nu++){
	C = zero;
	C = trace(Gamma(g[mu])*Prop*Gamma(g[nu])*Gamma(g[5])*(adj(Prop)*Gamma(g[5])));
	sliceSum(C,Ct,3);
	std::cout << "t     <" << list[mu] << " x "<< list[nu] << ">" << std::endl;
	for (int t=0;t<T;t++){
	  std::cout << " " << t << "     " << TensorRemove(Ct[t])  << std::endl;
	}
      }
    }

    for (mu=0;mu<6;mu++){
      for (nu=1;nu<5;nu++){
        C = zero;
        C =  trace(Gamma(g[mu])*Prop*Gamma(g[nu])*(adj(Prop)*Gamma(g[5])));
        sliceSum(C,Ct,3);
	std::cout << "t     <" << list[mu] << " x "<< list[nu+11] << ">" << std::endl;
        for (int t=0;t<T;t++){
	  std::cout << " " << t << "     " << TensorRemove(Ct[t])  << std::endl;
        }
      }
    }

  }    
  Grid_finalize();
}

