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

inline RealD Sign(int a, int b){ return (a==b) ? 0:(b-a)/fabs(b-a);}

typedef WilsonTwoIndexAntiSymmetricFermionR HirepFermionOp;
typedef typename HirepFermionOp::FermionField HirepFermionField;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout<<"FIXME: I AM NOT GAUGE INVARIANT!!!"<<std::endl;

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4});


  GridParallelRNG  pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);
  GridSerialRNG    sRNG; 

  LatticeGaugeField Umu(&Grid); 
  TwoIndexAntiSymmetricRepresentation HiRep(&Grid);


  RealD m_ud=-0.2;
  RealD m_s=-0.2;


#define POINT_SOURCE
  int Nperm=12; //Nc!/2   
  int ep[12][Nc];
  int em[12][Nc];
  int contp=0,contm=0;
  int S1,S2,S3,S4,S5,S6;
  // qqQ with q in fund SU(4) and Q in twoindex anti-sym SU(4)  
  // \epsilon_{ijkl}                                                                                                                                                                
  for (int i=0;i<Nc;i++){
    for (int j=0;j<Nc;j++){
      for (int k=0;k<Nc;k++){
	for (int l=0;l<Nc;l++){
	  S1=Sign(i,j);
	  S2=Sign(i,k);
	  S3=Sign(i,l);
	  S4=Sign(j,k);
	  S5=Sign(j,l);
	  S6=Sign(k,l);
	  if ((S1*S2*S3*S4*S5*S6)>0){
	    ep[contp][0]=i;
	    ep[contp][1]=j;
	    ep[contp][2]=k;
	    ep[contp][3]=l;
	    contp++;
	  } else if ((S1*S2*S3*S4*S5*S6)<0){
	    em[contm][0]=i;
	    em[contm][1]=j;
	    em[contm][2]=k;
	    em[contm][3]=l;
	    contm++;
	  }
	}
      }
    }
  }
  
  //Mapping (i,j) -> A , with i,j=1,...,4 and A=1,...,6
  int A[Nc][Nc];
  int cont=0;
  for (int a=1;a<Nc;a++){
    for (int b=0;b<a;b++){
      A[a][b]=cont;
      A[b][a]=cont;
      cont++;
    }
  }
  

  /*
  for (int i=1;i<Nc;i++){
    for (int j=0;j<i;j++){
  std::cout << "A[" << i << "]["<< j <<"]= " << A[i][j] << std::endl;
    }
  }


  for (int i=0;i<12;i++){
    for (int j=0;j<4;j++){
      std::cout << "ep[" << i << "][" << j << "]=" << ep[i][j] << std::endl;
    }
  }

  */

  NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(std::string("ckpoint_SU4_lat"),
                                                  std::string("ckpoint_SU4_rng"), 1);

  int CNFGSTART=33;
  int CNFGEND=33;
  int CNFGSTEP=1;


  const Gamma::GammaMatrix *g = Gamma::GammaMatrices;
  const char **list           = Gamma::GammaMatrixNames;
  int T = Grid.GlobalDimensions()[3];
  int X = Grid.GlobalDimensions()[0];
  int Y = Grid.GlobalDimensions()[1];
  int Z = Grid.GlobalDimensions()[2];

  std::cout << "lattice dim = " << X << "x"<< Y <<"x"<< Z <<"x"<< T << std::endl;


  for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){//loop over saved cnfg
    Checkpoint.CheckpointRestore(cnfg,Umu, sRNG, pRNG);

    //Gauge Transformation 
    /*
    LatticeColourMatrix Omega(&Grid);
    LatticeColourMatrix ShiftedOmega(&Grid);
    LatticeGaugeField Up(&Grid); Up=zero;
    LatticeColourMatrix U(&Grid); U=zero;
    LatticeColourMatrix Up_mu(&Grid); Up_mu=zero;
    SU<Nc>::LieRandomize(pRNG, Omega, 1.0);
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
    */
    //end transformation



 //Source Position
 std::vector<int> position(4);
 for (int i=0;i<3;i++) position[i]=0; position[3]=0;
 SpinColourMatrix one;
 one=1.0;
 typename TwoIndexRep<Nc, AntiSymmetric>::SpinColourMatrix Hirep_one;
 Hirep_one=1.0;    



 //Dirac Op Fund Rep    
 WilsonFermionR Dw_ud(Umu,Grid,RBGrid,m_ud);
 //Fundamental Point Source
 LatticePropagator src(&Grid); src=zero;
 LatticePropagator Prop_ud(&Grid); Prop_ud=zero;
 LatticeFermion src_vec(&Grid); src_vec=zero;
 LatticeFermion Mdagsrc_vec_ud(&Grid); Mdagsrc_vec_ud=zero;
 
 pokeSite(one,src,position);

 std::cout<< "Fund POINT SOURCE in "  << position  <<std::endl;
 for (unsigned int s = 0; s < Ns; s++){
   for (unsigned int c = 0; c < Nc; c++){
     PropToFerm<LatticePropagator,LatticeFermion,qcd>(src_vec,src,s,c);
     std::cout<< "s = " << s << "c = " << c << std::endl;
     ConjugateGradient<LatticeFermion> CG(1.0e-13, 1000000);
     LatticeFermion Psi_ud(&Grid);
     MdagMLinearOperator<WilsonFermionR, LatticeFermion> MdagM_ud(Dw_ud); //Fund-Fund prop (light-light)
     Psi_ud = zero;
     
     MdagM_ud.AdjOp(src_vec,Mdagsrc_vec_ud);
     CG(MdagM_ud, Mdagsrc_vec_ud, Psi_ud); 

     FermToProp<LatticePropagator,LatticeFermion,qcd>(Prop_ud,Psi_ud,s,c);
   }
 }  

 //Dirac Op Hirep
 HiRep.update_representation(Umu);
 typename TwoIndexRep<Nc, AntiSymmetric>::LatticeField Umu_Hirep = HiRep.U;
 HirepFermionOp Dw_s(Umu_Hirep,Grid,RBGrid,m_s);
 //Hirep Point Source
 typename TwoIndexRep<Nc, AntiSymmetric>::LatticePropagator src_Hirep(&Grid); 
 typename TwoIndexRep<Nc, AntiSymmetric>::LatticePropagator Prop_s(&Grid);    
 HirepFermionField src_vec_Hirep(&Grid); src_vec_Hirep=zero;
 HirepFermionField Mdagsrc_vec_Hirep(&Grid); Mdagsrc_vec_Hirep=zero;

 pokeSite(Hirep_one,src_Hirep,position);

 std::cout<< "Hirep POINT SOURCE in "  << position  <<std::endl;
 for (unsigned int s = 0; s < Ns; s++){
   for (unsigned int c = 0; c < Nc*(Nc-1)/2; c++){
     PropToFerm< TwoIndexRep<Nc, AntiSymmetric>::LatticePropagator, HirepFermionField,  TwoIndexRep<Nc, AntiSymmetric>::parms>(src_vec_Hirep,src_Hirep,s,c);
     std::cout<< "s = " << s << "[i,j] = " << c << std::endl;
     ConjugateGradient<HirepFermionField> Hirep_CG(1.0e-13, 1000000);
     HirepFermionField Psi_Hirep(&Grid);
     MdagMLinearOperator<HirepFermionOp,HirepFermionField> MdagM_Hirep(Dw_s); //Hirep-Hirep prop (heavy-heavy)  
     Psi_Hirep = zero;

     MdagM_Hirep.AdjOp(src_vec_Hirep,Mdagsrc_vec_Hirep);
     Hirep_CG(MdagM_Hirep, Mdagsrc_vec_Hirep, Psi_Hirep);

     FermToProp< TwoIndexRep<Nc, AntiSymmetric>::LatticePropagator, HirepFermionField, TwoIndexRep<Nc, AntiSymmetric>::parms>(Prop_s,Psi_Hirep,s,c);
   }
 }



 LatticeComplex Ctrtr(&Grid); Ctrtr=zero;
 LatticeComplex Ctr(&Grid); Ctr=zero;
 LatticeComplex C(&Grid); C=zero;
 std::vector <TComplex> Ct;


 SpinMatrix Ga,Gb,P;
 SpinMatrix Id; Id=zero;
 RealD unit = 1.0;
 for (int i=0;i<Ns;i++) Id()(i,i)()=1;
 
 // J^P=(1/2)^+
 Ga = Id;
 Gb = makeGammaProd(10) * makeGammaProd(15); // Gb = Gamma(g[4])*(Gamma(g[2])*Gamma(g[5])); = C*g5
 P = 0.5 * ( Id + makeGammaProd(8) ); 
 
 // std::cout<< Gb << "--" << (Ga*P*Ga) << std::endl;
 LatticeSpinMatrix Prop_ud1(&Grid); Prop_ud1 = zero;
 LatticeSpinMatrix Prop_ud2(&Grid); Prop_ud2 = zero;
 LatticeSpinMatrix Prop_s1(&Grid);  Prop_s1  = zero;

#define tr

#ifdef tr
 for (int i=0;i<12;i++){
   for (int j=0;j<12;j++){
     // std::cout << " " << i << " " << j << "so far so good" << std::endl;

     //positive colours permutations
     Prop_s1=peekColour(Prop_s,A[ep[i][1]][ep[i][2]],A[ep[j][1]][ep[j][2]]);
     Prop_ud1=peekColour(Prop_ud,ep[i][0],ep[j][0]);
     Prop_ud2=peekColour(Prop_ud,ep[i][3],ep[j][3]);

     Ctrtr += Sign(ep[i][2],ep[i][1])*Sign(ep[j][2],ep[j][1]) * trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 ); 
     Ctr += Sign(ep[i][2],ep[i][1])*Sign(ep[j][2],ep[j][1]) * trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );  
     //     Ctrtr +=  trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     //     Ctr +=  trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );

     Prop_s1  = zero;
     Prop_ud1 = zero;
     Prop_ud2 = zero;

     Prop_s1=peekColour(Prop_s,A[em[i][1]][em[i][2]],A[em[j][1]][em[j][2]]);
     Prop_ud1=peekColour(Prop_ud,em[i][0],em[j][0]);
     Prop_ud2=peekColour(Prop_ud,em[i][3],em[j][3]);
     Ctrtr += Sign(em[i][2],em[i][1])*Sign(em[j][2],em[j][1]) * trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     Ctr += Sign(em[i][2],em[i][1])*Sign(em[j][2],em[j][1]) * trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     //     Ctrtr +=  trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     //     Ctr +=  trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );

     Prop_s1  = zero;
     Prop_ud1 = zero;
     Prop_ud2 = zero;

     //negative colours permutations
     Prop_s1=peekColour(Prop_s,A[ep[i][1]][ep[i][2]],A[em[j][1]][em[j][2]]);
     Prop_ud1=peekColour(Prop_ud,ep[i][0],em[j][0]);
     Prop_ud2=peekColour(Prop_ud,ep[i][3],em[j][3]);
     Ctrtr -= Sign(ep[i][2],ep[i][1])*Sign(em[j][2],em[j][1]) * trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     Ctr -= Sign(ep[i][2],ep[i][1])*Sign(em[j][2],em[j][1]) * trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     //     Ctrtr -=  trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     //     Ctr -=  trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );

     Prop_s1  = zero;
     Prop_ud1 = zero;
     Prop_ud2 = zero;

     Prop_s1=peekColour(Prop_s,A[em[i][1]][em[i][2]],A[ep[j][1]][ep[j][2]]);
     Prop_ud1=peekColour(Prop_ud,em[i][0],ep[j][0]);
     Prop_ud2=peekColour(Prop_ud,em[i][3],ep[j][3]);
     Ctrtr -= Sign(em[i][2],em[i][1])*Sign(ep[j][2],ep[j][1]) * trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     Ctr -= Sign(em[i][2],em[i][1])*Sign(ep[j][2],ep[j][1]) * trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     //     Ctrtr -=  trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     //     Ctr -=  trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );

     Prop_s1  = zero;
     Prop_ud1 = zero;
     Prop_ud2 = zero;

   }
 }
 C = (Ctrtr + Ctr);                                                                                                                                                            
 sliceSum(C,Ct,3);
#endif

 #ifdef notr
 //check all spin indices (no traces)
 Complex GaA;
 Complex GaB;
 Complex GbA;
 Complex GbB;
 Complex Ppm;
 LatticeComplex Cc(&Grid);   Cc=zero;
 LatticeComplex Ds(&Grid);   Ds=zero;
 LatticeComplex DudA(&Grid); DudA=zero;
 LatticeComplex DudB(&Grid); DudB=zero;
 LatticeComplex DudC(&Grid); DudC=zero;
 LatticeComplex DudD(&Grid); DudD=zero;

 for (int ii=0;ii<3;ii++){
   for (int jj=0;jj<3;jj++){ 

     //peekColour
     Prop_s1=peekColour(Prop_s,ep[ii][1],ep[jj][1]);
     Prop_ud1=peekColour(Prop_ud,ep[ii][0],ep[jj][0]);  
     Prop_ud2=peekColour(Prop_ud,ep[ii][2],ep[jj][2]);
     for (int i=0;i<Ns;i++){
       for (int j=0;j<Ns;j++){
	 for (int k=0;k<Ns;k++){
	   for (int l=0;l<Ns;l++){
	     for (int m=0;m<Ns;m++){
	       for (int n=0;n<Ns;n++){
		 for (int o=0;o<Ns;o++){
		   for (int p=0;p<Ns;p++){
		     //peekSpin
		     GaA=Ga()(k,l)();
		     GaB=Ga()(m,n)();
		     Ppm=P()(l,m)();//lm
		     GbA=Gb()(i,j)();
		     GbB=Gb()(o,p)();
		     Ds=peekSpin(Prop_s1,p,j);
		     DudA=peekSpin(Prop_ud1,n,k);
		     DudB=peekSpin(Prop_ud2,o,i);
		     DudC=peekSpin(Prop_ud1,n,i);
		     DudD=peekSpin(Prop_ud2,o,k);
		     Cc += GaA*Ppm*GaB*GbA*GbB* Ds*(DudA*DudB + DudC*DudD); 
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }  
     C +=Cc;
     Cc=zero;
     //peekColour                                                                                                                   
     Prop_s1=peekColour(Prop_s,em[ii][1],em[jj][1]);
     Prop_ud1=peekColour(Prop_ud,em[ii][0],em[jj][0]);
     Prop_ud2=peekColour(Prop_ud,em[ii][2],em[jj][2]);
     for (int i=0;i<Ns;i++){
       for (int j=0;j<Ns;j++){
         for (int k=0;k<Ns;k++){
           for (int l=0;l<Ns;l++){
             for (int m=0;m<Ns;m++){
               for (int n=0;n<Ns;n++){
                 for (int o=0;o<Ns;o++){
                   for (int p=0;p<Ns;p++){
                     GaA=Ga()(k,l)();
                     GaB=Ga()(m,n)();
                     Ppm=P()(l,m)();
                     GbA=Gb()(i,j)();
                     GbB=Gb()(o,p)();
                     Ds=peekSpin(Prop_s1,p,j);
                     DudA=peekSpin(Prop_ud1,n,k);
                     DudB=peekSpin(Prop_ud2,o,i);
                     DudC=peekSpin(Prop_ud1,n,i);
                     DudD=peekSpin(Prop_ud2,o,k);
                     Cc += GaA*Ppm*GaB*GbA*GbB* Ds*(DudA*DudB + DudC*DudD);
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }
     C +=Cc;
     Cc=zero;
     
     //peekColour                                                                                                                                                              
     Prop_s1=peekColour(Prop_s,em[ii][1],ep[jj][1]);
     Prop_ud1=peekColour(Prop_ud,em[ii][0],ep[jj][0]);
     Prop_ud2=peekColour(Prop_ud,em[ii][2],ep[jj][2]);
     for (int i=0;i<Ns;i++){
       for (int j=0;j<Ns;j++){
	 for (int k=0;k<Ns;k++){
	   for (int l=0;l<Ns;l++){
	     for (int m=0;m<Ns;m++){
	       for (int n=0;n<Ns;n++){
		 for (int o=0;o<Ns;o++){
		   for (int p=0;p<Ns;p++){
		     GaA=Ga()(k,l)();
		     GaB=Ga()(m,n)();
		     Ppm=P()(l,m)();
		     GbA=Gb()(i,j)();
		     GbB=Gb()(o,p)();
		     Ds=peekSpin(Prop_s1,p,j);
		     DudA=peekSpin(Prop_ud1,n,k);
		     DudB=peekSpin(Prop_ud2,o,i);
		     DudC=peekSpin(Prop_ud1,n,i);
		     DudD=peekSpin(Prop_ud2,o,k);
		     Cc += GaA*Ppm*GaB*GbA*GbB* Ds*(DudA*DudB + DudC*DudD);
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }
     C -=Cc;
     Cc=zero;
     
     Prop_s1=peekColour(Prop_s,ep[ii][1],em[jj][1]);
     Prop_ud1=peekColour(Prop_ud,ep[ii][0],em[jj][0]);
     Prop_ud2=peekColour(Prop_ud,ep[ii][2],em[jj][2]);
     for (int i=0;i<Ns;i++){
       for (int j=0;j<Ns;j++){
	 for (int k=0;k<Ns;k++){
	   for (int l=0;l<Ns;l++){
	     for (int m=0;m<Ns;m++){
	       for (int n=0;n<Ns;n++){
		 for (int o=0;o<Ns;o++){
		   for (int p=0;p<Ns;p++){
		     GaA=Ga()(k,l)();
		     GaB=Ga()(m,n)();
		     Ppm=P()(l,m)();
		     GbA=Gb()(i,j)();
		     GbB=Gb()(o,p)();
		     Ds=peekSpin(Prop_s1,p,j);
		     DudA=peekSpin(Prop_ud1,n,k);
		     DudB=peekSpin(Prop_ud2,o,i);
		     DudC=peekSpin(Prop_ud1,n,i);
		     DudD=peekSpin(Prop_ud2,o,k);
		     Cc += GaA*Ppm*GaB*GbA*GbB* Ds*(DudA*DudB + DudC*DudD);
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }
     C -=Cc;
     Cc=zero;
   }
 } 
 sliceSum(C,Ct,3);
#endif 


 std::cout << "t     <CC>" << std::endl;
 for (int t=0;t<T;t++){
   std::cout << " " << t << "     " << TensorRemove(Ct[t]).real() << "     " << TensorRemove(Ct[t]).imag() << " "<< std::endl;
 }
 
  } 
  Grid_finalize();
}

