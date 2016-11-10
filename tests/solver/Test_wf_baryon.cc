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
  //SU3::TepidConfiguration(pRNG, Umu);  

  RealD m_ud=0.1;
  RealD m_s=0.1;


#define POINT_SOURCES

  
  NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(std::string("ckpoint_SU3new_lat"),
                                                  std::string("ckpoint_SU3new_rng"), 1);

  int CNFGSTART=29;
  int CNFGEND=29;
  int CNFGSTEP=1;

  
  for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){   //loop over saved cnfg                                                                                                
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

WilsonFermionR Dw_ud(Umu,Grid,RBGrid,m_ud);
WilsonFermionR Dw_s(Umu,Grid,RBGrid,m_s);    

#ifdef POINT_SOURCES
 LatticePropagator src(&Grid); src=zero;
 LatticePropagator Prop_ud(&Grid); Prop_ud=zero;
 LatticePropagator Prop_s(&Grid); Prop_s=zero;
 LatticeFermion src_vec(&Grid); src_vec=zero;
 LatticeFermion Mdagsrc_vec_ud(&Grid); Mdagsrc_vec_ud=zero;
 LatticeFermion Mdagsrc_vec_s(&Grid); Mdagsrc_vec_s=zero;

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
     ConjugateGradient<LatticeFermion> CG(1.0e-12, 1000000);
     LatticeFermion Psi_ud(&Grid);
     LatticeFermion Psi_s(&Grid);
     MdagMLinearOperator<WilsonFermionR, LatticeFermion> MdagM_ud(Dw_ud); //light prop
     MdagMLinearOperator<WilsonFermionR, LatticeFermion> MdagM_s(Dw_s); //heavy prop
     Psi_ud = zero;
     Psi_s = zero;
     
     MdagM_ud.AdjOp(src_vec,Mdagsrc_vec_ud);
     MdagM_s.AdjOp(src_vec,Mdagsrc_vec_s);

     CG(MdagM_ud, Mdagsrc_vec_ud, Psi_ud); 
     CG(MdagM_s, Mdagsrc_vec_s, Psi_s); 

     FermToProp<LatticePropagator,LatticeFermion,qcd>(Prop_ud,Psi_ud,s,c);
     FermToProp<LatticePropagator,LatticeFermion,qcd>(Prop_s,Psi_s,s,c);
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
       ConjugateGradient<LatticeFermion> CG(1.0e-6, 10);
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

 int Nperm=3; //Nc!/2

 int ep[3][Nc];
 int em[3][Nc];
 int contp=0,contm=0;
 int S1,S2,S3,S4,S5,S6;

 if (Nc==3){
  // qqq  with q in fund SU(3)
   // \epsilon_{ijk}
   for (int i=0;i<Nc;i++){
     for (int j=0;j<Nc;j++){
       for (int k=0;k<Nc;k++){
	 S1=Sign(i,j);
	 S2=Sign(j,k);
	 S3=Sign(i,k);
	 if ((S1*S2*S3)>0){
	   ep[contp][0]=i;
	   ep[contp][1]=j;
	   ep[contp][2]=k;
	   contp++;
	 } else if ((S1*S2*S3)<0){
	   em[contm][0]=i;
	   em[contm][1]=j;
	   em[contm][2]=k;
	   contm++;
	 }
       }
     }
   }
  } else if (Nc==4){
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
 }
 
 /*
 for (int i=0;i<Nperm;i++){
   for (int j=0;j<Nc;j++){
     std::cout<< "ep["<< i << "]["<< j << "]=" << ep[i][j] <<std::endl;
   }
 }
 std::cout<<"-----"<<std::endl;
 for (int i=0;i<Nperm;i++){
   for (int j=0;j<Nc;j++){
     std::cout<< "em["<< i << "]["<< j << "]=" << em[i][j] <<std::endl;
   }
 }
 */
 

 //Mapping (i,j) -> A , with i,j=1,...,4 and A=1,...,6
 int A[Nc][Nc];
 int cont=0; 
 for (int a=1;a<Nc;a++){
   for (int b=0;b<a;b++){
     A[a][b]=cont;
     cont++;
   }
 }


 for (int a=0;a<Nc;a++){
   for (int b=0;b<Nc;b++){
     std::cout << "A[" << a << "]["<<b<<"]= "<<A[a][b]<<std::endl;
   }
 }

 
 //Correlators
 /* Baryon (llh):                                                                                                                                                                   
  e_{abc} e_{a'b'c'}  ( Tr[(Ga Ppm Ga) S^{light}_{ca'}]Tr[(Gb S^{heavy}_{bc'} Gb^T)^T S^{light}_{ab'}] -
  + Tr[(Gb S^{heavy}_{bc'} Gb^T)^T S^{light}_{aa'} (Ga Ppm Ga) S^{light}_{bb'}] )
                                                                            
  Ppm = (Id±g4)/2; (parity projector)                                                                                                                                        
  C = g0g2;        (charge conjugation)

  J^P = (1/2)^+ -> (Ga,Gb) = (Id,Cg5),(g5,C),(Id,ig4Cg5)                                                                                                            
  J^P = (3/2)^+ -> (Ga,Gb) = (Id,Cgj), with any j;                                                                                                
 */

 const Gamma::GammaMatrix *g = Gamma::GammaMatrices;
 const char **list           = Gamma::GammaMatrixNames;


 
 int T = Grid.GlobalDimensions()[3];
 int X = Grid.GlobalDimensions()[0];
 int Y = Grid.GlobalDimensions()[1];
 int Z = Grid.GlobalDimensions()[2];
 
 std::cout << "lattice dim = " << X << "x"<< Y <<"x"<< Z <<"x"<< T << std::endl;

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

 // Gb = Gamma(g[4])*(Gamma(g[2])*Gamma(g[5])); = C*g5
  Gb = makeGammaProd(10) * makeGammaProd(15);
 
 P = 0.5 * ( Id + makeGammaProd(8) ); 
 
 // std::cout<< Gb << "--" << (Ga*P*Ga) << std::endl;
 LatticeSpinMatrix Prop_ud1(&Grid); Prop_ud1 = zero;
 LatticeSpinMatrix Prop_ud2(&Grid); Prop_ud2 = zero;
 LatticeSpinMatrix Prop_s1(&Grid);  Prop_s1  = zero;

#define tr

#ifdef tr
 for (int i=0;i<3;i++){
   for (int j=0;j<3;j++){

     //positive colours permutations
     Prop_s1=peekColour(Prop_s,ep[i][1],ep[j][1]);
     Prop_ud1=peekColour(Prop_ud,ep[i][0],ep[j][0]);
     Prop_ud2=peekColour(Prop_ud,ep[i][2],ep[j][2]);
     Ctrtr += trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 ); 
     Ctr += trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );  

     Prop_s1  = zero;
     Prop_ud1 = zero;
     Prop_ud2 = zero;

     Prop_s1=peekColour(Prop_s,em[i][1],em[j][1]);
     Prop_ud1=peekColour(Prop_ud,em[i][0],em[j][0]);
     Prop_ud2=peekColour(Prop_ud,em[i][2],em[j][2]);
     Ctrtr += trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     Ctr += trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );

     Prop_s1  = zero;
     Prop_ud1 = zero;
     Prop_ud2 = zero;

     //negative colours permutations
     Prop_s1=peekColour(Prop_s,ep[i][1],em[j][1]);
     Prop_ud1=peekColour(Prop_ud,ep[i][0],em[j][0]);
     Prop_ud2=peekColour(Prop_ud,ep[i][2],em[j][2]);
     Ctrtr -= trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     Ctr -= trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );

     Prop_s1  = zero;
     Prop_ud1 = zero;
     Prop_ud2 = zero;

     Prop_s1=peekColour(Prop_s,em[i][1],ep[j][1]);
     Prop_ud1=peekColour(Prop_ud,em[i][0],ep[j][0]);
     Prop_ud2=peekColour(Prop_ud,em[i][2],ep[j][2]);
     Ctrtr -= trace( (Ga*P*Ga) * Prop_ud1 ) * trace( transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );
     Ctr -= trace( (Ga*P*Ga) * Prop_ud1 * transpose(Gb*Prop_s1*transpose(Gb)) * Prop_ud2 );

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

