#include <dolfin.h>
#include <mshr.h>
#include "./BeautifiedUFL/Wrapper.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

using namespace dolfin;
using namespace std;

/**
 * @brief      Boundary on which environmental condition for liquid pressure is applied
 */
class TopSatLiquid: public SubDomain
{
public:
  /**
   * @brief      Test if position is correct and if seepage occurs
   *
   * @param[in]  x            { parameter_description }
   * @param[in]  on_boundary  On boundary
   *
   * @return     x is on a layer connected to air && liquid outflow is positive && Saturation occurs
   */
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    if((x[1]-4)>-5*DOLFIN_EPS ){
      dolfin::Array<double> values_pl((x).size());
      dolfin::Array<double> values_vl((x).size());
      pl.eval(values_pl, x);
      vl.eval(values_vl, x);
      return (on_boundary and values_vl[1] > eps_vl and  (values_pl[0])-p_atm >eps_pc);
    }else{
      return(false);
    }// top & saturé & débit sortant
  }
  /***
   * @brief      Sets the functions.
   *
   * @param[in]  pl_   Liquid pressure
   * @param[in]  pg_   Gaz pressure
   * @param[in]  vl_   Liquid velocity vector
   */
  void set_functions(dolfin::Function pl_, dolfin::Function pg_, dolfin::Function vl_)
  {
    pl = pl_;
    pg = pg_;
    vl = vl_;
    pl.set_allow_extrapolation(true);
    pg.set_allow_extrapolation(true);
  }
  void set_outflow_minimum(double eps_vl_){
    eps_vl = eps_vl_;
  }
  void set_pressure_minimum(double eps_pc_){
    eps_pc = eps_pc_;
  }
  void set_p_atm(double p_atm_){
    p_atm = p_atm_;
  }
private:
  dolfin::Function pl;
  dolfin::Function pg;
  dolfin::Function vl;
  double eps_vl;
  double eps_pc;
  double p_atm;
};


/**
 * @brief      Top Layer y = ymax(x)
 */
class Top: public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return ( on_boundary and (abs(x[1]-4)<2*DOLFIN_EPS) );
  }
};

/**
 * @brief      Bottom layer y = 0
 */
class Bottom: public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (on_boundary and x[1]<DOLFIN_EPS);
  }
};

/**
 * @brief      Bottom layer y = 0
 */
class BottomCenter: public SubDomain
{
public:
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (on_boundary and x[1]<DOLFIN_EPS and abs(2-x[0]) <radius);
  }
  void set_radius(double radius_){
    radius = radius_;
  }
private:
  double radius;
};

/**
 * @brief      Left layer x = 0
 */
class Left: public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (on_boundary and x[0]<DOLFIN_EPS);
  }
};

/**
 * @brief      Right layer x = xmax
 */
class Right: public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (on_boundary and x[0]>4-DOLFIN_EPS);
  }
};


/**
 * @brief      Adimensional analytical initial liquid pressure profile
 */
class LiquidPressureProfile : public Expression
{
public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] =p_atm+4.01-x[1];   // DRAINAGE
    // values[0] =12.5-x[1];   // IMBIBITION
 }
   void set_p_atm(double p_atm_){
    p_atm = p_atm_;
  }
private:
  double p_atm;
};

/**
 * @brief      Adimensional analytical initial gaz pressure profile
 */
class GazPressureProfile : public Expression
{
public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
    // if(x[1]<4.2){values[0] = 0.9999*((p_atm+4.5-x[1]));}else{values[0] = p_atm;} //DRAINAGE
    if(x[1]<4.01){values[0] = 0.9*((p_atm+4.01-x[1]));}else{values[0] = p_atm;} //IMBIBITION
  }
  void set_p_atm(double p_atm_){
    p_atm = p_atm_;
  }
private:
  double p_atm;
};

class LogGenerator
{
public:
  LogGenerator(std::string filename_, bool activated_)
  {
    filename = std::string(filename_);
    std::ofstream ofs;
    // std::cout<<filename<<std::endl;
    ofs.open(filename,std::ofstream::trunc);
    ofs << "LOGFILE\n" <<std::endl;
    ofs.close();
    activated = activated_;
  }
  void write(std::string text){
    ofs.open(filename,std::ofstream::app);
    ofs << text <<std::endl;
    ofs.close();
    if(activated){std::cout<< text <<std::endl;}
  }
  void writeandshout(std::string text){
    ofs.open(filename,std::ofstream::app);
    ofs << text <<std::endl;
    ofs.close();
    std::cout<< text <<std::endl;
  }
  
  void set_filename(std::string foo){
    filename = std::string(foo);
  }
private:
  std::string filename;
  std::ofstream ofs;
  bool activated;
};

/**
 * @brief      main program
 * 
 * Organized in multiple parts:
 * 
 * Save
 * Solver
 * Mesh and Forms
 * Constants
 * Initial profiles
 * Boundary conditions
 * 
 * Nested Algorithms: Euler, Alterned Scheme for BCs, Newton-Raphson
 */
int main(int argc, char *argv[])
{
  #ifdef HAS_PETSC

  // #### Save ####
  // Folder name for results
  // double modif = 1/1; 7/4 7/5 6/5 1/1 7/8 6/8 5/8

  static const std::string results_folder_name = "./Results/rate"+std::string(argv[1])+"-"+std::string(argv[2]);
  if(mkdir((results_folder_name).c_str(),0777)==-1){std::cout<<"WARNING: Results folder may exist\n"<<std::endl;}
  static const std::string log_file_name = results_folder_name +"/log.txt" ;
  std::ifstream  src("../main.cpp", std::ios::binary);
  std::ofstream  dst( results_folder_name +"/script.cpp",   std::ios::binary);
  dst << src.rdbuf();  
  static const std::string mesh_file_name ="../mesh/Sink3.xml" ;

  // Enable or Disable Log prompt
  static const bool ENABLE_LOG = false;
  LogGenerator log(log_file_name,ENABLE_LOG);

  
  
  log.writeandshout("Input: " + std::string(argv[1])+" - "+ std::string(argv[2])+"\n" );
  const int modif1 = ((int)std::stol(argv[1], nullptr, 0));
  const int modif2 = ((int)std::stol(argv[2], nullptr, 0));

  log.writeandshout("\n--- SAVE FILES ---\n");


  // Save files definition for each time step
  static const bool SAVE_COMMON = false;
  File plfile(results_folder_name +"/00_Pl.pvd");
  File pgfile(results_folder_name +"/00_Pg.pvd");
  File ufile(results_folder_name +"/00_U.pvd");

  //Enable or Disable saving Newton-Raphson transient data
  static const bool SAVE_NEWTON = false;
  File plfileK(results_folder_name +"/00_Pl_K.pvd");
  File pgfileK(results_folder_name +"/00_Pg_K.pvd");
  File ufileK(results_folder_name +"/00_U_K.pvd");
  File plfileDelta(results_folder_name +"/00_deltaPl.pvd");
  File pgfileDelta(results_folder_name +"/00_deltaPg.pvd");
  File ufileDelta(results_folder_name +"/00_deltaU.pvd");
  File pprojfileK(results_folder_name +"/00_Psatproj_K.pvd");
  
  // Enable or Disable saving stress tensor and velocity vector at each time step
  static const bool SAVE_DUAL = false;
  File stress(results_folder_name +"/00_sigma.pvd");
  File fluxL(results_folder_name +"/00_VL.pvd");
  File fluxG(results_folder_name +"/00_Vg.pvd");

  // Save file to follow subdomains definition
  File subdomainLiq_file(results_folder_name +"/subdomainsLiq.pvd");
  File subdomainGaz_file(results_folder_name +"/subdomainsGaz.pvd");

  static const bool SAVE_POINT = true;
  static const std::string eps_vpoint_file_name = results_folder_name +"/eps_vpoint.txt" ;
  LogGenerator eps_vpointfile(eps_vpoint_file_name,false);

  static const std::string plpoint_file_name = results_folder_name +"/plpoint.txt" ;
  LogGenerator plpointfile(plpoint_file_name,false);

  static const std::string pgpoint_file_name = results_folder_name +"/pgpoint.txt" ;
  LogGenerator pgpointfile(pgpoint_file_name,false);

  static const std::string uxpoint_file_name = results_folder_name +"/uxpoint.txt" ;
  LogGenerator uxpointfile(uxpoint_file_name,false);

  static const std::string uypoint_file_name = results_folder_name +"/uypoint.txt" ;
  LogGenerator uypointfile(uypoint_file_name,false);

   double p0[3] = {2.,0.,0. }; eps_vpointfile.write("("+to_string(p0[0]) + ";"+to_string(p0[1]) + ";"+ to_string(p0[2]) + ")" );
  pgpointfile.write("("+to_string(p0[0]) + ";"+to_string(p0[1]) + ";"+ to_string(p0[2]) + ")" );
  plpointfile.write("("+to_string(p0[0]) + ";"+to_string(p0[1]) + ";"+ to_string(p0[2]) + ")" );
  uxpointfile.write("("+to_string(p0[0]) + ";"+to_string(p0[1]) + ";"+ to_string(p0[2]) + ")" );
  uypointfile.write("("+to_string(p0[0]) + ";"+to_string(p0[1]) + ";"+ to_string(p0[2]) + ")" );
  const dolfin::Array<double> point0(sizeof(p0)/sizeof(double),p0);
  double p1[3] = {2.,1.,0. }; eps_vpointfile.write("("+to_string(p1[0]) + ";"+to_string(p1[1]) + ";"+ to_string(p1[2]) + ")" );
  pgpointfile.write("("+to_string(p1[0]) + ";"+to_string(p1[1]) + ";"+ to_string(p1[2]) + ")" );
  plpointfile.write("("+to_string(p1[0]) + ";"+to_string(p1[1]) + ";"+ to_string(p1[2]) + ")" );
  uxpointfile.write("("+to_string(p1[0]) + ";"+to_string(p1[1]) + ";"+ to_string(p1[2]) + ")" );
  uypointfile.write("("+to_string(p1[0]) + ";"+to_string(p1[1]) + ";"+ to_string(p1[2]) + ")" );
  const dolfin::Array<double> point1(sizeof(p1)/sizeof(double),p1);
  double p2[3] = {2.,2.,0. }; eps_vpointfile.write("("+to_string(p2[0]) + ";"+to_string(p2[1]) + ";"+ to_string(p2[2]) + ")" );
  pgpointfile.write("("+to_string(p2[0]) + ";"+to_string(p2[1]) + ";"+ to_string(p2[2]) + ")" );
  plpointfile.write("("+to_string(p2[0]) + ";"+to_string(p2[1]) + ";"+ to_string(p2[2]) + ")" );
  uxpointfile.write("("+to_string(p2[0]) + ";"+to_string(p2[1]) + ";"+ to_string(p2[2]) + ")" );
  uypointfile.write("("+to_string(p2[0]) + ";"+to_string(p2[1]) + ";"+ to_string(p2[2]) + ")" );
  const dolfin::Array<double> point2(sizeof(p1)/sizeof(double),p2);
  double p3[3] = {2.,3.,0. }; eps_vpointfile.write("("+to_string(p3[0]) + ";"+to_string(p3[1]) + ";"+ to_string(p3[2]) + ")" );
  pgpointfile.write("("+to_string(p3[0]) + ";"+to_string(p3[1]) + ";"+ to_string(p3[2]) + ")" );
  plpointfile.write("("+to_string(p3[0]) + ";"+to_string(p3[1]) + ";"+ to_string(p3[2]) + ")" );
  uxpointfile.write("("+to_string(p3[0]) + ";"+to_string(p3[1]) + ";"+ to_string(p3[2]) + ")" );
  uypointfile.write("("+to_string(p3[0]) + ";"+to_string(p3[1]) + ";"+ to_string(p3[2]) + ")" );
  const dolfin::Array<double> point3(sizeof(p1)/sizeof(double),p3);
  double p4[3] = {2.,4.,0. }; eps_vpointfile.write("("+to_string(p4[0]) + ";"+to_string(p4[1]) + ";"+ to_string(p4[2]) + ")" );
  pgpointfile.write("("+to_string(p4[0]) + ";"+to_string(p4[1]) + ";"+ to_string(p4[2]) + ")" );
  plpointfile.write("("+to_string(p4[0]) + ";"+to_string(p4[1]) + ";"+ to_string(p4[2]) + ")" );
  uxpointfile.write("("+to_string(p4[0]) + ";"+to_string(p4[1]) + ";"+ to_string(p4[2]) + ")" );
  uypointfile.write("("+to_string(p4[0]) + ";"+to_string(p4[1]) + ";"+ to_string(p4[2]) + ")" );
  const dolfin::Array<double> point4(sizeof(p1)/sizeof(double),p4);

  static const bool SAVE_PROFILINIT = true;
  static const std::string plinit_file_name = results_folder_name +"/plinit.txt" ;
  LogGenerator plinitfile(plinit_file_name,false);
  static const std::string pginit_file_name = results_folder_name +"/pginit.txt" ;
  LogGenerator pginitfile(pginit_file_name,false);

  static const bool SAVE_PROFILMID = true;
  static const std::string plmid_file_name = results_folder_name +"/plmid.txt" ;
  LogGenerator plmidfile(plmid_file_name,false);
  static const std::string pgmid_file_name = results_folder_name +"/pgmid.txt" ;
  LogGenerator pgmidfile(pgmid_file_name,false);

  static const bool SAVE_PROFILFIN = true;
  static const std::string plfin_file_name = results_folder_name +"/plfin.txt" ;
  LogGenerator plfinfile(plfin_file_name,false);
  static const std::string pgfin_file_name = results_folder_name +"/pgfin.txt" ;
  LogGenerator pgfinfile(pgfin_file_name,false);

  // #### Solver ####
  log.writeandshout("\n--- SOLVER ---\n");
  // Set PETSc solve type and preconditioner
  // note : https://petsc4py.readthedocs.io/en/stable/manual/ksp/
  PETScOptions::set("ksp_type", "preonly");
  PETScOptions::set("pc_type", "lu");
  PETScOptions::set("pc_factor_mat_solver_type", "mumps");
  PETScKrylovSolver solver;
  solver.set_from_options();

  // #### Mesh and Forms ####
  log.writeandshout("\n--- MESH and FORMS ---\n");
  auto mesh = std::make_shared<Mesh>(mesh_file_name); 

  // Basic geometry instanciation 
  auto top = std::make_shared<Top>();   
  auto bot = std::make_shared<Bottom>();
  //auto botCenter = std::make_shared<BottomCenter>(); botCenter->set_radius(0.3);
  auto right = std::make_shared<Right>(); 
  auto left = std::make_shared<Left>(); 

  // Common Function Space : (x,y) displacement, liquid pressure, fluid pressure
  auto W  = std::make_shared<Newton2D::FunctionSpace>(mesh);

  // Sum forms and matrices
  auto aSum = std::make_shared<Sum2D::BilinearForm>(W,W);
  auto LSum = std::make_shared<Sum2D::LinearForm>(W);
  auto ASum = std::make_shared<PETScMatrix>();
  auto bSum = std::make_shared<PETScVector>(); 

  // Newton-Raphson forms and matrices
  auto a = std::make_shared<Newton2D::BilinearForm>(W,W);
  auto L = std::make_shared<Newton2D::LinearForm>(W);
  auto A = std::make_shared<PETScMatrix>();
  auto b = std::make_shared<PETScVector>(); 

  // Initial field instantiation forms and matrices
  auto aInit = std::make_shared<FirstGuess2D::BilinearForm>(W,W);
  auto LInit = std::make_shared<FirstGuess2D::LinearForm>(W);
  auto AInit = std::make_shared<PETScMatrix>();
  auto bInit = std::make_shared<PETScVector>(); 

  // Initial field instantiation forms and matrices
  auto aInitN = std::make_shared<InitialNewton2D::BilinearForm>(W,W);
  auto LInitN = std::make_shared<InitialNewton2D::LinearForm>(W);
  auto AInitN = std::make_shared<PETScMatrix>();
  auto bInitN = std::make_shared<PETScVector>();

  // Dual forms and matrices
  auto WDual = std::make_shared<Dual2D::FunctionSpace>(mesh);
  auto aDual = std::make_shared<Dual2D::BilinearForm>(WDual,WDual);
  auto LDual = std::make_shared<Dual2D::LinearForm>(WDual);
  auto ADual = std::make_shared<PETScMatrix>();
  auto bDual = std::make_shared<PETScVector>(); 

  // Surface Pressure correction forms and matrices
  auto Pproj = std::make_shared<SatProj2D::FunctionSpace>(mesh);
  auto aproj = std::make_shared<SatProj2D::BilinearForm>(Pproj,Pproj);
  auto Lproj = std::make_shared<SatProj2D::LinearForm>(Pproj);
  auto Aproj = std::make_shared<PETScMatrix>();
  auto bproj = std::make_shared<PETScVector>(); 

  // Environmental Boundary surfaces forms
  auto formSatSwitch = std::make_shared<SurfInt2D::Form_f1>(mesh);
  auto formUnsatSwitch = std::make_shared<SurfInt2D::Form_f2>(mesh);

  // Environmental Boundary surfaces forms
  auto formQuasiStatLRel = std::make_shared<QuasiStat2D::Form_f1>(mesh);
  auto formQuasiStatGRel = std::make_shared<QuasiStat2D::Form_f2>(mesh);
  auto formQuasiStatLt = std::make_shared<QuasiStat2D::Form_f3>(mesh);
  auto formQuasiStatGt = std::make_shared<QuasiStat2D::Form_f4>(mesh);
  auto formQuasiStatLSpeedrel = std::make_shared<QuasiStat2D::Form_f7>(mesh);
  auto formQuasiStatGSpeedrel = std::make_shared<QuasiStat2D::Form_f8>(mesh);
  auto formQuasiStatLSpeedt = std::make_shared<QuasiStat2D::Form_f11>(mesh);
  auto formQuasiStatGSpeedt = std::make_shared<QuasiStat2D::Form_f12>(mesh);

  // Sum forms and matrices
  auto WVolD  = std::make_shared<VolDef::FunctionSpace>(mesh);
  auto aVolD = std::make_shared<VolDef::BilinearForm>(WVolD,WVolD);
  auto LVolD = std::make_shared<VolDef::LinearForm>(WVolD);
  auto AVolD = std::make_shared<PETScMatrix>();
  auto bVolD = std::make_shared<PETScVector>(); 

  // #### Constants ####
  log.writeandshout("\n--- CONSTANTS ---\n");
  // Solid density
  double rhoSvar = 2650.0;
  // Liquid density
  double rhoLvar = 1000.0;
  // Caracteristic length
  double cvar = 1;
  // Liquid Viscosity
  double Etalvar = 0.001;
  // Gaz viscosity
  double Etagvar = 0.00002;
  // Gravitational acceleration
  double gvar = 9.81; // 9.81 m.s-2

  // A Custom Gaz Constant for density to pressure conversion for an Ideal Gaz: Molar mass/ (R * Temperature)  
  double Cstegvar = 28.895*(1e-3)/(8.314*293.15); 

  // Initial porosity
  double phi0var = 0.35;
  // Solid Bulk Modulus
  double Kvar = 3571428.6;// 10 GPa
  // Solid Shear Modulus
  double Muvar = 1923076.9; // 10 GPa
  // Biot's coefficient
  double b_biotvar = 0.97;  
  double N_biotvar = (Kvar/(1-b_biotvar))/(b_biotvar-phi0var); // 500 GPa 

  // Intrinsic Permeability, based on water viscosity and water flow permeability of sand
  double Kivar = 1e-8*Etalvar; // 10-12 m.s //FIXME
  // Gaz permeability
  double Kigvar = Kivar;
  // Liquid permeability
  double Kilvar = Kivar;

  // van Genuchten's constant for capillary pressure to saturation degree conversion 
  // Capillary modulus
  double Mvar = 5000;
  // van Genuchten's constant for air-water-? mixture
  double mvar = 0.8;
  // entry pressure
  double pevar = 0;//gvar*rhoLvar*cvar;

  // Adimensional base constant
  // adimensional mass, based on liquid density and caracteristical length
  double masse_carac = rhoLvar*cvar*cvar*cvar;
  // adimensional pressure, based on a liquid head pressure for a caracteristical length
  double pression_carac = rhoLvar*gvar*cvar;//rhoLvar*gvar*cvar;
  // adimensional time, based on fully saturated liquid filtration phenomenom
  double temps_carac = Etalvar*cvar*cvar/(Kivar*pression_carac);  


 
  /** TEMPS - TIME -------------------------------------
   *
   * Criteria for low-cut filter come from Vermeer,Verruijt (1981) and is:
   * dt >= temps_carac * (1/6 )*(dx)²/(ki*K/eta)
   * 
   * It is implemented based on liquid flow in a the saturated porous medium
   * Gaz viscosity is 50 times less, gaz and liquid carateristic times differ with the same ratio
   * Relative permeabilities of each fluid slow down flows, there's a margin 
   * allowing lower time step than the one used in liquid saturated simulation 
   */ 
  // Time step
  double dt = temps_carac/100.;//1e-3//1e-2 
  // Total simulaiton duration
  double duration = 100000*dt; //600*dt;// 1.0*temps_carac;
  // Save frequency for Euler method results
  int save_freq = 1;

  // A force norm
  double omegavar = -0*1.0 *1000.0*9.81; // 1tonne sur 1 m2  (10KPa) //(premier nombre : masse en tonne sur 1 m²)
  // a gaz FILTRATION vector norm -- Pay attention to gaz density here
  double Vgvar = -0*(1e-0)*cvar/temps_carac;
  // a liquid velocity vector norm 
  double Vlvar = -0*(1e-0)*cvar/temps_carac; //  cvar/temps_carac = 16*100 dx/dt

  double p_widthvar = rhoLvar*gvar*(3*mesh->hmin());
  log.writeandshout("pH = "+ to_string(rhoLvar*gvar*(3*mesh->hmin())));

  // Adimensionalisation and insertion in forms attributes
  auto h = std::make_shared<Constant>(dt/temps_carac);                a->h = h; L->h = h;aInitN->h = h;LInitN->h = h; //Lproj->h=h;
  auto rhoS = std::make_shared<Constant>(rhoSvar/rhoLvar);            L->rhoS = rhoS;  LInitN->rhoS = rhoS;//a->rhoS = rhoS; aInitN->rhoS = rhoS; // LDual->rhoS = rhoS;
  auto rhoL = std::make_shared<Constant>(rhoLvar/rhoLvar);            L->rhoL = rhoL; a->rhoL = rhoL; LInitN->rhoL = rhoL; aInitN->rhoL = rhoL; LDual->rhoL = rhoL;
  auto T = std::make_shared<Constant>(0.0,omegavar/pression_carac);   L->T = T;
  auto Vg = std::make_shared<Constant>(0.0,Vgvar/(cvar/temps_carac)); //L->Vg = Vg;
  auto Vl = std::make_shared<Constant>(0.0,Vlvar/(cvar/temps_carac)); //L->Vl = Vl;
  auto K = std::make_shared<Constant>(Kvar/pression_carac);           a->K = K; L->K = K;  LDual->K = K; aInitN->K = K;  LInitN->K = K;  
  auto Mu = std::make_shared<Constant>(Muvar/pression_carac);         a->Mu = Mu; L->Mu = Mu; LDual->Mu = Mu;aInitN->Mu = Mu;LInitN->Mu = Mu;  
  auto N_biot = std::make_shared<Constant>(N_biotvar/pression_carac); //a->N = N_biot; L->N = N_biot; // LDual->N = N_biot; 
  auto b_biot = std::make_shared<Constant>(b_biotvar);                //a->b = b_biot; L->b = b_biot; LDual->b = b_biot;;
  auto phi0 = std::make_shared<Constant>(phi0var);                    a->phi0 = phi0; L->phi0 = phi0; aInitN->phi0 = phi0;LInitN->phi0 = phi0; // LDual->phi0 = phi0;
  auto Kil = std::make_shared<Constant>(Kilvar/Kivar);                a->Kil = Kil; L->Kil = Kil; LDual->Kil = Kil;LInitN->Kil = Kil;aInitN->Kil = Kil;
  auto Etal = std::make_shared<Constant>(Etalvar/Etalvar);             a->Etal = Etal; L->Etal = Etal; LDual->Etal = Etal;aInitN->Etal = Etal; LInitN->Etal = Etal;
  auto Kig = std::make_shared<Constant>(Kigvar/Kivar);                a->Kig = Kig; L->Kig = Kig; LDual->Kig = Kig; LInitN->Kig = Kig;aInitN->Kig = Kig;
  auto Etag = std::make_shared<Constant>(Etagvar/Etalvar);             a->Etag = Etag; L->Etag = Etag; LDual->Etag = Etag;aInitN->Etag = Etag; LInitN->Etag = Etag;
  auto Csteg = std::make_shared<Constant>(Cstegvar*pression_carac/rhoLvar); 
                                                                      // a->Csteg = Csteg; L->Csteg = Csteg; LDual->Csteg = Csteg; aInit->Csteg = Csteg; LInit->Csteg = Csteg;
  auto M = std::make_shared<Constant>(Mvar/pression_carac);           a->M = M; L->M = M; LDual->M = M; aInitN->M = M; LInitN->M = M;  
  auto m = std::make_shared<Constant>(mvar);                          a->m = m; L->m = m; LDual->m = m; aInitN->m = m; LInitN->m = m;  
  auto pe = std::make_shared<Constant>(pevar/pression_carac);         a->pe = pe; L->pe = pe; LDual->pe = pe; aInitN->pe = pe; LInitN->pe = pe;  
  auto g = std::make_shared<Constant>(0.0,(-gvar/(pression_carac/(rhoLvar*cvar))));           
                                                                      a->g = g; L->g = g; LDual->g = g; aInitN->g = g; LInitN->g = g;  
  auto p_sat = std::make_shared<Constant>(2*p_widthvar/pression_carac);    L->p_sat = p_sat;LInitN->p_sat=p_sat;aInitN->p_sat=p_sat;a->p_sat=p_sat;     //Lproj->epsilon_g = epsilon_g; // FIXME: TROP PETIT -> BRUIT SUR IMBI
  auto epsilon_l = std::make_shared<Constant>(0);  Lproj->epsilon_l = epsilon_l;
  // Initial function - FIXME: No Use
  auto w0 = std::make_shared<Function>(W);                            L->w0 = w0; LDual->w0=w0;LInitN->w0 = w0;
  auto omega = std::make_shared<Constant>(1); a->omega = omega;L->omega = omega;LInitN->omega=omega;aInitN->omega=omega;


  // Null pressure 
  auto p_null = std::make_shared<Constant>(0.0);
  // Unit adimenssional pressure
  auto p_unit = std::make_shared<Constant>(1.0);
  double p_atmvar = 1e5;  
  auto p_atm = std::make_shared<Constant>(p_atmvar/pression_carac);  Lproj->p_atm = p_atm; 
  // FIXME divided by 2
  auto p_width = std::make_shared<Constant>(p_widthvar/pression_carac); a->p_width = p_width; L->p_width = p_width;LInitN->p_width = p_width;aInitN->p_width = p_width;
  // Pressure increment
  double p_incr = (-0.01)/((double) modif1/(double) modif2); //modif = 1./1.
  int n_incr = (200*modif1) /modif2; //modif = 1./1.

  // std::cout<< p_incr << " "<<n_incr<<std::endl;

  auto p_ramp = std::make_shared<Constant>(p_incr);
  // Null displacement
  auto ui_null = std::make_shared<Constant>(0.0);



  // #### Initial profiles ####
  log.writeandshout("\n--- INITIAL PROFILES ---\n");

  // Solution at instant t
  auto wt = std::make_shared<Function>(W);     L->wt = wt;  LDual->w = wt; Lproj->wt = wt; 
  formQuasiStatLRel->wt = wt;formQuasiStatGRel->wt = wt;
  formQuasiStatLt->wt = wt;formQuasiStatGt->wt = wt;

  auto wk = std::make_shared<Function>(W); 
    a->wk = wk;  L->wk = wk;  LSum->wk = wk; LInitN->wk = wk;aInitN->wk = wk;
  formQuasiStatLRel->wk = wk; formQuasiStatGRel->wk = wk; LVolD->wk = wk;
  // Stress tensor and flux vectors
  auto wDual = std::make_shared<Function>(WDual);
  auto wDualt = std::make_shared<Function>(WDual);
  formQuasiStatLSpeedrel->wDual = wDual; formQuasiStatGSpeedrel->wDual = wDual; 
  formQuasiStatLSpeedrel->wDualt = wDualt; formQuasiStatGSpeedrel->wDualt = wDualt;formQuasiStatLSpeedt->wDualt = wDualt; formQuasiStatGSpeedt->wDualt = wDualt;
    // Newton-Raphson increment 
  auto delta_w = std::make_shared<Function>(W); 
    LSum->delta_w = delta_w;
  // Solution at k+1-th iteration of Newton-Raphson, approximating instant t+dt 
  auto wksol = std::make_shared<Function>(W); 
  // Pressures correction for atmospheric pressure boundary condition
  auto p_satproj = std::make_shared<Function>(Pproj);
  // Volumetric strain
  auto eps_v = std::make_shared<Function>(WVolD);

  // Initial liquid pressure profile
  auto lpp = std::make_shared<LiquidPressureProfile>();   lpp->set_p_atm(p_atmvar/pression_carac); LInit->pl0 = lpp; 
  // Initial gaz pressure profile
  auto gpp = std::make_shared<GazPressureProfile>();      gpp->set_p_atm(p_atmvar/pression_carac); LInit->pg0 = gpp; 

  // Instantiate BCs for deduced displacement field
  auto INITnoNormMoveBot = std::make_shared<DirichletBC>(W->sub(0)->sub(1), ui_null, bot);
  auto INITnoNormMoveRight= std::make_shared<DirichletBC>(W->sub(0)->sub(0), ui_null, right);
  auto INITnoNormMoveLeft= std::make_shared<DirichletBC>(W->sub(0)->sub(0), ui_null, left);
  std::vector<std::shared_ptr<const DirichletBC>> bcsNewtonInit = {{INITnoNormMoveBot, INITnoNormMoveRight, INITnoNormMoveLeft}};
  
  // Project initial fluid pressure profiles and deduce displacement field at equilibrium
  assemble(*AInit, *aInit);
  assemble(*bInit, *LInit);
  for (std::size_t i = 0; i < bcsNewtonInit.size(); i++){
    bcsNewtonInit[i]->apply(*AInit, *bInit);  bcsNewtonInit[i]->apply(*wt->vector());}

  log.write("Compute initial displacement field at equilibrium an project fluid pressure field");
  solver.solve(*AInit, *wt->vector(),*bInit); 
  
  *wk = Function(*wt);
  
  // Solve intial stress and fluxes field
  assemble(*ADual, *aDual);
  assemble(*bDual, *LDual);
  log.write("Compute wDual");
  solver.solve(*ADual, *wDual->vector(),*bDual); 

  // Solve pressure correctional terms
  assemble(*Aproj, *aproj);
  assemble(*bproj, *Lproj);
  log.write("Compute p_satproj");
  solver.solve(*Aproj, *p_satproj->vector(),*bproj);

  // First Saves
  if(SAVE_COMMON){
    pgfile << ((*wt)[2]);
    plfile << ((*wt)[1]);
    ufile << ((*wt)[0]);
  }
  if(SAVE_NEWTON){
    log.write("Save wk delta_w");
    pgfileDelta<< ((*delta_w)[2]);
    plfileDelta<< ((*delta_w)[1]);
    ufileDelta<< ((*delta_w)[0]);
  
    pgfileK << ((*wk)[2]);
    plfileK << ((*wk)[1]);
    ufileK << ((*wk)[0]);
  }

  if(SAVE_DUAL){
    log.write("Save wDual");
    stress << ((*wDual)[0]);
    fluxL << ((*wDual)[1]);
    fluxG << ((*wDual)[2]);
  }

  // #### Boundary Conditions #####
  log.writeandshout("\n--- BOUNDARY CONDITIONS ---\n");
  // Boundary Areas
  auto topSatLiquid = std::make_shared<TopSatLiquid>();

  // Area used to apply Neumann boundary condition
  auto subdomains = std::make_shared<MeshFunction<std::size_t>>(mesh, mesh->topology().dim()-1); 
  *subdomains = 0;
    L->ds = (subdomains);
  
  // Areas use to measure surface
  auto surfMeasLiq = std::make_shared<MeshFunction<std::size_t>>(mesh, mesh->topology().dim());
  *surfMeasLiq = 0;
    formSatSwitch->ds = surfMeasLiq;
  auto surfMeasGaz = std::make_shared<MeshFunction<std::size_t>>(mesh, mesh->topology().dim()-1);
  *surfMeasGaz = 0;
    formUnsatSwitch->ds = surfMeasGaz;
  double liquidSurf_cur(0.);
  double gazSurf_cur(0.);
  double liquidSurf_last(0.);
  double gazSurf_last(0.);

  // Set Connected to Air boundary attributes for liquid pressure
  (*topSatLiquid).set_functions(Function((*wk)[1]),Function((*wk)[2]),Function((*wDual)[1]));
  (*topSatLiquid).set_p_atm(p_atmvar/pression_carac*(1-*epsilon_l) ); 
  (*topSatLiquid).set_outflow_minimum(-0.);
  (*topSatLiquid).set_pressure_minimum(-5*0.000001*p_atmvar/pression_carac);
  
  // Link surfMeasLiq to topSatLiquid definition
  (*topSatLiquid).mark(*surfMeasLiq,1);

  auto noNormMoveBot = std::make_shared<DirichletBC>(W->sub(0)->sub(1), ui_null, bot);
  auto noNormMoveRight = std::make_shared<DirichletBC>(W->sub(0)->sub(0), ui_null, right);
  auto noNormMoveLeft = std::make_shared<DirichletBC>(W->sub(0)->sub(0), ui_null, left);

  auto PlNullBot = std::make_shared<DirichletBC>(W->sub(1), p_null, bot);//Center);
  auto PgNullBot = std::make_shared<DirichletBC>(W->sub(2), p_null, bot);
  auto PlRampBot = std::make_shared<DirichletBC>(W->sub(1), p_ramp, bot);//Center);

  auto PgNullTop = std::make_shared<DirichletBC>(W->sub(2), p_null, top);
  auto PlNullTop = std::make_shared<DirichletBC>(W->sub(1), p_null, top);

  // Copy pressure correctional values 
  log.write("Divide p_satproj");
  auto pl_satproj = std::make_shared<Function>((*p_satproj)[0]);
  auto pg_satproj = std::make_shared<Function>((*p_satproj)[1]);

  // Specific correction term to apply atmospheric pressure
  log.write("Insert p_satproj into Dirichlet BCs");
  auto PlTopSatLiquid = std::make_shared<DirichletBC>(W->sub(1), pl_satproj, topSatLiquid); 
  auto PgTopGaz = std::make_shared<DirichletBC>(W->sub(2), pg_satproj, top);
  auto PlNullTopSatLiquid = std::make_shared<DirichletBC>(W->sub(1), p_null, topSatLiquid); 
  auto PgNullTopGaz = std::make_shared<DirichletBC>(W->sub(2), p_null, top);

  // Wrap-up Dirichlet boundary condition for Newton-Raphson algorithm (t=0 && k>0)
  std::vector<std::shared_ptr<const DirichletBC>> bcsNewtonInitNiter = {{
    noNormMoveBot,
    noNormMoveRight,
    noNormMoveLeft,
    PlNullTopSatLiquid,
    PgNullTopGaz,
    PlNullBot}};

  // Wrap-up Dirichlet boundary condition for Newton-Raphson algorithm (t= 0 && k=0)
  std::vector<std::shared_ptr<const DirichletBC>> bcsNewtonInitN = {{
    noNormMoveBot,
    noNormMoveRight,
    noNormMoveLeft,
    PlTopSatLiquid,
    PgTopGaz,
    PlNullBot
  }};

  // Wrap-up Dirichlet boundary condition for Newton-Raphson algorithm (k>1)
  std::vector<std::shared_ptr<const DirichletBC>> bcsNewton = {{
    noNormMoveBot,
    noNormMoveRight,
    noNormMoveLeft,
    PlNullBot,
    PlNullTopSatLiquid,
    PgNullTopGaz}};

  // Wrap-up Dirichlet boundary condition for Newton-Raphson algorithm (t>1 && k=1)
  std::vector<std::shared_ptr<const DirichletBC>> bcsNewtonRamp = {{
    noNormMoveBot,
    noNormMoveRight,
    noNormMoveLeft,
    PlRampBot, 
    PlTopSatLiquid,
    PgTopGaz
  }};

  // Wrap-up Dirichlet boundary condition for Newton-Raphson algorithm (t=1 && k=1)
  std::vector<std::shared_ptr<const DirichletBC>> bcsNewtonIncr = {{
    noNormMoveBot,
    noNormMoveRight, 
    noNormMoveLeft,
    PlNullBot, 
    PlTopSatLiquid,
    PgTopGaz
  }};

  log.writeandshout("\n--- STARTING ALGORITHMS ---\n");

  // Alterned scheme convergence parameters
  auto tolerancySwitch = std::make_shared<Constant>(double((cvar*cvar/(cvar*cvar)) *0.01));
  double breakSwitch = true;
  auto iterSwitchmax = std::make_shared<Constant>(int(5));

  // Newton-Raphson convergence parameters
  auto tolerancyNewton = std::make_shared<Constant>(double(5e-5));
  auto iterNewtonmax = std::make_shared<Constant>(int(200));

  // Incremental Pressure is applied when psoeudo equilibrium is reached
  auto tolerancyQuasiStat = std::make_shared<Constant>(double(5e-4)); 
  auto tolerancyQuasiStatSpeedrelL = std::make_shared<Constant>(double(1e-4));
  auto tolerancyQuasiStatSpeedrelG = std::make_shared<Constant>(double(1e-4));
  auto tolerancyQuasiStatPressurerelL = std::make_shared<Constant>(double(1e-4));
  auto tolerancyQuasiStatPressurerelG = std::make_shared<Constant>(double(1e-4));

  double liquidQuasiStatRel(0.);double liquidQuasiStatt(0.);double liquidQuasiStatSpeedRel(0.);double liquidQuasiStatSpeedt(0.);
  double GasQuasiStatRel(0.);double GasQuasiStatt(0.);double GasQuasiStatSpeedRel(0.);double GasQuasiStatSpeedt(0.);
  bool critere = true;

  int cycle(0);
  int incr_count(0);

  // Loop counter
  int iterEuler(0);
  int iterSwitch(0);
  int iterNewton(0);
  int iterLagrange(0);
  double t = 0;
  bool breaker = false;
  int n_meas = 0;

  // Implicit Euler algorithm starts here
  while (t <= duration)
  {
    // Alterned BCs algorithm starts here
    iterSwitch = 0;

    *wDualt = Function(*wDual);
    do{
      // Update BCs
      assemble(*Aproj, *aproj);
      assemble(*bproj, *Lproj);
      log.write("Compute p_satproj");
      solver.solve(*Aproj, *p_satproj->vector(),*bproj);
      *pl_satproj = Function((*p_satproj)[0]);
      *pg_satproj = Function((*p_satproj)[1]);

      *PlNullBot = DirichletBC(W->sub(1), p_null, bot);//Center);
      *PlNullTopSatLiquid = DirichletBC(W->sub(1), p_null, topSatLiquid); 
      *PgNullTopGaz = DirichletBC(W->sub(2), p_null, top); 
      *PlTopSatLiquid = DirichletBC(W->sub(1), pl_satproj, topSatLiquid); 
      *PgTopGaz = DirichletBC(W->sub(2), pg_satproj, top); 
      *PlRampBot = DirichletBC(W->sub(1), p_ramp, bot);//Center);

      // Newton-Raphson algorithm starts here
      iterNewton = 0; 
      bool breaker = false;
       
      do{
        // Reset Sat./Part.sat interface definition
        *subdomains = 0;     

        // Update Matrices
        if(iterEuler==0){
          assemble(*AInitN, *aInitN);
          solver.set_operator(AInitN);
          assemble(*bInitN, *LInitN);        
        }else{
          assemble(*A, *a);
          solver.set_operator(A);
          assemble(*b, *L);
        }

        //Apply the correct set of BCs
        if (iterNewton ==0 and iterEuler == 0){
          // Heaviside
          log.write("Initial BC");
          for (std::size_t i = 0; i < bcsNewtonInitN.size(); i++){
            bcsNewtonInitN[i]->apply(*AInitN, *bInitN);
            bcsNewtonInitN[i]->apply(*delta_w->vector());
          }
        }
        else if (iterNewton >0 and iterEuler == 0){
          // Heaviside
          log.write("Initial Newton BC");
          for (std::size_t i = 0; i < bcsNewtonInitNiter.size(); i++){
            bcsNewtonInitNiter[i]->apply(*AInitN, *bInitN);
            bcsNewtonInitNiter[i]->apply(*delta_w->vector());
          }
        }
        else if (iterNewton ==0 and iterEuler == 1){
          // Heaviside
          log.write("Heaviside BC");
          for (std::size_t i = 0; i < bcsNewtonIncr.size(); i++){
            bcsNewtonIncr[i]->apply(*A, *b);
            bcsNewtonIncr[i]->apply(*delta_w->vector());
          }
        }
        else if (iterNewton ==0 and iterEuler>1){
          // Time-discrete Ramp
          log.write("Ramp BC");
          for (std::size_t i = 0; i < bcsNewtonRamp.size(); i++){
            bcsNewtonRamp[i]->apply(*A, *b);
            bcsNewtonRamp[i]->apply(*delta_w->vector());
          }
        }
        else{
          // Newton-Raphson BCs for multiple iterations
          log.write("Newton's algo BC");
          for (std::size_t i = 0; i < bcsNewton.size(); i++){
            bcsNewton[i]->apply(*A, *b);
            bcsNewton[i]->apply(*delta_w->vector());  
          }
        }
        // Find increment for in Newton-Raphson algorithm
        
        if(iterEuler==0){
          log.writeandshout("Compute delta_w for initial profiles");
          solver.solve(*AInitN, *delta_w->vector(),*bInitN);
        }else{
          log.write("Compute delta_w");
          solver.solve(*A, *delta_w->vector(),*b);}

        // Add the increment to the previous Newton term
        assemble(*ASum, *aSum);
        assemble(*bSum, *LSum);
        log.write("Compute Sum wk+delta_w");
        solver.solve(*ASum, *wksol->vector(),*bSum);
        *wk = *wksol;

        /**Deduce dual variables at k-th iteration of Newton-Rahpson algorithm
        * Note that LDual->w is wk and not wt as vl_k is used to update final
        * Boundary conditions
        */
        LDual->w = wk;
        assemble(*ADual, *aDual);
        assemble(*bDual, *LDual);
        log.write("Compute wDual");
        solver.solve(*ADual, *wDual->vector(),*bDual); //FIXME: custom solver

        // Save 
        if(SAVE_NEWTON){
          log.write("Save wk delta_w p_satproj");
          pgfileDelta << ((*delta_w)[2]);
          plfileDelta << ((*delta_w)[1]);
          ufileDelta  << ((*delta_w)[0]);
          pgfileK << ((*wk)[2]);
          plfileK << ((*wk)[1]);
          ufileK  << ((*wk)[0]);
          pprojfileK << (*p_satproj);
        }

        iterNewton++; 
        if( iterEuler == 0){
          log.write("Residual : " +std::to_string( (*bInitN).norm("l2") ) + " / Tol. : " +std::to_string( (double)*tolerancyNewton));
          breaker =((*bInitN).norm("l2") < *tolerancyNewton or iterNewton >= *iterNewtonmax);
        }
        else if(iterEuler > 0){
          log.write("Residual : " +std::to_string( (*b).norm("l2") ) + " / Tol. : " +std::to_string( (double)*tolerancyNewton));
          breaker =((*b).norm("l2") < *tolerancyNewton or iterNewton >= *iterNewtonmax);
        }
   
      // Newton-Raphson loop
      }while(not(breaker));
      log.writeandshout("k = "+ std::to_string(iterNewton)+ " /"+ std::to_string(double(*iterNewtonmax)));

      // Reset connected to air BCs
      (*topSatLiquid).set_functions(Function((*wk)[1]),Function((*wk)[2]),Function((*wDual)[1]));
      
      // Mark surMeas for convergence criterion
      *surfMeasLiq = 0;
      *surfMeasGaz = 0;
      (*topSatLiquid).mark(*surfMeasLiq,1);

      // measure area of the newly found BCs
      liquidSurf_cur = assemble(*formSatSwitch);
      gazSurf_cur = assemble(*formUnsatSwitch);
      log.write("Liquid area :" +to_string(assemble(*formSatSwitch))+" | Gaz area :" +to_string(assemble(*formSatSwitch)));
      log.write("test : "+to_string(abs(liquidSurf_cur-liquidSurf_last)+abs(gazSurf_cur-gazSurf_last))+ " < " +to_string( (*tolerancySwitch) )); 

      // Alterned Scheme for connected to air BCs criteria
      breakSwitch = iterSwitch>(*iterSwitchmax) or (iterSwitch>0 and (abs(liquidSurf_cur-liquidSurf_last)+abs(gazSurf_cur-gazSurf_last)<(*tolerancySwitch)));// and (not(abs(stddevSwitch_cur-stddevSwitch_last)<1e-6) or not(nbNodeSwitch_cur==nbNodeSwitch_last)));
      if(not(breakSwitch)){  
        *wk = Function(*wt);
        if(iterSwitch==0){
          log.writeandshout("Switching BCs - First Guess");
        }else{
          log.writeandshout("Switching BCs - Refining");
        }
      }

      liquidSurf_last = (double) liquidSurf_cur;
      gazSurf_last = (double) gazSurf_cur;
      
    // Alterned Scheme loop
    iterSwitch++;
    }while(not(breakSwitch));

    // Save
    if(iterEuler% save_freq==0){
      if(SAVE_COMMON){
        pgfile << ((*wt)[2]);
        plfile << ((*wt)[1]);
        ufile << ((*wt)[0]);
        subdomainLiq_file << *surfMeasLiq;
        subdomainGaz_file << *surfMeasGaz;
      }

      if(SAVE_DUAL){  
        log.write("Save Dual");
        stress << ((*wDual)[0]);
        fluxL << ((*wDual)[1]);
        fluxG << ((*wDual)[2]);
      }
      
      if(SAVE_POINT){
        dolfin::Array<double> values_temp(1);
        std::string buffer_res = "";
        assemble(*AVolD, *aVolD);
        assemble(*bVolD, *LVolD);
        log.write("Compute eps_v");
        solver.solve(*AVolD, *eps_v->vector(),*bVolD); //FIXME: custom solver
        log.write("Save some eps_v values");

        std::ostringstream streamObj;

        eps_v->eval(values_temp,point0); streamObj.str("");streamObj.clear(); streamObj << values_temp[0];
        buffer_res = buffer_res + streamObj.str();
        buffer_res = buffer_res + "; ";
        eps_v->eval(values_temp,point1); streamObj.str("");streamObj.clear(); streamObj << values_temp[0];
        buffer_res = buffer_res + streamObj.str();
        buffer_res = buffer_res + "; ";
        eps_v->eval(values_temp,point2); streamObj.str("");streamObj.clear(); streamObj << values_temp[0];
        buffer_res = buffer_res + streamObj.str();
        buffer_res = buffer_res + "; ";
        eps_v->eval(values_temp,point3); streamObj.str("");streamObj.clear(); streamObj << values_temp[0];
        buffer_res = buffer_res + streamObj.str();
        buffer_res = buffer_res + "; ";
        eps_v->eval(values_temp,point4); streamObj.str("");streamObj.clear(); streamObj << values_temp[0];
        buffer_res = buffer_res + streamObj.str();
        
        eps_vpointfile.write(buffer_res);

        buffer_res = "";
        std::string buffer_res_bis = "";

        log.write("Save some pl and pg values");

        Function plpointtemp((*wk)[1]);    plpointtemp.set_allow_extrapolation(true);
        Function pgpointtemp((*wk)[2]);    pgpointtemp.set_allow_extrapolation(true);
        plpointtemp.eval(values_temp,point0); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        pgpointtemp.eval(values_temp,point0); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
        buffer_res = buffer_res + "; "; buffer_res_bis = buffer_res_bis + "; ";
        plpointtemp.eval(values_temp,point1); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        pgpointtemp.eval(values_temp,point1); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
        buffer_res = buffer_res + "; "; buffer_res_bis = buffer_res_bis + "; ";
        plpointtemp.eval(values_temp,point2); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        pgpointtemp.eval(values_temp,point2); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
        buffer_res = buffer_res + "; "; buffer_res_bis = buffer_res_bis + "; ";
        plpointtemp.eval(values_temp,point3); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        pgpointtemp.eval(values_temp,point3); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
        buffer_res = buffer_res + "; "; buffer_res_bis = buffer_res_bis + "; ";
        plpointtemp.eval(values_temp,point4); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        pgpointtemp.eval(values_temp,point4); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
                
        plpointfile.write(buffer_res);
        pgpointfile.write(buffer_res_bis);
        
        buffer_res = "";
        buffer_res_bis = "";

        log.write("Save displacement values");

        Function uxpointtemp((*wk)[0][0]);    uxpointtemp.set_allow_extrapolation(true);
        Function uypointtemp((*wk)[0][1]);    uypointtemp.set_allow_extrapolation(true);
        uxpointtemp.eval(values_temp,point0); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        uypointtemp.eval(values_temp,point0); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
        buffer_res = buffer_res + "; "; buffer_res_bis = buffer_res_bis + "; ";
        uxpointtemp.eval(values_temp,point1); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        uypointtemp.eval(values_temp,point1); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
        buffer_res = buffer_res + "; "; buffer_res_bis = buffer_res_bis + "; ";
        uxpointtemp.eval(values_temp,point2); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        uypointtemp.eval(values_temp,point2); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
        buffer_res = buffer_res + "; "; buffer_res_bis = buffer_res_bis + "; ";
        uxpointtemp.eval(values_temp,point3); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        uypointtemp.eval(values_temp,point3); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
        buffer_res = buffer_res + "; "; buffer_res_bis = buffer_res_bis + "; ";
        uxpointtemp.eval(values_temp,point4); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res = buffer_res + streamObj.str();
        uypointtemp.eval(values_temp,point4); streamObj.str("");streamObj.clear(); 
        streamObj << values_temp[0];buffer_res_bis = buffer_res_bis + streamObj.str();
                
        uxpointfile.write(buffer_res);
        uypointfile.write(buffer_res_bis);

      }
      log.writeandshout("t = " +to_string(t/temps_carac) + " / "+ to_string(duration/temps_carac)+ "\n\n");
    }

    liquidQuasiStatRel = assemble(*formQuasiStatLRel);GasQuasiStatRel = assemble(*formQuasiStatGRel);
    liquidQuasiStatt = assemble(*formQuasiStatLt);GasQuasiStatt = assemble(*formQuasiStatGt);
    liquidQuasiStatSpeedRel = assemble(*formQuasiStatLSpeedrel);GasQuasiStatSpeedRel = assemble(*formQuasiStatGSpeedrel);
    liquidQuasiStatSpeedt = assemble(*formQuasiStatLSpeedt);GasQuasiStatSpeedt = assemble(*formQuasiStatGSpeedt);
    log.writeandshout("liqPresssureL2rel :" + to_string(sqrt(liquidQuasiStatRel))    + " < " + to_string(*tolerancyQuasiStatPressurerelL) + "\ngasPresssureL2rel :" + to_string(sqrt(GasQuasiStatRel))     + " < " + to_string(*tolerancyQuasiStatPressurerelG)+" ?");
    log.writeandshout("liqSpeedL2rel :" + to_string(sqrt(liquidQuasiStatSpeedRel)) + " < " + to_string(*tolerancyQuasiStatSpeedrelL)+ "\ngasPressureL2rel :" + to_string(sqrt(GasQuasiStatSpeedRel))+ " < " + to_string(*tolerancyQuasiStatSpeedrelG)+" ?");
    
    //Stop criterion for whole simulation
    critere = sqrt(liquidQuasiStatRel)<*tolerancyQuasiStatPressurerelL;
    critere = critere and sqrt(GasQuasiStatRel) < *tolerancyQuasiStatPressurerelG;
    critere = critere and sqrt(liquidQuasiStatSpeedRel) < (*tolerancyQuasiStatSpeedrelL);
    critere = critere and sqrt(GasQuasiStatSpeedRel) < (*tolerancyQuasiStatSpeedrelG);

    log.writeandshout("critere : " + to_string(critere));

    if(SAVE_PROFILINIT and ((cycle==0 and incr_count == 0 and n_meas ==0) or ( cycle==1 and incr_count ==n_incr and n_meas ==1)) ){
      n_meas++;
      log.write("Save vertical pressure profiles - n_meas = " + to_string(n_meas));
      std::ostringstream streamObj_init;
      dolfin::Array<double> values_init(1);
      double p_init[3] = {2.,0.,0. }; 
      dolfin::Array<double> point_init(sizeof(p_init)/sizeof(double),p_init);

      Function plinitprofile((*wk)[1]);    plinitprofile.set_allow_extrapolation(true);
      Function pginitprofile((*wk)[2]);    pginitprofile.set_allow_extrapolation(true);
      for(double deep = 0; deep<4.001;deep += 0.004){
        point_init[1] = deep;
        plinitprofile.eval(values_init,point_init);
        streamObj_init.str("");streamObj_init.clear();streamObj_init << values_init[0];
        plinitfile.write("("+to_string(p_init[0]) + ","+to_string(p_init[1]) + ","+ to_string(p_init[2]) + "); " + streamObj_init.str() );
        pginitprofile.eval(values_init,point_init);
        streamObj_init.str("");streamObj_init.clear();streamObj_init << values_init[0];
        pginitfile.write("("+to_string(p_init[0]) + ","+to_string(p_init[1]) + ","+ to_string(p_init[2]) + "); " + streamObj_init.str() );
      }
      plinitfile.write("--");
      pginitfile.write("--");
    }

    if(SAVE_PROFILMID and ((cycle==1 and incr_count == n_incr and n_meas ==2 and (critere)) or ( cycle==2 and incr_count ==0 and n_meas ==3)) ){
      n_meas++; 
      log.write("Save vertical pressure profiles - n_meas = " + to_string(n_meas) +" cycle "+ to_string(cycle));
      std::ostringstream streamObj_mid;
      dolfin::Array<double> values_mid(1);
      double p_mid[3] = {2.,0.,0. }; 
      dolfin::Array<double> point_mid(sizeof(p_mid)/sizeof(double),p_mid);

      Function plmidprofile((*wk)[1]);    plmidprofile.set_allow_extrapolation(true);
      Function pgmidprofile((*wk)[2]);    pgmidprofile.set_allow_extrapolation(true);
      for(double deep = 0; deep<4.001;deep += 0.004){
        point_mid[1] = deep;
        plmidprofile.eval(values_mid,point_mid);
        streamObj_mid.str("");streamObj_mid.clear();streamObj_mid << values_mid[0];
        plmidfile.write("("+to_string(p_mid[0]) + ","+to_string(p_mid[1]) + ","+ to_string(p_mid[2]) + "); " + streamObj_mid.str() );
        pgmidprofile.eval(values_mid,point_mid);
        streamObj_mid.str("");streamObj_mid.clear();streamObj_mid << values_mid[0];
        pgmidfile.write("("+to_string(p_mid[0]) + ","+to_string(p_mid[1]) + ","+ to_string(p_mid[2]) + "); " + streamObj_mid.str() );
      }
      plmidfile.write("--");
      pgmidfile.write("--");
    }

    if(cycle==2 and incr_count ==0 and critere){
      if(SAVE_PROFILFIN and cycle==2 and incr_count ==0){
        n_meas++;
        log.write("Save vertical pressure profiles - n_meas = " + to_string(n_meas));
        std::ostringstream streamObj_fin;
        dolfin::Array<double> values_fin(1);
        double p_fin[3] = {2.,0.,0. }; 
        dolfin::Array<double> point_fin(sizeof(p_fin)/sizeof(double),p_fin);

        Function plfinprofile((*wk)[1]);    plfinprofile.set_allow_extrapolation(true);
        Function pgfinprofile((*wk)[2]);    pgfinprofile.set_allow_extrapolation(true);
        for(double deep = 0; deep<4.001;deep += 0.004){
          point_fin[1] = deep;
          plfinprofile.eval(values_fin,point_fin);
          streamObj_fin.str("");streamObj_fin.clear();streamObj_fin << values_fin[0];
          plfinfile.write("("+to_string(p_fin[0]) + ","+to_string(p_fin[1]) + ","+ to_string(p_fin[2]) + "); " + streamObj_fin.str() );
          pgfinprofile.eval(values_fin,point_fin);
          streamObj_fin.str("");streamObj_fin.clear();streamObj_fin << values_fin[0];
          pgfinfile.write("("+to_string(p_fin[0]) + ","+to_string(p_fin[1]) + ","+ to_string(p_fin[2]) + "); " + streamObj_fin.str() );
        }
      }
    }

    if(cycle==2 and incr_count ==0 and (critere)){break;}

    if((incr_count == n_incr or incr_count == 0) and not(critere)){
      *p_ramp = *p_null;
    }else{
      if(cycle%2 == 0 ){incr_count++;}else{incr_count--;}
      *p_ramp = p_incr;      
    }

    if(incr_count == n_incr and cycle%2 == 0 ){p_incr = - p_incr;cycle++;}

    if(incr_count == 0 and cycle%2 == 1){p_incr = - p_incr;cycle++;}

    *wt = Function(*wk);

    // Implicit Euler loop
    iterEuler++;
    t += dt;
  
  }

  log.writeandshout("JOB DONE!");

  log.writeandshout("La simulation temps réel dure : " +std::to_string(t) +" s.");

  #else
  std::cout << "This requires DOLFIN to be configured with PETSc." << std::endl;
  #endif

  return 0;
};
