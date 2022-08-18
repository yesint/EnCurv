/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EnCurv method for maintaining membrane curvature.
    
    (c) Semen Yesylevskyy, 2020. yesint4@gmail.com
    
    Supposed to be compiled with PLUMED v2.5 or higher
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "core/Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

namespace PLMD {
namespace colvar {


inline void smooth(double edge0, double edge1, double x, double& s, double& ds) {
    x = M_PI * (x - edge0) / (edge1 - edge0);
    s = 0.5*(1.0-cos(x));
    ds = 0.5*sin(x)* M_PI / (edge1 - edge0);
}


struct Bin {
    double ang;
    double r_mean;
    double N;
    double wm,Z;
    double P;

    Bin(): N(0), wm(0), Z(0) {}

    void clear(){
        N=0; wm=0; Z=0;
    }
};


class EnCurv:
  public Colvar
{
private:
  Vector center;  
  double r_target;
  double cap_size;
  double x_span;

  unsigned N;  
  vector<double> angles;
  vector<double> radii;
  vector<Vector> vectors, tangents;

  unsigned Nbins;
  vector<Bin> bins;
  vector<array<unsigned,2>> atom_bin;
  vector<array<double,2>> atom_s;
  vector<array<double,2>> atom_sd;

  bool no_phi_force;
  bool is_tube;

public:
  explicit EnCurv(const ActionOptions&ao);
  void calculate();
  static void registerKeywords( Keywords& keys );
  int bending_axis;
  int X_ax, Z_ax;
};


PLUMED_REGISTER_ACTION(EnCurv,"ENCURV")


void EnCurv::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms");
  keys.add("compulsory","R","Desired radius.");
  keys.add("compulsory","AXIS","1","Bending axis X=0,Y=1,Z=2. Default: 1.");
  keys.add("compulsory","NBINS","50","Number of bins.");
  keys.addFlag("NO_PHI_FORCE",false,"Exclude tangential forces for equilibration.");
  keys.add("compulsory","CAP_SIZE","2.5","Caps to skip at the ends of the membrane in nm.");
  keys.add("compulsory","XSPAN","0","Defines sector for biasing. Disables CAP_SIZE. Useful for periodic bilayers.");
  keys.addFlag("TUBE",false,"Tubular geometry.");
  
  keys.addOutputComponent("val","val","Value");
  keys.addOutputComponent("rmsd","rmsd","RMSD");
  keys.addOutputComponent("angle","angle","ANGLE");
}


EnCurv::EnCurv(const ActionOptions&ao):
    PLUMED_COLVAR_INIT(ao)
{  
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("at least one atom should be specified");
  
  parse("R",r_target);
  parse("NBINS",Nbins);
  parse("AXIS",bending_axis);
  parse("CAP_SIZE",cap_size);
  parse("XSPAN",x_span);
  parseFlag("NO_PHI_FORCE",no_phi_force);
  parseFlag("TUBE",is_tube);

  if(no_phi_force){
      log << "  PHI forces are disabled by the user!\n";
  } else {
      log << "  PHI forces are enabled.\n";
  }  

  checkRead();
  
  // Geometry is defined for bending in XZ plane
  // Mapping to actual system geometry is done by defining indexes corresponding to real coordinates
  if(bending_axis==1){
    X_ax = 0;
    Z_ax = 2;
  } else if(bending_axis==2) {
    X_ax = 0;
    Z_ax = 1;  
  } else if(bending_axis==0) {
    X_ax = 1;
    Z_ax = 2;  
  }

  addComponentWithDerivatives("val"); componentIsNotPeriodic("val");
  addComponent("rmsd"); componentIsNotPeriodic("rmsd");
  addComponentWithDerivatives("angle"); componentIsNotPeriodic("angle");

  N = atoms.size();

  angles.resize(N);
  radii.resize(N);
  vectors.resize(N);
  tangents.resize(N);

  bins.resize(Nbins);
  atom_bin.resize(N);
  atom_s.resize(N);
  atom_sd.resize(N);

  requestAtoms(atoms);
}


void bin_add(vector<Bin>& v1,const vector<Bin>& v2){
    for(unsigned i=0;i<v1.size();++i){
        v1[i].wm += v2[i].wm;
        v1[i].Z += v2[i].Z;
    }
}


void init_bins(vector<Bin>& v1){
    for(unsigned i=0;i<v1.size();++i) v1[i].clear();
}


void EnCurv::calculate() {
    // Set center to first atom
    center = getPosition(0);
    center[bending_axis] = 0.0;

    // In case of predefined sector set X to box center
    // to accomodate for box changes
    if(x_span){
        center[X_ax] = 0.5*getBox()(X_ax,X_ax);
    }

    double min_ang = 1e10, max_ang=-1e10;
    double mean_ang = 0.0; // Average angle
    double total_mass = 0.0;

    Vector axis_vector(0,0,0);
    if(bending_axis==1){
        axis_vector[bending_axis] = 1.0;
    } else {
        axis_vector[bending_axis] = -1.0;
    }

    for(unsigned i=1; i<N; i++){
        Vector p = getPosition(i);
        p[bending_axis]=0.0;
        // Vector from center to atom
        vectors[i] = delta(center,p);
        double r = vectors[i].modulo();
        vectors[i] /= r; // Normalize vectors

        // Angle relative to Z axis from -pi to pi
        double ang = atan2(vectors[i][X_ax],vectors[i][Z_ax]);

        mean_ang += getMass(i)*ang;
        total_mass += getMass(i);

        angles[i] = ang;
        radii[i] = r;

        // tangent vector. Note -1, it matters to get correct direction!
        // tangent is to right (in direction of phi increase)
        tangents[i] = crossProduct(vectors[i],axis_vector);

        if(ang<min_ang) min_ang = ang;
        if(ang>max_ang) max_ang = ang;
    }

    if(x_span){
        double a = asin(0.5*x_span/r_target);
        min_ang = -a;
        max_ang =  a;
        // Mean angle is assumed zero
        mean_ang = 0.0;
    } else {
        min_ang += cap_size/r_target;
        max_ang -= cap_size/r_target;
        // Get mass-weigted average angle
        mean_ang /= total_mass;
    }
    
    // In case of tube set angles to +/-pi
    if(is_tube){
        min_ang = -M_PI;
        max_ang = +M_PI;
    }


    // Divide into bins and compute radial centers
    for(unsigned i=0; i<Nbins; i++) bins[i].clear();
    double d_ang = (max_ang-min_ang)/float(Nbins);
    
    // Set bin angles
    for(unsigned i=0; i<Nbins; i++) bins[i].ang = min_ang+d_ang*(i+0.5);

    // Each atom contributes to two adjucent bins
    for(unsigned i=1; i<N; i++){
        // Current bin
        unsigned b,b1,b2;
        
        b = floor((angles[i]-min_ang)/d_ang);
        if(x_span){
            // For predefined sector ignore atoms outside the sector
            if(angles[i]<min_ang || angles[i]>max_ang) continue;
        } else if(is_tube) {
            // For tube wrap bins around periodically
            if(angles[i]<min_ang) b = Nbins-1;
            if(angles[i]>max_ang) b = 0;
        } else {
            if(angles[i]<min_ang) b=0;
            if(angles[i]>max_ang) b = Nbins-1;
        }
        
        // Set adjucent bins        
        double side = angles[i]-bins[b].ang;
        if(is_tube){
            // For tube wrap around
            if(side<=0){
                b1 = (b>0) ? b-1 : Nbins-1;
                b2 = b;
            } else {
                b1 = b;
                b2 = (b<Nbins-1) ? b+1 : 0;
            }
            // Order bins b1<b2
            if(b1>b2) std::swap(b1,b2);
        } else {
            if(side<=0){
                b1 = (b>0) ? b-1 : b;
                b2 = b;
            } else {
                b1 = b;
                b2 = (b<Nbins-1) ? b+1 : b;
            }
        }
        
        atom_bin[i][0] = b1;
        atom_bin[i][1] = b2;
        
        double m = getMass(i);
        double s,ds;

        if(b1==b2){
            if(!x_span){
                // This is the edges, so just apply for the current bin without interpolation
                atom_s[i][0] = atom_s[i][1] = 1.0;
                bins[b].Z += m/radii[i];
                bins[b].wm += m;
                bins[b].N += 1;
                // Angular component is zero
                atom_sd[i][0] = atom_sd[i][1] = 0.0;
            } else {
                // For predifined sector last half-bind should be ignored
                // For left bin (1-c) is applied, for right bin (c) is applied
                if(side<=0){
                    smooth(bins[b1].ang-d_ang, bins[b1].ang, angles[i], s, ds);
                    bins[b1].Z += s * m / radii[i];
                    bins[b1].wm +=  s * m;
                    bins[b1].N +=  s;
                    atom_s[i][0] = 0.0;
                    atom_s[i][1] = s;
                    atom_sd[i][0] = 0.0;
                    atom_sd[i][1] = +ds;
                } else {
                    smooth(bins[b1].ang, bins[b1].ang+d_ang, angles[i], s, ds);
                    bins[b1].Z += (1.0-s) * m / radii[i];
                    bins[b1].wm +=  (1.0-s) * m;
                    bins[b1].N +=  (1.0-s);
                    atom_s[i][0] = (1.0-s);
                    atom_s[i][1] = 0.0;
                    atom_sd[i][0] = -ds;
                    atom_sd[i][1] = 0.0;
                }

            }
        } else {
            // In the middle interpolate
            smooth(bins[b1].ang, bins[b2].ang, angles[i], s, ds);
            // For left bin (1-c) is applied, for right bin (c) is applied
            bins[b1].Z += (1.0-s) * m / radii[i];
            bins[b1].wm +=  (1.0-s) * m;
            bins[b1].N +=  (1.0-s);

            bins[b2].Z += s * m / radii[i];
            bins[b2].wm +=  s * m;
            bins[b2].N +=  s;

            atom_s[i][0] = (1.0-s);
            atom_s[i][1] = s;

            atom_sd[i][0] = -ds;
            atom_sd[i][1] = +ds;
        }
    }

    // Compute r_com for all bins
    for(unsigned b=0; b<Nbins; b++){
        if(bins[b].Z>0){
            bins[b].r_mean = bins[b].wm/bins[b].Z;
            bins[b].P = (bins[b].r_mean-r_target)/bins[b].Z;
        } else {
            bins[b].r_mean = r_target;
            bins[b].P = 0.0;
        }
    }

    Value* v_ptr=getPntrToComponent("val");

    // For each atom find deviation from target_r
    for(unsigned i=1; i<N; i++){
        unsigned b1 = atom_bin[i][0];
        unsigned b2 = atom_bin[i][1];

        double rv,ra;
        double m = getMass(i);

        if(b1!=b2){
            // Midlle of bicelle with interpolating between two bins
            // Radial component
            rv  = (m/(radii[i]*radii[i])) * (  bins[b1].P * bins[b1].r_mean * atom_s[i][0]
                                             + bins[b2].P * bins[b2].r_mean * atom_s[i][1] ) ;

            // Angular component	    
            
            //ra  = m * (  bins[b1].P * atom_sd[i][0] * (1.0-bins[b1].r_mean/radii[i])/bins[b1].r_mean
            //           + bins[b2].P * atom_sd[i][1] * (1.0-bins[b2].r_mean/radii[i])/bins[b2].r_mean );
                       
            ra  = m * (  bins[b1].P * atom_sd[i][0] * (1.0-bins[b1].r_mean/radii[i])/radii[i]
                       + bins[b2].P * atom_sd[i][1] * (1.0-bins[b2].r_mean/radii[i])/radii[i] );
            
        } else {
            // Half-bins at the ends of bicelle. No interpolation.
            rv  = bins[b1].P * bins[b1].r_mean * m * atom_s[i][0] / pow(radii[i],2);
            // Angular component is zero here
            ra = 0.0;
        }

        if(!no_phi_force){
            setAtomsDerivatives(v_ptr,i, vectors[i]*rv - tangents[i]*ra);
        } else {
            setAtomsDerivatives(v_ptr,i, vectors[i]*rv);
        }
    }

    // Set derivatives for zero for pivot atom
    setAtomsDerivatives(v_ptr,0,Vector(0,0,0));
    // Set box derivs
    setBoxDerivativesNoPbc(v_ptr);
    // Set value
    v_ptr->set(r_target+1.0); // Gives force at 1 nm if bias is set to r_target

    // Angular potential    
    if(!x_span){ // Not used for predefined sector
        Value* ang_ptr=getPntrToComponent("angle");
        // Set value
        ang_ptr->set(mean_ang);
        // Set derivatives
        for(unsigned i=1; i<N; i++){
            setAtomsDerivatives(ang_ptr,i,  -tangents[i] * getMass(i) / total_mass);
        }
        // Zero for pivot atom
        setAtomsDerivatives(ang_ptr,0,Vector(0,0,0));
        // Set box derivs
        setBoxDerivativesNoPbc(ang_ptr);
    }

    // RMSD
    double rmsd = 0.0;
    for(unsigned b=0; b<Nbins; b++){
        rmsd += std::pow(bins[b].r_mean-r_target,2.0);
    }
    Value* rmsd_ptr=getPntrToComponent("rmsd");
    rmsd_ptr->set( sqrt(rmsd/double(Nbins)) );

    // Dump some data to log file
    long int t = getStep();
    if(t%1000==0){
        // Radii
        log << t << ": ";
        for(unsigned b=0; b<Nbins; b++){
           log << bins[b].r_mean << " ";
        }
        log << "\n";

        // Number of atoms per bin
        double Nmean = 0.0;
        log << "Atoms per bin: ";
        for(unsigned b=0; b<Nbins; b++){
           log << bins[b].N << " ";
           Nmean += bins[b].N;
        }
        log << "\n";
        log << "Mean atoms per bin: " << Nmean/double(Nbins) << "\n";
    }
  }
}

}

