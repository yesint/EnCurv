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


class MembranePore:
        public Colvar
{
private:

    double radius, force_constant;
    unsigned N;
    int axis;

public:
    explicit MembranePore(const ActionOptions&ao);
    void calculate();
    static void registerKeywords( Keywords& keys );
};


PLUMED_REGISTER_ACTION(MembranePore,"MEMBRANEPORE")


void MembranePore::registerKeywords(Keywords& keys) {
    Colvar::registerKeywords( keys );
    keys.add("atoms","ATOMS","Atoms to force");
    keys.add("compulsory","R","Desired radius.");
    keys.add("compulsory","AXIS","Pore along this axis. X=0,Y=1,Z=2");

    keys.addOutputComponent("val","val","Value");
}


MembranePore::MembranePore(const ActionOptions&ao):
    PLUMED_COLVAR_INIT(ao)
{  
    vector<AtomNumber> atoms;
    parseAtomList("ATOMS",atoms);
    if(atoms.size()==0) error("at least one atom should be specified");

    parse("R",radius);
    parse("AXIS",axis);
    checkRead();

    addComponentWithDerivatives("val"); componentIsNotPeriodic("val");

    requestAtoms(atoms);
    N = atoms.size();
}

void MembranePore::calculate() {
    // Set desired center to first atom
    Vector center = getPosition(0);

    Value* v_ptr=getPntrToComponent("val");

    // Cycle over atoms and find those inside the pore radius
    for(unsigned i=1; i<N; i++){
        Vector p = getPosition(i);
        Vector d = pbcDistance(center,p);
        double factor = 0.0;

        if(d[axis]<0){ // Below center
            d[axis] = 0.0;
            double r = d.modulo();
            if(r<radius){
                // Inside the pore
                factor = 1.0-r/radius;
            }
        }
        setAtomsDerivatives(v_ptr,i, d*factor);
    }

    // Set derivatives for zero for pivot atom
    setAtomsDerivatives(v_ptr,0,Vector(0,0,0));
    // Set box derivs
    setBoxDerivativesNoPbc(v_ptr);
    // Set value
    v_ptr->set(radius); // Gives force at 1 nm if bias is set to radius
}

}

}

