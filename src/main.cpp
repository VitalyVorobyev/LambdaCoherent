#include "PodioReader.h"
#include "CohMSqLL.h"

#include <string>
#include <iostream>

using std::string;
using std::cout;
using std::endl;

int main(int argc, char** argv) {
    string fname = "/home/vitaly/jpsipppipi.root";
    PodioReader pr(fname);
    CohMsqLL m(1., 1., 0.5*3.1415, 0.75, 0.75);
    // LeviCivita lc;
    // cout << "Levi: " << lc(0, 1, 2, 3) << endl;

    // for (size_t i = 0; i < pr.nEvt(); i++) {
    for (size_t i = 0; i < 100; i++) {
        if (!(i % 1000)) {
            cout << i << " events" << endl;
        }
        pr.event();
        cout << m(pr.p4p(), pr.p4pbar(), pr.p4pip(), pr.p4pin()) << endl;
    }

    return 0;
}
