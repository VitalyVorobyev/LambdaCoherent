#include "PodioReader.h"
#include "CohMSqLL.h"
#include "LLKine.h"

#include <string>
#include <iostream>
#include <fstream>

using std::string;
using std::cout;
using std::endl;

using std::ofstream;

void writeEvent(const PodioReader& pr, ofstream& os) {
    static CohMsqLL mp( 1., 1., 0., 0.75, 0.75);
    static CohMsqLL mn(-1., 1., 0., 0.75, 0.75);
    static CohMsqLL mz( 0., 1., 0., 0.75, 0.75);

    os << "Event" << endl
       << pr.p4p() << endl
       << pr.p4pbar() << endl
       << pr.p4pip()  << endl
       << pr.p4pin() << endl;
    auto xi = LLKine::xi(pr.p4p(), pr.p4pin(), pr.p4pbar(), pr.p4pip());
    for (auto x : xi) {
        os << x << " ";
    }
    os << endl
       << mp(pr.p4p(), pr.p4pbar(), pr.p4pip(), pr.p4pin()) << " "
       << mn(pr.p4p(), pr.p4pbar(), pr.p4pip(), pr.p4pin()) << " "
       << mz(pr.p4p(), pr.p4pbar(), pr.p4pip(), pr.p4pin()) << " "
       << endl;
}

int main(int argc, char** argv) {
    string fname = "/home/vitaly/jpsipppipi.root";
    PodioReader pr(fname);
    ofstream os("ll.dat", ofstream::out);

    for (size_t i = 0; i < pr.nEvt(); i++) {
    // for (size_t i = 0; i < 100000; i++) {
        if (!(i+1 % 10000)) {
            cout << i+1 << " events" << endl;
        }
        pr.event();
        writeEvent(pr, os);
        // cout << m(pr.p4p(), pr.p4pbar(), pr.p4pip(), pr.p4pin()) << endl;
    }

    return 0;
}
