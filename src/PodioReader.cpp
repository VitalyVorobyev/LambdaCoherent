#include "PodioReader.h"

// Data model
#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"

#include <iostream>
#include <cmath>
#include <cassert>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using std::make_unique;

using linal::LVect;
using linal::Vect;

const string PodioReader::c_gpclCol("allGenParticles");
const string PodioReader::c_gvtxCol("allGenVertices");

PodioReader::PodioReader(const std::string& ifile) {
    m_store = make_unique<podio::EventStore>();
    m_reader = make_unique<podio::ROOTReader>();
    m_reader->openFile(ifile);
    m_store->setReader(m_reader.get());
    m_nevents = m_reader->getEntries();

    cout << "### PodioReader Initialized ###" << endl
         << "  " << m_nevents << " events" << endl;
}

double energy(double mass, const Vect& p3) {
    return sqrt(mass*mass + p3.r2());
}

void PodioReader::event() {
    m_p4map.clear();
    m_gammaEnergy = 0.;

    const fcc::MCParticleCollection* particles;
    if (!m_store->get(c_gpclCol, particles)) {
        cerr << "MCParticles" << endl;
        assert("");
    }

    for (const auto& pcl : *particles) {
        Vect p3(pcl.p4().px, pcl.p4().py, pcl.p4().pz);
        if (pcl.pdgId() == 22) {
            m_gammaEnergy += p3.r();
        } else {
            m_p4map.emplace(pcl.pdgId(), LVect(energy(pcl.p4().mass, p3), p3));
            // cout << pcl.pdgId() << ": " << LVect(energy(pcl.p4().mass, p3), p3) << endl;
        }
    }
    // cout << "EGamma: " << m_gammaEnergy << endl;

    const fcc::GenVertexCollection* vertices;
    if (!m_store->get(c_gvtxCol, particles)) {
        cerr << "GenVertices" << endl;
        assert("");
    }

    m_store->clear();
    m_reader->endOfEvent();
}

const LVect& PodioReader::getp4(int pdgid) const {
    auto it = m_p4map.find(pdgid);
    if (it == m_p4map.end()) {
        cerr << "PodioReader: can't find particle id " << pdgid << endl;
        assert("");
    }
    return it->second;
}
