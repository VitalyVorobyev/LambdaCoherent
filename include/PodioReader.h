/** **/

#ifndef PODIOREADER_H__
#define PODIOREADER_H__

// libLinal
#include "lvect.h"

// podio
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include <memory>
#include <map>
#include <string>

class PodioReader {
    /** Event Store **/
    std::unique_ptr<podio::EventStore> m_store;
    /** ROOT backend for event reading **/
    std::unique_ptr<podio::ROOTReader> m_reader;
    /** Momentum table **/
    std::map<int, linal::LVect> m_p4map;

    /** Number of events in store **/
    size_t m_nevents;

    /** Energy caried by FSR photons **/
    double m_gammaEnergy;

    /** Find particle in the momentum table **/
    const linal::LVect& getp4(int) const;

    /** Collections names **/
    static const std::string c_gpclCol;
    static const std::string c_gvtxCol;

 public:
    PodioReader(const std::string& ifile);

    /** Read event and fill momentum map **/
    void event();

    /** Proton momentum **/
    inline const linal::LVect& p4p() const {
        return getp4(2212);
    }
    /** anti-Proton momentum **/
    inline const linal::LVect& p4pbar() const {
        return getp4(-2212);
    }
    /** pi+ momentum **/
    inline const linal::LVect& p4pip() const {
        return getp4(211);
    }
    /** pi- momentum **/
    inline const linal::LVect& p4pin() const {
        return getp4(-211);
    }
    /** Lambda0 momentum **/
    inline const linal::LVect& p4lam() const {
        return getp4(3122);
    }
    /** anti-Lambda0 momentum **/
    inline const linal::LVect& p4lbar() const {
        return getp4(-3122);
    }

    size_t nEvt() const {return m_nevents;}
};

#endif  // PODIOREADER_H__
