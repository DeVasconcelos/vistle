#ifndef VISTLE_READENSIGHT_GEOGOLDBIN_H
#define VISTLE_READENSIGHT_GEOGOLDBIN_H

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CLASS    GeoGoldBin
//
// Description: Abstraction of EnSight Gold geometry Files
//
// Initial version: 08.08.2003
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// (C) 2002 / 2003 by VirCinity IT Consulting
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//  Changes:
//

#include "EnFile.h"

#include <string>

class GeoGoldBin: public EnFile {
public:
    // creates file-rep. and opens the file
    GeoGoldBin(ReadEnsight *mod, const std::string &name, CaseFile::BinType binType = CaseFile::CBIN);
    ~GeoGoldBin();

    // read the file
    vistle::Object::ptr read(int timestep, int block, EnPart *part) override;

    // get part info
    bool parseForParts() override;

private:
    // read header
    bool readHeader(FILE *in);

    // read EnSight part information
    bool readPart(FILE *in, EnPart &actPart);

    // read part connectivities (EnSight Gold only)
    bool readPartConn(FILE *in, EnPart &actPart);

    // read bounding box (EnSight Gold)
    bool readBB(FILE *in);

    size_t m_numCoords = 0; // number of coordinates
    int m_actPartNum = -1;
    bool partFound = false;
};
#endif
