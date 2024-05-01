#pragma once

#include <iostream>
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolygonalLibrary {

bool ImportCell0Ds(const string &filename, PolygonalMesh& mesh);

bool ImportCell1Ds(const string &filename, PolygonalMesh& mesh);

bool ImportCell2Ds(const string &filename, PolygonalMesh& mesh);

bool ImportMesh(const string &filepath, PolygonalMesh& mesh);
}
