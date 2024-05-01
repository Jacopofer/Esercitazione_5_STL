#include <iostream>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

int main()
{
    PolygonalMesh mesh;
    string filepath = "PolygonalMesh";

    if (!ImportMesh(filepath, mesh))
    {
        return 1;
    }

    return 0;
}
