#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <Eigen/Eigen>
#include <cmath>

namespace PolygonalLibrary {

bool ImportCell0Ds(const string &filename, PolygonalMesh& mesh)
{
    ifstream file;
    file.open(filename);

    if (file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
    {
        replace(line.begin(), line.end(), ';', ' ');
        listLines.push_back(line);
    }

    listLines.pop_front();

    mesh.NumberCell0D = listLines.size();
    if (mesh.NumberCell0D == 0)
    {
        cerr << "Non ci sono celle 0D" << endl;
        return false;
    }

    mesh.Cell0DId.reserve(mesh.NumberCell0D);
    mesh.Cell0DCoordinates.reserve(mesh.NumberCell0D);

    for (const string& line : listLines)
    {
        istringstream converter(line);
        unsigned int id;
        unsigned int marker;
        Vector2d coord;

        converter >> id >> marker >> coord(0) >> coord(1);

        mesh.Cell0DId.push_back(id);
        mesh.Cell0DCoordinates.push_back(coord);

        if (marker != 0)
        {
            auto ret = mesh.Cell0DMarkers.insert({marker, {id}});
            if (!ret.second)
                (ret.first)->second.push_back(id);
        }
    }

    file.close();
    return true;
}



bool ImportCell1Ds(const string &filename, PolygonalMesh& mesh)
{
    ifstream file;
    file.open(filename);

    if (file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
    {
        replace(line.begin(), line.end(), ';', ' ');
        listLines.push_back(line);
    }

    listLines.pop_front();

    mesh.NumberCell1D = listLines.size();

    if (mesh.NumberCell1D == 0)
    {
        cerr << "Non ci sono celle 1D" << endl;
        return false;
    }

    mesh.Cell1DId.reserve(mesh.NumberCell1D);
    mesh.Cell1DVertices.reserve(mesh.NumberCell1D);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2i vertices;

        converter >> id >> marker >> vertices(0) >> vertices(1);

        mesh.Cell1DId.push_back(id);
        mesh.Cell1DVertices.push_back(vertices);

        if (marker != 0)
        {
            auto ret = mesh.Cell1DMarkers.insert({marker, {id}});
            if (!ret.second)
                (ret.first)->second.push_back(id);
        }
    }

    file.close();
    return true;
}



bool ImportCell2Ds(const string &filename, PolygonalMesh& mesh)
{
    ifstream file;
    file.open(filename);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
    {
        replace(line.begin(), line.end(), ';', ' ');
        listLines.push_back(line);
    }

    listLines.pop_front();

    mesh.NumberCell2D = listLines.size();

    if (mesh.NumberCell2D == 0)
    {
        cerr << "Non ci sono celle 2D" << endl;
        return false;
    }

    mesh.Cell2DId.reserve(mesh.NumberCell2D);
    mesh.Cell2DVertices.reserve(mesh.NumberCell2D);
    mesh.Cell2DEdges.reserve(mesh.NumberCell2D);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        converter >> id;
        converter >> marker;

        unsigned int NumVertices;
        converter >> NumVertices;

        vector<unsigned int> vertices(NumVertices);
        for (unsigned int i = 0; i < NumVertices; i++)
            converter >> vertices[i];

        unsigned int NumEdges;
        converter >> NumEdges;

        vector<unsigned int> edges(NumEdges);
        for (unsigned int i = 0; i < NumEdges; i++)
            converter >> edges[i];

        mesh.Cell2DId.push_back(id);
        mesh.Cell2DVertices.push_back(vertices);
        mesh.Cell2DEdges.push_back(edges);
    }

    file.close();
    return true;
}


bool ImportMesh(const string &filepath, PolygonalMesh& mesh)
{

    if (!ImportCell0Ds(filepath + "/Cell0Ds.csv", mesh))
    {
        return false;
    }
    else
    {
        cout << "Cell0D marker:" << endl;
        for (auto it = mesh.Cell0DMarkers.begin(); it != mesh.Cell0DMarkers.end(); it++)
        {
            cout << "key:\t" << it->first << "\t values:";
            for (const unsigned int id : it->second)
                cout << "\t" << id;

            cout << endl;
        }
    }


    if (!ImportCell1Ds(filepath + "/Cell1Ds.csv", mesh))
    {
        return false;
    }
    else
    {
        cout << "Cell1D marker:" << endl;
        for (auto it = mesh.Cell1DMarkers.begin(); it != mesh.Cell1DMarkers.end(); it++)
        {
            cout << "key:\t" << it->first << "\t values:";
            for (const unsigned int id : it->second)
                cout << "\t" << id;

            cout << endl;
        }

        // TEST LUNGHEZZA LATI
        for (unsigned int c = 0; c < mesh.NumberCell1D; c++)
        {
            Vector2i vertices = mesh.Cell1DVertices[c]; // prendo il vettore che contiene i vertici del c-esimo lato
            // primo vertice: coordinate
            double x1 = mesh.Cell0DCoordinates[vertices[0]][0];
            double y1 = mesh.Cell0DCoordinates[vertices[0]][1];
            // secondo vertice
            double x2 = mesh.Cell0DCoordinates[vertices[1]][0];
            double y2 = mesh.Cell0DCoordinates[vertices[1]][1];
            // differenza tra le x e le y
            double dx = x2 - x1;
            double dy = y2 - y1;

            double dist = sqrt(dx*dx + dy*dy);

            if (dist == 0)
            {
                cerr << "Sono stati trovati poligoni con lati degeneri (lunghezza nulla)" << endl;
                return 2;
            }
        }
    }

    if (!ImportCell2Ds(filepath + "/Cell2Ds.csv", mesh))
    {
        return false;
    }
    else
    {
        // TEST SU AREA POLIGONI
        for (unsigned int c = 0; c < mesh.NumberCell2D; c++)
        {
            vector<unsigned int> vertices = mesh.Cell2DVertices[c]; // vettore che contiene i vertici del c-esimo poligono
            unsigned int n = vertices.size();
            double area = 0.0;

            for (unsigned int i=0; i<n; i++)
            {
                unsigned int j = (i+1)%n;
                area += mesh.Cell0DCoordinates[vertices[i]][0] * mesh.Cell0DCoordinates[vertices[j]][1];
                area -= mesh.Cell0DCoordinates[vertices[j]][0] * mesh.Cell0DCoordinates[vertices[i]][1];
            }

            area = abs(area)/2.0;

            if (area == 0)
            {
                cerr << "Sono stati trovati poligoni degeneri (area nulla)" << endl;
                return 3;
            }
        }
    }

    return true;
}

}
