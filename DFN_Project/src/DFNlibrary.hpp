#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace DFNlibrary {

    //int NumFracture = 0; ///< number of fracture
struct Fracture
{
    int id = 0; ///< numero identificativo della frattura
    int NumVertices = 0; ///< NumVertices = numero di vertici
    vector<double> coordx = {}; ///< coordx[i] = coordinata x del vertice i
    vector<double> coordy = {}; ///< coordy[i] = coordinata y del vertice i
    vector<double> coordz = {}; ///< coordz[i] = coordinata z del vertice i
};




struct Traces
{   int id = 0; ///< numero identificativo della traccia
    int FractureID1 = -1; ///< numero identificativo della frattura 1
    int FractureID2 = -1; ///< numero identificativo della frattura 2
    Vector3d P1 = {}; ///< P1 = estremo 1
    Vector3d P2 = {}; ///< P2 = estremo 2
    bool Tips1 = {}; ///< Tips1 falso se la traccia è passante rispetto alla frattura ID1, vero se è non-passante
    bool Tips2 = {}; ///< Tips2 falso se la traccia è passante rispetto alla frattura ID2, vero se è non-passante
    double Length = {}; ///< lunghezza della traccia
};

// struttura utilizzata per ordinare le tracce con quicksort
// restringe le informazioni della traccia ai soli parametri
// id e lenght, inoltre contiene una modalità di inizializzazione
// ed un operatore di confronto
struct ShortTraces{
    int id = -1;
    double Length = -1;

    ShortTraces(int a,double b)    // inizializzazione tramite passaggio di parametri
    {   id = a;
        Length = b;}

    friend bool operator >=(const ShortTraces& s1, const ShortTraces& s2)
    {
        bool R;
        if (s1.Length>=s2.Length){
            R = true;
        }
        else{
            R = false;
        }
        return R;
    }
};



struct edges
{
    Vector3d P = {}; ///< P = origine del segmento
    Vector3d t = {}; ///< t = tangente del segmento (non normalizzata)
    vector<double> intersection = {}; ///< intersection = insiemi di valori apparteneti a (0,1)
                                      ///< in cui il segnmento è stato tagliato
};



struct PolygonalMesh
{
    int NumberCell0D = 0; ///< numero di Cell0D, ovvero nodi
    vector<int> Cell0DId = {}; ///< Cell0D id, size 1 x NumberCell0D
    vector<Vector3d> Cell0DCoordinates = {}; ///< coordinate Cell0D, size 3 x NumberCell0D
    vector<vector<int>> InvolvedEdges = {}; ///< InvolvedEdges[i] = lista contenente Cell1D id
                                            ///< appartenenti ad ogni lato avente il nodo
                                            ///< Cell0DId[i] come estremo

    int NumberCell1D = 0; ///< numero di Cell1D
    vector<int> Cell1DId = {}; ///< Cell1D id, size 1 x NumberCell1D
    vector<Vector2i> Cell1DVertices = {}; ///< indici dei vertici delle Cell1D,
                                          ///< size 2 x NumberCell1D (fromId,toId)

    int NumberCell2D = 0; ///< numero di Cell2D
    map<int,int> Cell2DId = {}; ///< codice identificativo e dimensione (numero lati)
                                ///< delle Cell2D {id,dimensione}; size 2 x NumberCell2D
    /// note: dimensione = NumVertcies = NumEdges (quando sono uguali, altrimenti dimensione = 0)

    vector<vector<int>> Cell2DVertices = {}; ///< insieme degli indici dei vertici delle Cell2D,
                                             ///< Cell2DVertices[i] si riferisce alla cella di
                                             ///< indice i (non all'id, che potrebbe essere diverso)
                                             ///< (la dimensione è dinamica, non è una tabella "n x m")

    vector<vector<int>> Cell2DEdges = {}; ///< insieme degli indici dei lati delle Cell2D,
                                          ///< Cell2DEdges[i] si riferisce alla cella di
                                          ///< indice i (non all'id, che potrebbe essere diverso)
                                          ///< (la dimensione è dinamica, non è una tabella "n x m")
};

}
