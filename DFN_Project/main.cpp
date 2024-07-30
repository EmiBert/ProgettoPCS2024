#include <iostream>
#include "Utils.hpp"
#include "DFNlibrary.hpp"
#include "Eigen/Eigen"

// progetto di Emilio Enea Bertello, 295579
// unico componente del gruppo

using namespace std;
using namespace Eigen;
using namespace DFNlibrary;

int main()
{
    vector<Fracture> contenitoreFratture;
    int numeroFratture = 0;
    string filepath = "DFN/FR3_data.txt";

    if(!ImportaFratture(filepath, contenitoreFratture, numeroFratture)){
        return 1;
    }

    vector<Traces> contenitoreTracce = {};

    for(int i =0; i<numeroFratture-1; i++){
        for(int j =i+1; j<numeroFratture; j++){
            CercaTracce(contenitoreFratture[i],contenitoreFratture[j],contenitoreTracce);
        }
    }

    StampaTracce(contenitoreTracce,numeroFratture);

    map<int, vector<int>> traccePassantiOrdinate = {};
    map<int, vector<int>> tracceNonPassantiOrdinate = {};

    StampaTracceOrdinate(contenitoreTracce,numeroFratture,
                         traccePassantiOrdinate,tracceNonPassantiOrdinate);


    vector<PolygonalMesh> contenitoreMesh = {};

    // caricamento mesh
    for(int i = 0; i<numeroFratture; i++){
        vector<edges> latiBordo = {};
        vector<edges> latiInterni = {};

        TagliaFratture(contenitoreFratture[i], contenitoreTracce,
                       traccePassantiOrdinate, tracceNonPassantiOrdinate,
                       latiBordo,latiInterni);

        PolygonalMesh mesh = {};
        vector<int> idBordo = {};
        vector<int> idInterno = {};
        CaricamentoCell0e1D(latiBordo,latiInterni, mesh, idBordo, idInterno);
        CaricamentoCell2D (mesh, idBordo, idInterno);
        contenitoreMesh.push_back(mesh);
    }

    int idF = 0;
    cout<<"elaborazione mesh dal file: "<<filepath<<endl;
    for(auto& m: contenitoreMesh){
        cout<<"----Frattura numero "<<idF<<"----"<<endl<<endl;
        cout<<"numero nodi: "<<m.NumberCell0D<<endl;
        cout<<"Nodi: id; coordinate; id archi uscenti:"<<endl;
        for(int indiceN = 0; indiceN<m.NumberCell0D; indiceN++){
            cout<<m.Cell0DId[indiceN]<<"; ";
            cout<<m.Cell0DCoordinates[indiceN][0]<<"; "<<m.Cell0DCoordinates[indiceN][1]<<"; "<<m.Cell0DCoordinates[indiceN][2];
            for(auto& invE: m.InvolvedEdges[indiceN]){
                cout<<"; "<<invE;
            }
            cout<<endl;
        }
        cout<<endl;

        cout<<"numero archi: "<<m.NumberCell1D<<endl;
        cout<<"Archi: id; vertice1; vertice2 "<<endl;
        for(int indiceA = 0; indiceA<m.NumberCell1D; indiceA++){
            cout<<m.Cell1DId[indiceA]<<"; ";
            cout<<m.Cell1DVertices[indiceA][0]<<"; "<<m.Cell1DVertices[indiceA][1]<<endl;
        }
        cout<<endl;

        cout<<"numero celle: "<<m.NumberCell2D<<endl<<endl;
        for(int indiceC = 0; indiceC<m.NumberCell2D; indiceC++){
            cout<<"Celle: id; dimensione"<<endl;
            cout<<indiceC<<"; "<<m.Cell2DId[indiceC]<<endl;
            cout<<"id vertici: ";
            for(auto& v: m.Cell2DVertices[indiceC]){
                cout<< v<<"; ";
            }
            cout<<endl;
            cout<<"id lati: ";
            for(auto& l: m.Cell2DEdges[indiceC]){
                cout<< l<<"; ";
            }
            cout<<endl<<endl;
        }
        cout<<endl;
        idF++;
    }

return 0;
}
