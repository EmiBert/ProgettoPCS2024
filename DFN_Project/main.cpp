#include <iostream>
#include "Utils.hpp"
#include "DFNlibrary.hpp"
#include "Eigen/Eigen"
#include<iomanip>

using namespace std;
using namespace Eigen;
using namespace DFNlibrary;

int main()
{
    vector<Fracture> contenitoreFratture;
    int numeroFratture = 0;
    string filepath = "DFN/FR10_data.txt";

    if(!ImportaFratture(filepath, contenitoreFratture, numeroFratture)){
        return 1;
    }

    vector<Traces> contenitoreTracce = {};

    for(int i =0; i<numeroFratture-1; i++){
        for(int j =i+1; j<numeroFratture; j++){
            CercaTracce(contenitoreFratture[i],contenitoreFratture[j],contenitoreTracce);
        }
    }

/*
    cout<<" # Number of Fractures"<<endl;
    cout<<numeroFratture<<endl;

    for(auto f: contenitoreFratture){
        cout<<"# FractureId; NumVertices"<<endl;
        cout<<f.id<<"; "<<f.NumVertices<<endl;
        for(auto x: f.coordx){
            cout<<x<<"; ";
        }
        cout<<endl;
        for(auto y: f.coordy){
            cout<<y<<"; ";
        }
        cout<<endl;
        for(auto z: f.coordz){
            cout<<z<<"; ";
        }
        cout<<endl;
    }




    cout<<" # TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2; Tips1; Tips2; length"<<endl;
    for(auto t: contenitoreTracce){
        cout<< t.id<<";";
        cout<< t.FractureID1<<";";
        cout<< t.FractureID2<<";";
        cout<< t.P1[0]<<";";
        cout<< t.P1[1]<<";";
        cout<< t.P1[2]<<";";
        cout<< t.P2[0]<<";";
        cout<< t.P2[1]<<";";
        cout<< t.P2[2]<<";";
        cout<< t.Tips1<<";";
        cout<< t.Tips2<<";";
        cout<< t.Length<<endl;
    }
*/


    StampaTracce(contenitoreTracce,numeroFratture);

    map<int, vector<int>> traccePassantiOrdinate = {};
    map<int, vector<int>> tracceNonPassantiOrdinate = {};

    StampaTracceOrdinate(contenitoreTracce,numeroFratture,
                         traccePassantiOrdinate,tracceNonPassantiOrdinate);


/*
    for(int i =0; i<numeroFratture; i++){
        if(traccePassantiOrdinate.find(i) != traccePassantiOrdinate.end()){
            cout<<"tracce frattura "<<i<<" passanti"<<endl;
            for(auto& e: traccePassantiOrdinate[i]){
                cout<<e<<endl;
            }
        }
        if(tracceNonPassantiOrdinate.find(i) != tracceNonPassantiOrdinate.end()){
            cout<<"tracce frattura "<<i<<" non passanti"<<endl;
            for(auto& e: tracceNonPassantiOrdinate[i]){
                cout<<e<<endl;
            }
        }
    }
*/
    int indiceFrattura = 0;
    vector<edges> latiBordo = {};
    vector<edges> latiInterni = {};

    TagliaFratture(contenitoreFratture[indiceFrattura], contenitoreTracce,
                   traccePassantiOrdinate, tracceNonPassantiOrdinate,
                   latiBordo,latiInterni);

    //stampa bordo e interni
    cout<<"lati bordo"<<endl;
    for(auto& e: latiBordo){
        cout<<"P: "<<endl<<e.P<<endl;
        cout<<"t: "<<endl<<e.t<<endl;
        cout<<"intersection: ";
        for(auto& i: e.intersection){
            cout<<setprecision(8)<<i<<";";
        }
        cout<<endl;
    }
    cout<<endl<<"------------------------"<<endl;
    cout<<"lati interni"<<endl;
    for(auto& e: latiInterni){
        cout<<"P: "<<endl<<e.P<<endl;
        cout<<"t: "<<endl<<e.t<<endl;
        cout<<"intersection: ";
        for(auto& i: e.intersection){
            cout<<setprecision(8)<<i<<";";
        }
        cout<<endl;
    }

    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    PolygonalMesh mesh = {};
    vector<int> idBordo = {};
    vector<int> idInterno = {};
    CaricamentoCell0e1D(latiBordo,latiInterni, mesh, idBordo,idInterno);

    cout<<"NumberCell0D: "<<mesh.NumberCell0D<<endl;
    cout<<"Cell0id: "<<endl;
    for(auto& e: mesh.Cell0DId){
        cout<<e<<endl;
    }
    cout<<"Cell0DCoordinates: "<<endl;
    for(auto& e: mesh.Cell0DCoordinates){
        cout<<e[0]<<"  "<<e[1]<<"  "<<e[2]<<endl;
    }
    cout<<"InvolvedEdges: "<<endl;
    for(auto& e: mesh.InvolvedEdges){
        cout<<"*";
        for(auto& i: e){
            cout<<i<<"  ";
        }
        cout<<endl;
    }

    cout<<"NumberCell1D: "<<mesh.NumberCell1D<<endl;
    cout<<"Cell1id: "<<endl;
    for(auto& e: mesh.Cell1DId){
        cout<<e<<endl;
    }
    cout<<"Cell1DVertices: "<<endl;
    for(auto& e: mesh.Cell1DVertices){
        cout<<e[0]<<"  "<<e[1]<<endl;
    }



    CaricamentoCell2D (mesh, idBordo, idInterno);

    cout<<"*****************************************"<<endl;
    cout<<"mesh.Cell2DId:"<<endl;
    for (auto& x : mesh.Cell2DId){
        cout<<"id: "<<x.first<<"  dim: "<<x.second<<endl;
    }
    cout<<"mesh.Cell2DVertices:"<<endl;
    for(auto& c: mesh.Cell2DVertices){
        for(auto& v: c){
            cout<<v<<" ";
        }
        cout<<endl;
    }
    cout<<"mesh.Cell2DEdges:"<<endl;
    for(auto& c: mesh.Cell2DEdges){
        for(auto& e: c){
            cout<<e<<" ";
        }
        cout<<endl;
    }
    cout<<"yyyyyyyyyyyyyyyyyyyyyyyyyyyy"<<endl;
    cout<<"mesh.NumberCell2D: "<<mesh.NumberCell2D<<endl;
return 0;
}
