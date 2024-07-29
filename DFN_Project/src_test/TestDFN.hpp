#include <iostream>
#include <gtest/gtest.h>
#include<iomanip>
#include "Eigen/Eigen"
#include "Utils.hpp"
#include "DFNlibrary.hpp"


using namespace std;
using namespace Eigen;

namespace DFNlibrary {
const double tau = 1e-6;
//********************************
TEST(FRACTURETEST, TestImportaFratture){

    vector<Fracture> contenitoreFratture;
    int numeroFratture = 0;
    string filepath = "DFN/FR3_data.txt";
    bool testFlag = ImportaFratture(filepath, contenitoreFratture, numeroFratture);

    EXPECT_EQ(testFlag, true);
    EXPECT_EQ(numeroFratture, 3);
}

TEST(FRACTURETEST, TestCalcoloNormale){
    Fracture F;
    F.NumVertices = 4;
    F.coordx = {0,1,1,0};
    F.coordy = {0,0,1,1};
    F.coordz = {0,0,0,0};
    F.NumVertices = 4;
    F.id = 0;
    Vector3d N = CalcoloNormale(F);
    Vector3d N1 = {0,0,1};
    Vector3d N2 = {0,0,-1};
    bool testFlag = (N == N1 || N == N2);
    EXPECT_EQ(testFlag, true);
}

TEST(FRACTURETEST, TestCercaIntersezioni){
    Fracture F;
    F.NumVertices = 4;
    F.coordx = {0,1,1,0};
    F.coordy = {0,0,1,1};
    F.coordz = {0,0,0,0};
    F.id = 0;
    Vector3d P = {2,0.5,0};
    Vector3d t1 = {0,1,0};
    Vector3d t2 = {-1,0,0};
    double c1 = 0;
    double c2 = 0;
    bool testFlag1 = CercaIntersezioni(P, t1, F, c1, c2);
    bool testFlag2 = CercaIntersezioni(P, t2, F, c1, c2);
    EXPECT_EQ(testFlag1, false);
    EXPECT_EQ(testFlag2, true);
    EXPECT_EQ(c1, 1);
    EXPECT_EQ(c2, 2);
}

TEST(TRACESTEST, TestInserisciTraccia){
    double alpha0 = 4;
    double beta0 = 3;
    double gamma0 = 2;
    double delta0 = 1;
    double alpha1 = 4;
    double beta1 = 1;
    double gamma1 = 3;
    double delta1 = 2;
    double alpha5 = 4;
    double beta5 = 2;
    double gamma5 = 3;
    double delta5 = 1;

    vector<Traces> tracesContainer = {};
    Vector3d P = {0,0,0};
    Vector3d t = {1,0,0};

    int F0id1 = 10;
    int F0id2 = 11;
    int F1id1 = 0;
    int F1id2 = 1;
    int F5id1 = 2;
    int F5id2 = 3;
    InserisciTraccia(alpha0, beta0, gamma0, delta0, tracesContainer, P, t, F0id1, F0id2);

    EXPECT_EQ(tracesContainer.size(), 0);

    InserisciTraccia(alpha1, beta1, gamma1, delta1, tracesContainer, P, t, F1id1, F1id2);

    EXPECT_EQ(tracesContainer.size(), 1);
    Traces T1 = tracesContainer[0];
    EXPECT_EQ(T1.Length, 1);
    EXPECT_EQ(T1.Tips1, true);
    EXPECT_EQ(T1.Tips2, false);

    InserisciTraccia(alpha5, beta5, gamma5, delta5, tracesContainer, P, t, F5id1, F5id2);

    EXPECT_EQ(tracesContainer.size(), 2);
    Traces T2 = tracesContainer[1];
    EXPECT_EQ(T2.Length, 1);
    EXPECT_EQ(T2.Tips1, true);
    EXPECT_EQ(T2.Tips2, true);
}

TEST(FRACTURETEST, TestCercaTracce){
    Fracture F1 = {};
    F1.NumVertices = 4;
    F1.coordx = {0,1,1,0};
    F1.coordy = {0,0,1,1};
    F1.coordz = {0,0,0,0};
    F1.id = 1;

    Fracture F2 = {};
    F2.NumVertices = 4;
    F2.coordx = {2,2,2,2};
    F2.coordy = {-1,-1,1,1};
    F2.coordz = {-1,1,1,-1};
    F2.id = 2;


    Fracture F3 = {};
    F3.NumVertices = 4;
    F3.coordx = {0.5,0.5,0.5,0.5};
    F3.coordy = {-1,-1,1,1};
    F3.coordz = {-1,1,1,-1};
    F3.id = 3;

    vector<Traces> tracesContainer = {};

    CercaTracce(F1, F2, tracesContainer);
    EXPECT_EQ(tracesContainer.size(), 0);

    CercaTracce(F1, F3, tracesContainer);
    EXPECT_EQ(tracesContainer.size(), 1);
}

TEST(PRINTTEST, TestStampaTracce){
    vector<Traces> tracesContainer = {};
    Traces T = {};
    T.FractureID1 = 0;
    T.FractureID2 = 1;
    T.Length = 3;
    T.P1 = {0,0,0};
    T.P2 = {3,0,0};
    T.Tips1 = true;
    T.Tips2 = true;
    T.id = 0;
    tracesContainer.push_back(T);
    int numF = 42; // non è importante ai fine del test che tale numero non sia coerente
    StampaTracce(tracesContainer, numF);
}

TEST(SORTTEST,TestQuickSort){
    vector<int> V = {3,2,4,5,1};
    QuickSort(V);
    vector<int> SortedV = {5,4,3,2,1};
    EXPECT_EQ(V, SortedV);
}

TEST(PRINTTEST, TestStampaTracceOrdinate){
    vector<Traces> tracesContainer = {};
    Traces T = {};
    T.FractureID1 = 0;
    T.FractureID2 = 1;
    T.Length = 3;
    T.P1 = {0,0,0};
    T.P2 = {3,0,0};
    T.Tips1 = true;
    T.Tips2 = true;
    T.id = 0;
    tracesContainer.push_back(T);
    int numF = 42; // non è importante ai fine del test che tale numero non sia coerente
    map<int, vector<int>> sortedPassanti = {};
    map<int, vector<int>> sortedNonPassanti = {};
    StampaTracceOrdinate(tracesContainer, numF, sortedPassanti, sortedNonPassanti);
}

TEST(EDGETEST, TestIntersezioneEdges){
    edges L1 = {};
    L1.P = {0,0,0};
    L1.t = {6,0,0};
    L1.intersection = {};

    edges L2 = {};
    L2.P = {1,1,0};
    L2.t = {1,0,0};
    L2.intersection = {};

    edges L3 = {};
    L3.P = {3,-3,0};
    L3.t = {0,4,0};
    L3.intersection = {};

    double alpha1 = 42;
    double beta1 = 42;
    double alpha2 = 42;
    double beta2 = 42;

    IntersezioneEdges(L1, L2, alpha1, beta1);
    EXPECT_EQ(alpha1, 42);
    EXPECT_EQ(beta1, 42);

    IntersezioneEdges(L1, L3, alpha2, beta2);
    EXPECT_EQ(alpha2, 0.5);
    EXPECT_EQ(beta2, 0.75);
}


bool VerificaUguagliazaVettori(vector<double> v1, vector<double> v2){
    if(v1.size() != v2.size()){
        return false;
    }
    double difSum = 0;
    for(int i=0; i<v1.size(); i++){
        difSum += abs((v1[i]-v2[i])/v2[i]);
    }
    if(difSum < tau){
        return true;
    }
    else{
        return false;
    }
}


TEST(FRACTURETEST,TestTagliaFratture){

    Fracture F = {};
    F.NumVertices = 4;
    F.coordx = {0,1,1,0};
    F.coordy = {0,0,1,1};
    F.coordz = {0,0,0,0};

    vector<Traces> contenitoreTracce = {};
    Traces t = {};
    t.id = 0;
    t.P1 = {0.2,0,0};
    t.P2 = {0.2,1,0};
    // traccia passante
    contenitoreTracce.push_back(t);

    t.id = 1;
    t.P1 = {0.6,1,0};
    t.P2 = {0.6,0,0};
    // traccia passante
    contenitoreTracce.push_back(t);

    t.id = 2;
    t.P1 = {0.4,0,0};
    t.P2 = {0.8,1,0};
    // traccia passante
    contenitoreTracce.push_back(t);

    t.id = 3;
    t.P1 = {0.8,0,0};
    t.P2 = {0.9,0.2,0};
    // traccia non-passante
    contenitoreTracce.push_back(t);

    t.id = 4;
    t.P1 = {0.1,0.4,0};
    t.P2 = {0,0.2,0};
    // traccia non-passante
    contenitoreTracce.push_back(t);

    vector<edges> latiBordo = {};
    vector<edges> latiInterni = {};
    map<int, vector<int>> traccePassantiOrdinate {
        {0, {2,0,1}}
    };
    map<int, vector<int>> tracceNonPassantiOrdinate {
        {0, {3,4}}
    };

    TagliaFratture(F, contenitoreTracce, traccePassantiOrdinate, tracceNonPassantiOrdinate,
                                      latiBordo, latiInterni);

    vector<double> intB0 = {0.4, 0.2, 0.6, 0.8};
    vector<double> intB1 = {0.4};
    vector<double> intB2 = {0.2, 0.8, 0.4};
    vector<double> intB3 = {0.8};
    vector<double> intI0 = {0.5};
    vector<double> intI1 = {0.6};
    vector<double> intI2 = {0.5};
    vector<double> intI3 = {};
    vector<double> intI4 = {};

    EXPECT_TRUE(VerificaUguagliazaVettori(latiBordo[0].intersection, intB0));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiBordo[1].intersection, intB1));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiBordo[2].intersection, intB2));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiBordo[3].intersection, intB3));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiInterni[0].intersection, intI0));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiInterni[1].intersection, intI1));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiInterni[2].intersection, intI2));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiInterni[3].intersection, intI3));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiInterni[4].intersection, intI4));
}

TEST(FRACTURETEST,TestDoppiaEstensione){

    Fracture F = {};
    F.NumVertices = 4;
    F.coordx = {0,1,1,0};
    F.coordy = {0,0,1,1};
    F.coordz = {0,0,0,0};

    vector<Traces> contenitoreTracce = {};
    Traces t = {};
    t.id = 0;
    t.P1 = {0,0.2,0};
    t.P2 = {1,0.2,0};
    // traccia passante
    contenitoreTracce.push_back(t);

    t.id = 1;
    t.P1 = {0,0.5,0};
    t.P2 = {0.5,1,0};
    // traccia passante
    contenitoreTracce.push_back(t);

    t.id = 2;
    t.P1 = {0.2,0.5,0};
    t.P2 = {0.2,0.8,0};
    // traccia non-passante
    contenitoreTracce.push_back(t);

    vector<edges> latiBordo = {};
    vector<edges> latiInterni = {};
    map<int, vector<int>> traccePassantiOrdinate {
        {0, {0,1}}
    };
    map<int, vector<int>> tracceNonPassantiOrdinate {
        {0, {2}}
    };

    TagliaFratture(F, contenitoreTracce, traccePassantiOrdinate, tracceNonPassantiOrdinate,
                   latiBordo, latiInterni);

    vector<double> intB0 = {};
    vector<double> intB1 = {0.2};
    vector<double> intB2 = {0.5,0.8};
    vector<double> intB3 = {0.8,0.5};
    vector<double> intI0 = {0.2};
    vector<double> intI1 = {0.4};
    vector<double> intI2 = {0.375}; //  3/8

    EXPECT_TRUE(VerificaUguagliazaVettori(latiBordo[0].intersection, intB0));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiBordo[1].intersection, intB1));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiBordo[2].intersection, intB2));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiBordo[3].intersection, intB3));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiInterni[0].intersection, intI0));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiInterni[1].intersection, intI1));
    EXPECT_TRUE(VerificaUguagliazaVettori(latiInterni[2].intersection, intI2));
}

TEST(MESHTEST, TestCercaEstremo){
    vector<Vector3d> CoordinateNodi = {{0,0,0}, {2,3,4}, {5.5,6.6,7.7}};
    vector<int> idNodi = {0,1,2};
    Vector3d P1 = {7,8,9};
    Vector3d P2 = {2,3,4};
    bool flag1;
    bool flag2;
    int id1 = -1;
    int id2 = -1;

    CercaEstremo(P1, id1, CoordinateNodi, idNodi, flag1);
    CercaEstremo(P2, id2, CoordinateNodi, idNodi, flag2);

    EXPECT_EQ(flag1, false);
    EXPECT_EQ(flag2, true);
    EXPECT_EQ(id2, 1);
}


TEST(MESHTEST, TestCaricamentoCell0e1D){

    vector<edges> latiBordo = {};

    edges E = {};
    E.P = {0,0,0};
    E.t = {1,0,0};
    E.intersection = {0.4, 0.2, 0.6, 0.8};
    latiBordo.push_back(E);

    E.P = {1,0,0};
    E.t = {0,1,0};
    E.intersection = {0.4};
    latiBordo.push_back(E);

    E.P = {1,1,0};
    E.t = {-1,0,0};
    E.intersection = {0.2, 0.8, 0.4};
    latiBordo.push_back(E);

    E.P = {0,1,0};
    E.t = {0,-1,0};
    E.intersection = {0.8};
    latiBordo.push_back(E);

    vector<edges> latiInterni = {};

    E.P = {0.4,0,0};
    E.t = {0.4,1,0};
    E.intersection = {0.5};
    latiInterni.push_back(E);

    E.P = {0.2,0,0};
    E.t = {0,1,0};
    E.intersection = {0.6};
    latiInterni.push_back(E);

    E.P = {0.6,1,0};
    E.t = {0,-1,0};
    E.intersection = {0.5};
    latiInterni.push_back(E);

    E.P = {0.8,0,0};
    E.t = {0.2,0.4,0};
    E.intersection = {};
    latiInterni.push_back(E);

    E.P = {0,0.2,0};
    E.t = {0.2,0.4,0};
    E.intersection = {};
    latiInterni.push_back(E);

    PolygonalMesh mesh = {};
    vector<int> idBordo = {};
    vector<int> idInterno = {};
    CaricamentoCell0e1D(latiBordo, latiInterni, mesh, idBordo, idInterno);

    EXPECT_EQ(mesh.NumberCell0D, 15);
    EXPECT_EQ(mesh.NumberCell1D, 21);
}

TEST(MESHTEST, TestCalcoloNormaleMesh){
    PolygonalMesh mesh = {};
    mesh.NumberCell0D = 4;
    mesh.Cell0DCoordinates = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}};
    Vector3d N = CalcoloNormaleMesh(mesh);
    Vector3d N1 = {0,0,1};
    Vector3d N2 = {0,0,-1};
    bool testFlag = (N == N1 || N == N2);
    EXPECT_EQ(testFlag, true);
}

TEST(MESHTEST, TestLatoSuccessivo){
    PolygonalMesh mesh = {};
    mesh.NumberCell0D = 15;
    mesh.Cell0DId = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    mesh.Cell0DCoordinates = {{0,0,0},{0.2,0,0},{0.4,0,0},{0.6,0,0},{0.7,0,0},
                              {1,0,0},{1,0.4,0},{1,1,0},
                              {0.8,1,0},{0.6,1,0},{0.2,1,0},
                              {0,1,0},{0,0.2,0},
                              {0.6,0.5,0}, {0.2,0.6,0}};
    mesh.InvolvedEdges = {{0,12},{0,1,15},{1,2,13},{2,3,18},{3,4,19},{4,5},{5,6,19},{6,7},
                          {7,8,14},{8,9,17},{9,10,16},{10,11},{11,12,20},{13,14,17,18},{15,16,20}};

    mesh.NumberCell1D = 21;
    mesh.Cell1DId = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    mesh.Cell1DVertices = {{0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},{7,8},{8,9},{9,10},{10,11},{11,12},{12,0},
                           {2,13},{13,8},{1,14},{14,10},{9,13},{13,3},{4,6},{12,14}};

    Vector3d CurrentEdgeTan = {0.2,0.5,0};
    int CurrentNode = 13;
    int CurrentEdgeId = 13;
    vector<int> nodiPoly = {2,13};
    vector<int> latiPoly = {13};
    int inverti = 2;
    bool chiuso = false;
    Vector3d N = {0,0,1};
    vector<int> idBordo = {0,1,2,3,4,5,6,7,8,9,10,11,12};
    vector<int> idInterno = {13,14,15,16,17,18,19,20};
    LatoSuccessivo(CurrentEdgeTan, CurrentNode, CurrentEdgeId,
                   nodiPoly, latiPoly, inverti, chiuso,
                   mesh, N, idBordo, idInterno);

    Vector3d NuovaTangente = {0,0.5,0};

    EXPECT_EQ(CurrentEdgeTan, NuovaTangente);
    EXPECT_EQ(CurrentNode, 9);
    EXPECT_EQ(CurrentEdgeId, 17);
    EXPECT_EQ(inverti, 2);
    EXPECT_EQ(chiuso, false);
}


TEST(MESHTEST, TestCaricamentoCell2D){
    PolygonalMesh mesh = {};
    mesh.NumberCell0D = 15;
    mesh.Cell0DId = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    mesh.Cell0DCoordinates = {{0,0,0},{0.2,0,0},{0.4,0,0},{0.6,0,0},{0.7,0,0},
                              {1,0,0},{1,0.4,0},{1,1,0},
                              {0.8,1,0},{0.6,1,0},{0.2,1,0},
                              {0,1,0},{0,0.2,0},
                              {0.6,0.5,0}, {0.2,0.6,0}};
    mesh.InvolvedEdges = {{0,12},{0,1,15},{1,2,13},{2,3,18},{3,4,19},{4,5},{5,6,19},{6,7},
                          {7,8,14},{8,9,17},{9,10,16},{10,11},{11,12,20},{13,14,17,18},{15,16,20}};

    mesh.NumberCell1D = 21;
    mesh.Cell1DId = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    mesh.Cell1DVertices = {{0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},{7,8},{8,9},{9,10},{10,11},{11,12},{12,0},
                           {2,13},{13,8},{1,14},{14,10},{9,13},{13,3},{4,6},{12,14}};

    vector<int> idBordo = {0,1,2,3,4,5,6,7,8,9,10,11,12};
    vector<int> idInterno = {13,14,15,16,17,18,19,20};
    CaricamentoCell2D (mesh, idBordo, idInterno);

    vector<int> Cell2D1Vertices = {0,1,14,12};
    vector<int> Cell2D1Edges = {0,15,20,12};

    EXPECT_EQ(mesh.NumberCell2D, 7);
    EXPECT_EQ(mesh.Cell2DVertices[0], Cell2D1Vertices);
    EXPECT_EQ(mesh.Cell2DEdges[0], Cell2D1Edges);

}


}
