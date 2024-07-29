#pragma once

#include <iostream>
#include "DFNlibrary.hpp"

using namespace std;
namespace DFNlibrary{

// Funzione che legge e memorizza dal file di input
// le informazioni sulle fratture considerate
//
// filepath = percorso del file di input da leggere
// farct = vettore di Fratture in cui vengono memorizzare le informazioni
// numF = numero di fratture complessivo
//
// return = true se il processo di lettura e caricamento va a
// buon fine, false altrimenti (se il file non si apre o se
// non sono presenti fratture)
bool ImportaFratture(const string& filepath,
                     vector<Fracture>& fract,
                     int& numF);


// Funzione che calcola la normale al piano
// contenete la frattura considerata
//
// F = frattura considerata
//
// return = vettore normale
Vector3d CalcoloNormale(const Fracture F);


// Funzione che cerca le intersezioni tra una frattura ed
// una retta ad essa complanare. Se vi è intersezione
// vengono memorizzati i due estremi del segmento di intersezione
// (si suppongono fratture convesse e che non avvengano mai
// intersezioni degeneri con un solo punto in comune)
//
// P = punto appartenente alla retta complanare alla frattura
// t = tangente della retta complanare alla frattura
// F = la frattura considerata
// c1 = intersezione 1 (a valori in (0,1), il punto si
// ricava come P + c1 * t)
// c2 = intersezione 2 (a valori in (0,1), il punto si
// ricava come P + c2 * t)
//
// return = true se c'è intersezione, false altrimenti
bool CercaIntersezioni(Vector3d P,
                       Vector3d t,
                       Fracture F,
                       double& c1,
                       double& c2);


// Conoscendo le posizioni di due segmenti sulla stessa retta
// la funzione ne valuta la posizione reciproca ed aggiunge
// le informazioni della traccia giacente sulla retta (relativa
// alle due fratture in esame) al contenitore delle tracce.
//
// alpha = coefficiente 1 del segmento 1
// beta = coefficiente 2 del segmento 1
// gamma = coefficiente 1 del segmento 2
// delta = coefficiente 2 del segmento 2
// tracesContainer = contenitore delle tracce
// P = punto della retta
// t = tangente della retta
// Fid1 = frattura id 1
// Fid2 = frattura id 2
void InserisciTraccia(double alpha,
                      double beta,
                      double gamma,
                      double delta,
                      vector<Traces>& tracesContainer,
                      Vector3d P,
                      Vector3d t,
                      int Fid1,
                      int Fid2);



// Funzione che valuta se due fratture si interecano
// ed eventualmente ne calcola la traccia e la aggiunge
// al contenitore "tracesContainer"
//
// F1 = frattura 1
// F2 = frattura 2
// tracesContainer = contenitore delle tracce
void CercaTracce(const Fracture F1,
                 const Fracture F2,
                 vector<Traces>& tracesContainer);



// Funzione che stampa in un file di output denominato
// "traces_FRX_data.csv", con X = numero fratture nel file considerato
// tutte le tracce presenti.
//
// tracesContainer = contenitore delle tracce
// numF = numero di fratture
void StampaTracce(vector<Traces> tracesContainer,
                  int numF);



// Le funzioni "Scambia", "Distribuzione" e "Quicksort"
// sono definite e descritte in questo file hpp  e non
// in quello cpp per evitare l'errore di molteplice
// definizione della stessa funzione.
// Queste funzioni corrispondono a quelle classiche
// utilizzate nell'implementare l'algoritmo di quicksort
//
// A = vettore di variabili T (template) da ordinare
//
// nota: il vettore sarà ordinato in senso decrescente
template<typename T>
void Scambia(vector<T>& A, int i, int j){
    T temp = A[j];
    A[j] = A[i];
    A[i] = temp;
}


template<typename T>
int Distribuzione(vector<T>& A, int sinistra, int destra){
    T x = A[destra];
    int i = sinistra-1;
    for(int j=sinistra; j<destra; j++){
        if(A[j] >= x){
            i++;
            Scambia(A,i,j);
        }
    }
    Scambia(A,i+1,destra);
    return i+1;
}



//    2 pre: 0≤sinistra,destra≤n−1
template<typename T>
void QuickSort(vector<T>& A, int sinistra, int destra){
    if (sinistra < destra){
        // il pivot è l'ultimo indice, "destra"
        int rango = Distribuzione(A, sinistra, destra);
        QuickSort(A, sinistra, rango-1);
        QuickSort(A, rango+1, destra);
    }
}


template<typename T>
void QuickSort(vector<T>& A){
    int sinistra = 0;
    int destra = A.size()-1;
    QuickSort(A, sinistra, destra);
}

// Funzione che stampa in un file di output
// l'elenco delle tracce ordinate in senso decrescente
// (e divise precedentemente in passanti e non-passanti).
// Inoltre carica le mappe "sortedPassanti" e
// "sortedNonPassanti".
//
// tracesContainer = contenitore delle tracce
// numFracture = numero di fratture
// sortedPassanti = mappa avente come chiave l'id della frattura
// e come valore l'elenco delle tracce passanti per quella frattura
//
// sortedNoNPassanti = mappa avente come chiave l'id della frattura
// e come valore l'elenco delle tracce non-passanti per quella frattura
void StampaTracceOrdinate(vector<Traces> tracesContainer,
                          int numFracture,
                          map<int, vector<int>>& sortedPassanti,
                          map<int, vector<int>>& sortedNonPassanti);


// Funzione che controlla se le due rette su cui giacciono
// i due lati (complanari) L1 e L2 si intersecano, e nel
// caso ne calcola l'intersezione.
// nota: anche se le due rette si intersecano non è detto
// che lo facciano pure i lati.
//
// L1 = primo lato
// L2 = secondo lato
// alpha = coefficiente intersezione primo lato
// beta = coefficiente intersezione secondo lato
void IntersezioneEdges(edges L1,
                       edges L2,
                       double& alpha,
                       double& beta);

// Funzione che a partire dall'informazione sulla frattura
// considerata e dalla tracce associate calcola tutti i lati
// (interni ed esterni) di tale frattura con le relative intersezioni.
// Inoltre per quanto riguarda le tracce non-passanti viene
// applicato il processo di estensione di tale tracce, seguendo
// l'ordine decrescente sulle lunghezze.
//
// F = frattura considerata
// contenitoreTracce = contenitore di tutte le tracce
// tPO = mappa delle tracce passanti ordinate, divise per frattura
// tNPO = mappa delle tracce non-passanti ordinate, divise per frattura
// latiBordo = insieme dei lati al bordo della frattura
// latiInterni  = insieme dei lati interni della frattura
// nota: con "lati" non si intendono le celle 1D ma le tracce
// con le proprie intersezioni
void TagliaFratture(Fracture F,
                    vector<Traces> contenitoreTracce,
                    map<int, vector<int>> tPO,
                    map<int, vector<int>> tNPO,
                    vector<edges>& latiBordo,
                    vector<edges>& latiInterni);

// Funzione che dato un punto dello spazio P valuta
// se tale punto è presente nel contenitore "CoordinateNodi"
// e in tal caso imposta flag = true e salva l'indice
// del nodo nel contenitore in id.
//
// P = punto da cercare
// id = id del nodo pari a P
// CoordinateNodi = contenitore dei nodi
// idNodi = identificativi dei nodi
// flag = booleano pari a true se il nodo è stato trovato, false altrimenti
//
// nota: la funzione è di tipo "void" e non "bool" perchè
// quando si vuole cercare un punto che sappiamo essere presente
// tra l'insieme di nodi non siamo interessati al valore
// della flag (che sarebbe sempre true)
void CercaEstremo(Vector3d P,
                  int& id,
                  vector<Vector3d> CoordinateNodi,
                  vector<int> idNodi,
                  bool& flag);


// Versione di CercaEstremo senza la flag tra gli argomenti
void CercaEstremo(Vector3d P,
                  int& id,
                  vector<Vector3d> CoordinateNodi,
                  vector<int> idNodi);

// Funzione che a partire dalle informazioni sui
// lati del bordo e interni carica nella mesh
// tutte le informazioni sui nodi e sugli archi
// (rispettivamente Cell0D e Cell1D).
// Inoltre sono caricati gli elenchi degli id degli
// archi al bordo e degli archi interni all'interno
// di appositi contenitori.
//
// latiBordo = contenitore dei lati al bordo della frattura
// latiInterni= contenitore dei lati interni della frattura
// mesh = struttura sui cui vengono caricate tutte le informazioni
// relative alla frattura tagliata
// idBordo = elenco degli id degli archi al bordo
// idInterno = elenco degli id degli archi interni
void CaricamentoCell0e1D(vector<edges>& latiBordo,
                         vector<edges>& latiInterni,
                         PolygonalMesh& mesh,
                         vector<int>& idBordo,
                         vector<int>& idInterno);



// Funzione che data una mesh ne calcola il vettore normale
//
// mesh = mesh considerata
//
// return = vettore normale
Vector3d CalcoloNormaleMesh(PolygonalMesh mesh);


// !!!ATTENZIONE!!!
// da questo momento in avanti ci si riferirà, con
// un abuso di notazione, agli archi con il
// termine "lati"




// Funzione che a partire dal un lato del poligono
// in costruzione ricerca il lato successivo, e
// aggiorna i dati del poligono, segnalandone
// l'eventuale chiusura (completamento).
//
// CurrentEdgeTan = tangente del lato corrente
// CurrentNode = nodo corrente (da cui escono i
// possibili lati successivi)
// CurrentEdgeId = id del lato corrente
// nodiPoly = elenco dei nodi del poligono in costruzione
// latiPoly = elenco dei lati del poligono in costruzione
// inverti = valore che indica la condizione del poligono
// chiuso = flag che indica se il poligono è completo o no
// mesh = mesh a cui appartiene il poligono in costruzione
// N = normale della mesh
// idBordo = elenco degli id dei lati al bordo
// idInterno = elenco degli id dei lati interni
void LatoSuccessivo(Vector3d& CurrentEdgeTan,
                    int& CurrentNode,
                    int& CurrentEdgeId,
                    vector<int>& nodiPoly,
                    vector<int>& latiPoly,
                    int& inverti,
                    bool& chiuso,
                    PolygonalMesh mesh,
                    Vector3d N,
                    vector<int>& idBordo,
                    vector<int>& idInterno);

// Funzione che a partire dalle informazioni della
// mesh relative alle Cell0D e alle Cell1D calcola
// tutte le Cell2D (poligoni) per la frattura considerata
//
// mesh = mesh con informazioni caricate sulle Cell0D e Cell1D
// idBordo = elenco degli id dei lati al bordo
// idInterno = elenco degli id dei lati interni
void CaricamentoCell2D (PolygonalMesh& mesh,
                       vector<int>& idBordo,
                       vector<int>& idInterno);

}
