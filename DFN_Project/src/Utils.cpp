#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

namespace DFNlibrary {
const double tau = 1e-6;

bool ImportaFratture(const string& filepath,
                     vector<Fracture>& fract,
                     int& numF){

    ifstream file;
    file.open(filepath);

    if(file.fail()){
        cerr << "il file di input non è stato aperto correttamente" << endl;
        return false;
    }

    vector<string> inputLines;
    string line;
    while (getline(file, line))
        inputLines.push_back(line);

    file.close();


    string s = inputLines[1];

    numF = stoi(s); // numero di fratture

    if (numF == 0)
    {
        cerr << "non ci sono fratture" << endl;
        return false;
    }

    fract.reserve(numF);

    int id;
    int numV;

    for (int i=0; i<numF; i++){

        Fracture F;
        int pos = i*6+3; // idicatore della posizione della frattura corrente i
                         // ogni frattura occupa 6 righe nel file di input,
                         // inoltre sono presenti 3 righe con informazioni generali all'inizio
        istringstream converter(inputLines[pos]);
        getline(converter, s, ';');
        id = stoi(s); // FractureId
        F.id = id;

        getline(converter, s);
        numV = stoi(s); // NumVertices
        F.NumVertices = numV;

        F.coordx.reserve(numV);
        F.coordy.reserve(numV);
        F.coordz.reserve(numV);

        pos+=2;
        istringstream converterX(inputLines[pos]); // coordinate asse x dei vertici della frattura

        for(int k =0; k<numV; k++){
            if(k == numV-1){
                getline(converterX, s);
                F.coordx.push_back(stod(s));
            }
            else{
                getline(converterX, s, ';');
                F.coordx.push_back(stod(s));
            }
        }

        pos++;
        istringstream converterY(inputLines[pos]); // coordinate asse y dei vertici della frattura

        for(int k =0; k<numV; k++){
            if(k == numV-1){
                getline(converterY, s);
                F.coordy.push_back(stod(s));
            }
            else{
                getline(converterY, s, ';');
                F.coordy.push_back(stod(s));
            }
        }

        pos++;
        istringstream converterZ(inputLines[pos]); // coordinate asse z dei vertici della frattura

        for(int k =0; k<numV; k++){
            if(k == numV-1){
                getline(converterZ, s);
                F.coordz.push_back(stod(s));
            }
            else{
                getline(converterZ, s, ';');
                F.coordz.push_back(stod(s));
            }
        }

        fract.push_back(F);

    }

    file.close();
    return true;
}

Vector3d CalcoloNormale(const Fracture F){

    // calcolo del vettore ortonormale N al piano contenente F

    Vector3d v; // vettore (F vertice 1 - F vertice 0)
    v[0] = F.coordx[1] - F.coordx[0]; // coordinata x di v
    v[1] = F.coordy[1] - F.coordy[0]; // coordinata y di v
    v[2] = F.coordz[1] - F.coordz[0]; // coordinata z di v

    Vector3d u; // vettore (F vertice 2 - F vertice 0)
    u[0] = F.coordx[2] - F.coordx[0]; // coordinata x di u
    u[1] = F.coordy[2] - F.coordy[0]; // coordinata y di u
    u[2] = F.coordz[2] - F.coordz[0]; // coordinata z di u

    Vector3d N; // v x u (prodotto scalare)
    N[0] = v[1]*u[2] - v[2]*u[1];
    N[1] = v[2]*u[0] - v[0]*u[2];
    N[2] = v[0]*u[1] - u[0]*v[1];

    // normalizzazione
    N = N/ N.norm();

    return N;
}




bool CercaIntersezioni(Vector3d P, Vector3d t, Fracture F, double& c1, double& c2){
    bool interseca = false; // "interseca" indica se è già stata trovata una delle due intersezioni
    for(int i = 0; i< F.NumVertices; i++){
        // k = vettore vertice i+1 - vertice i
        Vector3d k = {F.coordx[(i+1)%F.NumVertices]-F.coordx[i], F.coordy[(i+1)%F.NumVertices]-F.coordy[i], F.coordz[(i+1)%F.NumVertices]-F.coordz[i]};
        // il modulo serve ad associare all'ultima iterata i vertici con indice pari a F.NumVertices-1 e 0 (F.NumVertices % F.NumVertices)
        // ovvero per un quadrilatero il vertice 3 e il vertice 0

        // valCond = prodotto scalare t*k - prodotto delle norme di t e k
        // valCond = 0 se e solo se i due vettori t e k sono paralleli
        // (valCond è una differenza relativa)
        double valCond = (abs(t[0]*k[0] + t[1]*k[1] + t[2]*k[2]) - t.norm()*k.norm())/(t.norm()*k.norm());
        if(abs(valCond) < tau) // k e t sono paralleli, non puo' esserci intersezione
        {
            continue; // se i due vettori sono paralleli si passa all'iterata successiva
        }

        Vector3d b = {F.coordx[i]-P[0], F.coordy[i]-P[1], F.coordz[i]-P[2]};
        MatrixXd A(3, 2);
        A.col(0) = t;
        A.col(1) = k;
        Vector2d sol = A.fullPivLu().solve(b); // sistema lineare per l'intersezione
        // solo se -sol(1) appartiene a (0,1) vi è intersezione tra la retta ed il segmento
        if((0< -sol(1)) && (-sol(1)<1)){
            if (interseca){
                c2 = sol(0);
                break; // dopo aver trovato la seconda intersezione il ciclo si conclude
                       // (si suppone che tutte le fratture siano convesse)
            }
            else{
                interseca = true; // dopo aver trovato la prima intersezione "interseca" diventa vera
                c1 = sol(0);
            }
        }
    }

    return interseca;
}



void InserisciTraccia(double alpha, double beta, double gamma, double delta,
                      vector<Traces>& tracesContainer, Vector3d P, Vector3d t, int Fid1, int Fid2){
    double max1 = max(alpha,beta);
    double min1 = min(alpha,beta);
    double max2 = max(gamma,delta);
    double min2 = min(gamma,delta);

    // vengono esaminati 6 possibili scenari di intersezione tra i segmenti
    // suddivisi in casi "A" e "B" quando sono speculari
    // non occore relativizzare le grandezze mini e maxi con i=1,2
    // in quanto tali valori appartengono a (0,1) e sono già normalizzati

    // caso 0 A e B
    // le due fratture non si intersecano
    if( (min1-max2)> -tau || (min2-max1)> -tau ){
        return;
    }
    Vector3d estremo1 = {};
    Vector3d estremo2 = {};
    // Tips1 e Tips2:
    bool T1;
    bool T2;

    // caso 1 A
    // il segmento 1 contiene il segmento 2
    if( (max1-max2)>tau && (min2-min1)>tau){
        estremo1 = P + max2*t;
        estremo2 = P + min2*t;
        T1 = true;
        T2 = false;
    }

    // caso 1 B
    // il segmento 2 contiene il segmento 1 nella parte interna
    else if((max2-max1)>tau && (min1-min2)>tau){
        estremo1 = P + max1*t;
        estremo2 = P + min1*t;
        T1 = false;
        T2 = true;
    }

    // caso 2
    // i due segmenti coincidono
    else if(abs(max1-max2)<tau && abs(min1-min2)<tau){
        estremo1 = P + max1*t;
        estremo2 = P + min1*t;
        T1 = false;
        T2 = false;
    }

    // caso 3
    // i due segmenti hanno l'estremo superiore in comune, ma non quello inferiore
    else if(abs(max1-max2)<tau){
        //caso 3 A
        if(min2>min1){
            estremo1 = P + max1*t;
            estremo2 = P + min2*t;
            T1 = true;
            T2 = false;
        }
        //caso 3 B
        else{
            estremo1 = P + max1*t;
            estremo2 = P + min1*t;
            T1 = false;
            T2 = true;
        }
    }

    // caso 4
    // i due segmenti hanno l'estremo inferiore in comune, ma non quello superiore
    else if(abs(min1-min2)<tau){
        //caso 4 A
        if(max1>max2){
            estremo1 = P + max2*t;
            estremo2 = P + min1*t;
            T1 = true;
            T2 = false;
        }
        //caso 4 B
        else{
            estremo1 = P + max1*t;
            estremo2 = P + min1*t;
            T1 = false;
            T2 = true;
        }
    }

    // caso 5 A
    // i due segmenti hanno solo una porzione in comune
    else if( (max1-max2)>tau ){
        estremo1 = P + max2*t;
        estremo2 = P + min1*t;
        T1 = true;
        T2 = true;
    }

    // caso 5 B
    // i due segmenti hanno solo una porzione in comune
    else{
        estremo1 = P + max1*t;
        estremo2 = P + min2*t;
        T1 = true;
        T2 = true;
    }

    Traces T = {};
    T.id = tracesContainer.size();
    T.FractureID1 = Fid1;
    T.FractureID2 = Fid2;
    T.P1 = estremo1;
    T.P2 = estremo2;
    T.Tips1 = T1;
    T.Tips2 = T2;
    T.Length = (estremo1-estremo2).norm();
    // la traccia viene memorizzata ed aggiunta al contenitore "tracesContainer"
    tracesContainer.push_back(T);
}


void CercaTracce(const Fracture F1, const Fracture F2, vector<Traces>& tracesContainer){

    // calcolo dei vettori ortonormali ai piani F1, F2
    Vector3d N1 = CalcoloNormale(F1);
    Vector3d N2 = CalcoloNormale(F2);

    // calcolo direzione tangente t della retta intersezione dei due piani,
    // prodotto vettoriale tra N1 e N2
    Vector3d t;
    t[0] = N1[1]*N2[2] - N1[2]*N2[1];
    t[1] = N1[2]*N2[0] - N1[0]*N2[2];
    t[2] = N1[0]*N2[1] - N2[0]*N1[1];

    if (t.norm() < tau) // condizione: t == 0, in algebra finita con tolleranza tau
                        // (non serve relativizzare, dato che N1 e N1 sono già normali)
    {
        return; // i due piani sono paralleli, non possono esserci tracce
    }

    // prodotto scalare tra N1 e F1 vertice 0
    double d1 = N1[0]*F1.coordx[0] + N1[1]*F1.coordy[0] + N1[2]*F1.coordz[0];

    // prodotto scalare tra N2 e F2 vertice 0
    double d2 = N2[0]*F2.coordx[0] + N2[1]*F2.coordy[0] + N2[2]*F2.coordz[0];

    // creazione matrice A = {N1, N2, t}
    Matrix3d A = {};
    A.row(0) = N1;
    A.row(1) = N2;
    A.row(2) = t;

    // creazione vettore termine noto b = {d1, d2, 0}'
    Vector3d b = {d1, d2, 0};

    // risoluzione del sistema tramite la decomposizione PALU,
    // P punto appartente alla retta di intersezione dei due piani
    Vector3d P = A.fullPivLu().solve(b);

    double alpha = 0;
    double beta = 0;

    // se anche solo una delle due fratture non ha intersezioni con la retta, la funzione termina

    if(!CercaIntersezioni(P,t,F1,alpha,beta))
        return;


    double gamma = 0;
    double delta = 0;

    if(!CercaIntersezioni(P,t,F2,gamma,delta))
        return;

    // se entrambe le fratture intersecano la retta intersezione dei piani
    // potrebbe esserci una traccia
    InserisciTraccia(alpha, beta, gamma, delta, tracesContainer, P, t, F1.id, F2.id);

}


void StampaTracce(vector<Traces> tracesContainer, int numF){
    // fileName = "traces_FRX_data.csv", X = numero fratture nel file considerato
    string fileName = "traces_FR"+to_string(numF)+"_data.csv";
    ofstream outputTraces (fileName);
    outputTraces << "# Number of Traces"<<endl;    //intestazione
    outputTraces << tracesContainer.size()<<endl;
    outputTraces << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;

    for(auto& t: tracesContainer){
        outputTraces<<t.id<<";"<<t.FractureID1<<";"<<t.FractureID2<<";";
        outputTraces<<t.P1[0]<<";"<<t.P1[1]<<";"<<t.P1[2]<<";"<<t.P2[0]<<";"<<t.P2[1]<<";"<<t.P2[2]<<endl;
    }

    outputTraces.close();
}



void StampaTracceOrdinate(vector<Traces> tracesContainer, int numFracture,
                          map<int, vector<int>>& sortedPassanti,
                          map<int, vector<int>>& sortedNonPassanti){

    // fileName = "sorted_traces_FRX_data.csv", X = numero fratture nel file considerato
    string fileName = "sorted_traces_FR"+to_string(numFracture)+"_data.csv";
    ofstream outputSortedTraces (fileName);

    for(int i=0; i<numFracture; i++){
        // contenitori di tracce in forma breve, Vector2d: [TraceId, Length]
        vector<ShortTraces> tPassanti = {};
        vector<ShortTraces> tNonPassanti = {};

        // divisione delle tracce tra passanti e non-passanti relativamente alle fratture
        for(auto& t: tracesContainer){

            if(t.FractureID1 == i){
                ShortTraces shortT(t.id, t.Length);
                if(t.Tips1 == true){
                    tNonPassanti.push_back(shortT);
                }
                else{
                    tPassanti.push_back(shortT);
                }
            }

            if(t.FractureID2 == i){
                ShortTraces shortT(t.id, t.Length);
                if(t.Tips2 == true){
                    tNonPassanti.push_back(shortT);
                }
                else{
                    tPassanti.push_back(shortT);
                }
            }
        }

        // ordinamento decrescente (sulla lunghezza) delle tracce passanti
        if(tPassanti.size() != 0){
            QuickSort(tPassanti);
        }

        // ordinamento decrescente (sulla lunghezza) delle tracce non-passanti
        if(tNonPassanti.size() != 0){
            QuickSort(tNonPassanti);
        }

        vector<int> idTP = {}; // id tracce passanti
        vector<int> idTNP = {}; // id tracce non passanti

        // stampa delle tracce ordinate sul file di output più caricamento
        // dei contenitori "sortedPassanti" e "sortedNonPassanti"
        if(tPassanti.size() + tNonPassanti.size()!= 0){
            outputSortedTraces << "# FractureId; NumTraces"<<endl;
            outputSortedTraces << i<<";"<< tPassanti.size() + tNonPassanti.size()<<endl;
            outputSortedTraces << "# TraceId; Tips; Length"<<endl;

            if(!tPassanti.empty()){
                for(auto& t: tPassanti){
                    outputSortedTraces<<t.id<<";0;"<<t.Length<<endl;
                    idTP.push_back(t.id);
                }
                sortedPassanti.insert({i, idTP});
            }

            if(!tNonPassanti.empty()){
                for(auto& t: tNonPassanti){
                    outputSortedTraces<<t.id<<";1;"<<t.Length<<endl;
                    idTNP.push_back(t.id);
                }
                sortedNonPassanti.insert({i, idTNP});
            }
        }
    }
    outputSortedTraces.close();
}



void IntersezioneEdges(edges L1, edges L2, double& alpha, double& beta){

    alpha = 42; // valore sentinella
    beta = 42;  // valore sentinella

    // valCond = prodotto scalare t*k - prodotto delle norme di t e k
    // valCond = 0 se e solo se i due vettori t e k sono paralleli
    // (valCond è una differenza relativa)
    double valCond = (abs(L1.t[0]*L2.t[0] + L1.t[1]*L2.t[1] + L1.t[2]*L2.t[2]) - L1.t.norm()*L2.t.norm()) / (L1.t.norm()*L2.t.norm());
    if(abs(valCond)<tau) // k e t sono paralleli, non puo' esserci intersezione
    {
        return; // se i lati sono paralleli non può esserci intersezione e quindi la funzione termina
    }
    Vector3d b = L2.P - L1.P;
    MatrixXd A(3, 2);
    A.col(0) = L1.t;
    A.col(1) = L2.t;
    Vector2d sol = A.fullPivLu().solve(b); // sistema lineare per l'intersezione
    // alpha e beta sono i coefficienti dell'intersezione:
    // I (intersezione) = L1.P + alpha*L1.t = L2.P + beta.L2.t
    alpha = sol[0];
    beta = -sol[1];
}


void TagliaFratture(Fracture F, vector<Traces> contenitoreTracce,
                    map<int, vector<int>> tPO, map<int, vector<int>> tNPO,
                    vector<edges>& latiBordo, vector<edges>& latiInterni){
    // tPO = tracce Passanti Ordinate, tNPO = tracce Non-Passanti Ordinate

    // caricamento latiBordo
    for(int i=0; i<F.NumVertices; i++){
        Vector3d P = {F.coordx[i], F.coordy[i], F.coordz[i]};
        // il modulo serve a trovare l'ultimo lato, che collega l'ultimo nodo col primo
        Vector3d t = {F.coordx[(i+1)%F.NumVertices]-F.coordx[i], F.coordy[(i+1)%F.NumVertices]-F.coordy[i], F.coordz[(i+1)%F.NumVertices]-F.coordz[i]};
        edges e = {};
        e.P = P;
        e.t = t;
        latiBordo.push_back(e);
    }

    if(tPO.find(F.id) != tPO.end()){
        // caricamento delle tracce passanti e delle
        // relative intersezioni con i lati del bordo e
        // tra le stesse tracce passanti
        for(int i=0; i<tPO[F.id].size(); i++){
            edges l = {};
            int idTraccia = tPO[F.id][i];
            // estrazione dati della traccia
            for(auto& T: contenitoreTracce){
                if(T.id == idTraccia){
                    l.P = T.P1;
                    l.t = T.P2 - T.P1;
                    break;
                }
            }

            bool flagStart = false;
            bool flagEnd = false;

            // ricerca dei lati del bordo intersecati dalle tracce passanti
            // e delle relative intersezioni
            for(auto& b: latiBordo){
                // alpha e beta sono già valori relativizzati secondo
                // la lunghezza di b.t, quindi si possono ammettere
                // disugualianze assolute
                double alpha;
                double beta;
                if(flagStart && flagEnd){
                    break; // una volta trovare entrambe le intersezioni il ciclo finisce
                }

                IntersezioneEdges(b,l,alpha,beta);
                if(abs(beta) < tau){ // condizione beta = 0
                    if(0<alpha && alpha<1){
                        b.intersection.push_back(alpha);
                        flagStart = true;
                    }
                }
                else if(abs(beta-1) < tau){ // condizione beta = 1
                    if(tau<alpha && alpha<1-tau){
                        b.intersection.push_back(alpha);
                        flagEnd = true;
                    }                    
                }
            }

            // caricamento delle intersezioni con i lati interni già caricati
            for(auto& i: latiInterni){
                double alpha;
                double beta;
                IntersezioneEdges(i,l,alpha,beta);
                if(tau<alpha && alpha<1-tau && tau<beta && beta<1-tau){
                    i.intersection.push_back(alpha);
                    l.intersection.push_back(beta);
                }
            }

            latiInterni.push_back(l); // caricamento del lato interno
        }
    }

    // Caricamento delle tracce non-passanti e delle
    // relative intersezioni con i lati del bordo e
    // con i lati interni già caricati.
    // Le tracce non-passanti vanno estese prima di
    // essere caricate.
    if(tNPO.find(F.id) != tNPO.end()){
        // best[0] = indice lato che viene intersecato, best[1] = intersezione l, best[2] = intersezione b oppure i
        Vector3d best = {-1,0,0};
        double j = -1; // contatore progressivo per determinatre best[0]
        for(int i=0; i<tNPO[F.id].size(); i++){
            edges l = {};
            int idTraccia = tNPO[F.id][i];
            for(auto& T: contenitoreTracce){
                if(T.id == idTraccia){
                    l.P = T.P1;
                    l.t = T.P2 - T.P1;
                    break;
                }
            }
            bool doppiaEstensione = true;
            // doppiaEstensione = true se la traccia ha entrambi gli estremi non giacenti su un lato
            // doppiaEstensione = false se la traccia ha uno dei due estremi giacenti su un lato

            // valuto orientazione
            for(auto& b: latiBordo){
                double alpha;
                double beta;
                IntersezioneEdges(b,l,alpha,beta);
                if (abs(beta)<tau){ // caso già orientato
                    if(tau<alpha && alpha<1-tau){
                        b.intersection.push_back(alpha);
                    }
                    doppiaEstensione = false;
                    break;
                }
                else if(abs(beta-1)<tau){ // caso con orientazione da invertire
                    if(tau<alpha && alpha<1-tau){
                        b.intersection.push_back(alpha);
                    }
                    // inversione della traccia, ovvero l'estremo considerato viene
                    // posto come l'altro rimanente, e la tangente viene invertita di senso
                    l.P += l.t;
                    l.t = -l.t;
                    doppiaEstensione = false;
                    break;
                }
            }

            if(doppiaEstensione == true){ // se non è ancora stato trovato il lato su cui giace uno dei due estremi
                for(auto& i: latiInterni){
                    double alpha;
                    double beta;
                    IntersezioneEdges(i,l,alpha,beta);
                    if (abs(beta)<tau){ // caso già orientato
                        if(tau<alpha && alpha<1-tau){
                            i.intersection.push_back(alpha);
                        }
                        doppiaEstensione = false;
                        break;
                    }
                    else if(abs(beta-1)<tau){ // caso con orientazione da invertire
                        if(tau<alpha && alpha<1-tau){
                            i.intersection.push_back(alpha);
                        }
                        doppiaEstensione = false;
                        l.P += l.t;
                        l.t = -l.t;
                        break;
                    }
                }
            }

            // valuto intersezioni prima dell'estensione
            // nota: tali intersezioni possono avvenire solo con lati interni
            // dato che l'intersezione col lato al bordo è già stata caricata
            for(auto& i: latiInterni){
                double alpha;
                double beta;
                IntersezioneEdges(i,l,alpha,beta);
                if(tau<alpha && alpha<1-tau && tau<beta && beta<1-tau){ // alpha, beta appartenenti ad (0,1)
                    i.intersection.push_back(alpha);
                    l.intersection.push_back(beta);
                }
            }

            // ricerca dell'estremo mancante, ovvero completamneto
            // della traccia non-passante, fino a farla diventare passante
            best = {-1,0,0};
            j = -1;
            for(auto& b: latiBordo){
                j++;
                double alpha;
                double beta;
                IntersezioneEdges(b,l,alpha,beta);
                if(1+tau<beta && tau<alpha && alpha<1-tau){// alpha appartenente ad (0,1), 1<beta
                    if(beta<best[1] || best[1]==0){
                        best = {j,beta,alpha};
                    }
                }

            }

            for(auto& i: latiInterni){
                j++;
                double alpha;
                double beta;
                IntersezioneEdges(i,l,alpha,beta);
                if(1+tau<beta && tau<alpha && alpha<1-tau){// alpha appartenente ad (0,1), 1<beta
                    if(beta<best[1] || best[1]==0){
                        best = {j,beta,alpha};
                    }
                }
            }

            // aggiornamento di l
            l.t *= best[1];
            for(auto& inter: l.intersection){
                inter = inter/best[1];
            }


            // caso doppia enstensione
            if(doppiaEstensione){
                // caricamento intersezione nel contenitore latiBordo o latiInterni
                if(best[0] >= latiBordo.size()){
                    best[0] -= latiBordo.size();
                    latiInterni[best[0]].intersection.push_back(best[2]);
                }
                else{
                    latiBordo[best[0]].intersection.push_back(best[2]);
                }

                // inversione della traccia
                l.P += l.t;
                l.t = -l.t;

                // ricerca dell'estremo mancante (a questo punto è diventata una traccia non passante singola)
                best = {-1,0,0};
                j = -1;
                for(auto& b: latiBordo){
                    j++;
                    double alpha;
                    double beta;
                    IntersezioneEdges(b,l,alpha,beta);
                    if(1+tau<beta && tau<alpha && alpha<1-tau){// alpha appartenente ad (0,1), 1<beta
                        if(beta<best[1] || best[1]==0){
                            best = {j,beta,alpha};
                        }
                    }
                }

                for(auto& i: latiInterni){
                    j++;
                    double alpha;
                    double beta;
                    IntersezioneEdges(i,l,alpha,beta);
                    if(1+tau<beta && tau<alpha && alpha<1-tau){// alpha appartenente ad (0,1), 1<beta
                        if(beta<best[1] || best[1]==0){
                            best = {j,beta,alpha};
                        }
                    }
                }

                // aggiornamento di l
                l.t *= best[1];
                // le intersezioni vanno prima ribaltate rispetto agli estremi
                for(auto& inter: l.intersection){
                    inter = (1-inter)/best[1];
                }
            }

            latiInterni.push_back(l);

            // caricamento dell'intersezione nel contenitore latiBordo o latiInterni
            if(best[0] >= latiBordo.size()){
                best[0] -= latiBordo.size();
                latiInterni[best[0]].intersection.push_back(best[2]);
            }
            else{
                latiBordo[best[0]].intersection.push_back(best[2]);
            }
        }
    }
}


void CercaEstremo(Vector3d P, int& id, vector<Vector3d> CoordinateNodi,
                  vector<int> idNodi, bool& flag){
    flag = false;
    int indice = -1;
    for(auto& n: CoordinateNodi){
        indice++;
        if(((n-P).norm()/P.norm())<tau){ // condizione n = P
            flag = true;
            id = idNodi[indice]; // gli id dei nodi possono non essere progressivi
            break;
        }
    }
}

// funzione senza la flag tra gli argomenti
void CercaEstremo(Vector3d P, int& id, vector<Vector3d> CoordinateNodi, vector<int> idNodi){
    bool flag;
    CercaEstremo(P,id,CoordinateNodi,idNodi,flag);
}




void CaricamentoCell0e1D(vector<edges>& latiBordo, vector<edges>& latiInterni, PolygonalMesh& mesh,
                                                  vector<int>& idBordo, vector<int>& idInterno){
    // contenitori temporanei
    vector<int> tempCell0DId = {};
    vector<Vector3d> tempCell0DCoordinates = {};
    vector<vector<int>> tempInvolvedEdges = {};
    vector<int> tempCell1DId = {};
    vector<Vector2i> tempCell1DVertices = {};

    // caricamento tempCell0D e tempCell0DCoordinates per i lati del bordo
    int n = -1;
    for(auto& b: latiBordo){
        n++;
        tempCell0DId.push_back(n);
        tempCell0DCoordinates.push_back(b.P);
        if(!b.intersection.empty()){
            // ordinamento decrescente
            QuickSort(b.intersection);
            for(int j = b.intersection.size()-1; j>-1; j--){
                n++;
                tempCell0DId.push_back(n);
                tempCell0DCoordinates.push_back(b.P+b.intersection[j]*b.t);
            }
        }
    }

    // caricamento tempInvolvedEdges e tempCell1DVertices per i lati del bordo
    tempInvolvedEdges.push_back({n,0});  // primo vertice
    tempCell1DVertices.push_back({0,1}); // primo lato
    for(int j = 1; j<n; j++){
        tempInvolvedEdges.push_back({j-1,j});
        tempCell1DVertices.push_back({j,j+1});
    }
    tempInvolvedEdges.push_back({n-1,n}); // ultimo vertice
    tempCell1DVertices.push_back({n,0});  // ultimo lato

    // caricamento identificativi lati bordo
    // (l'id è progressivo)
    for(int j =0; j<n+1; j++){
        tempCell1DId.push_back(j);
    }


    int l = n;
    for(auto& i: latiInterni){
        int idEstremo1 = -1;
        int idEstremo2 = -1;
        CercaEstremo(i.P, idEstremo1, tempCell0DCoordinates, tempCell0DId);
        if(!i.intersection.empty()){
            // ordinamento decrescente
            QuickSort(i.intersection);
            for(int j = i.intersection.size()-1; j>-1; j--){
                l++;
                bool trovato;
                CercaEstremo(i.P+i.intersection[j]*i.t, idEstremo2, tempCell0DCoordinates, tempCell0DId, trovato);
                if(!trovato){
                    n++;
                    tempCell0DId.push_back(n);
                    tempCell0DCoordinates.push_back(i.P+i.intersection[j]*i.t);
                    idEstremo2 = n;
                }
                tempCell1DVertices.push_back({idEstremo1,idEstremo2});

                if(n>=tempInvolvedEdges.size()){
                    tempInvolvedEdges.resize(n*2);
                }
                tempInvolvedEdges[idEstremo1].push_back(l);
                tempInvolvedEdges[idEstremo2].push_back(l);
                idEstremo1 = idEstremo2;
            }
        }
        CercaEstremo(i.P+i.t, idEstremo2, tempCell0DCoordinates, tempCell0DId);
        l++;
        tempCell1DVertices.push_back({idEstremo1,idEstremo2});
        tempInvolvedEdges[idEstremo1].push_back(l);
        tempInvolvedEdges[idEstremo2].push_back(l);
    }

    // caricamento identificativi lati interni
    // (l'id è progressivo)
    idBordo = tempCell1DId;
    for(int j = idBordo.size(); j<l+1; j++){
        tempCell1DId.push_back(j);
        idInterno.push_back(j);
    }

    // caricamento della mesh
    mesh.NumberCell0D = tempCell0DId.size();
    mesh.Cell0DId = tempCell0DId;
    mesh.Cell0DCoordinates = tempCell0DCoordinates;
    tempInvolvedEdges.resize(mesh.NumberCell0D);
    mesh.InvolvedEdges = tempInvolvedEdges;
    mesh.NumberCell1D = tempCell1DId.size();
    mesh.Cell1DId = tempCell1DId;
    mesh.Cell1DVertices = tempCell1DVertices;
}



Vector3d CalcoloNormaleMesh(PolygonalMesh mesh){
    Vector3d origine = mesh.Cell0DCoordinates[0];
    Vector3d v1 = mesh.Cell0DCoordinates[1]-origine;
    Vector3d v2 = {};
    for(int i = 2; i<mesh.NumberCell0D; i++){
        v2 = mesh.Cell0DCoordinates[i]-origine;
        // valCond = prodotto scalare v1*v2 - prodotto delle norme di v1 e v2
        // valCond = 0 se e solo se i due vettori v2 e v2 sono paralleli
        // (valCond è una differenza relativa)
        double valCond = abs(abs(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]) - v1.norm()*v2.norm())/(v1.norm()*v2.norm());
        if(valCond > tau){
            // i due vettori non sono paralleli
            break;
        }
    }
    // prodotto vettoriale v1 x v2 = N
    double Nx = v1[1]*v2[2] - v1[2]*v2[1];
    double Ny = v1[2]*v2[0] - v1[0]*v2[2];
    double Nz = v1[0]*v2[1] - v1[1]*v2[0];
    Vector3d N = {Nx, Ny, Nz};
    N = N/N.norm();
    return N;
}


// !!!ATTENZIONE!!!
// da questo momento in avanti ci si riferirà, con
// un abuso di notazione, agli archi con il
// termine "lati"


void LatoSuccessivo(Vector3d& CurrentEdgeTan, int& CurrentNode, int& CurrentEdgeId,
                    vector<int>& nodiPoly, vector<int>& latiPoly, int& inverti, bool& chiuso,
                    PolygonalMesh mesh, Vector3d N, vector<int>& idBordo, vector<int>& idInterno){
    // inverti è un intero che può assumere tre valori:
    // 0: condizione di default, nessuna operazione associata
    // 1: il lato è da invertire
    // 2: il lato è già stato invertito, non può succedere
    // che si debba invertire un'altra volta per lo stesso poligono
    if(inverti != 2){
        inverti = 1;
    }
    double bestValue = -2;
    int idBestLato = -1;
    int idBestNodoEsterno = -1;
    Vector3d bestluTan = {}; // miglior tangente lato uscente
    int idNodoEsterno = -1;

    // prodotto vettoriale N x CurrentEdgeTan = versoAntiorario
    double vax = N[1]*CurrentEdgeTan[2] - N[2]*CurrentEdgeTan[1];
    double vay = N[2]*CurrentEdgeTan[0] - N[0]*CurrentEdgeTan[2];
    double vaz = N[0]*CurrentEdgeTan[1] - N[1]*CurrentEdgeTan[0];
    Vector3d versoAntiorario = {vax, vay, vaz};
    versoAntiorario = versoAntiorario/versoAntiorario.norm();
    for(auto& lu: mesh.InvolvedEdges[CurrentNode]){ //lu = id lati uscenti
        if(lu == CurrentEdgeId){
            continue; // il lato successivo non può essere il lato corrente
        }
        // controllo che il lato lu sia presente tra i lati al bordo o interni
        bool trovatoBordo = (find(idBordo.begin(),idBordo.end(), lu) != idBordo.end());
        bool trovatoInterno = (find(idInterno.begin(),idInterno.end(), lu) != idInterno.end());
        if(!(trovatoBordo || trovatoInterno)){
            // lu non trovato
            continue;
        }

        // assegnazione della posizione rispetto al lato
        // del secondo nodo (ovvero quello non corrente)
        // del lato, denominato NodoEsterno
        if(CurrentNode == mesh.Cell1DVertices[lu][0]){
            idNodoEsterno = mesh.Cell1DVertices[lu][1];
        }
        else{
            idNodoEsterno = mesh.Cell1DVertices[lu][0];
        }
        Vector3d luTan = mesh.Cell0DCoordinates[idNodoEsterno] - mesh.Cell0DCoordinates[CurrentNode];
        Vector3d normluTan = luTan/luTan.norm();
        // condizione: prodotto calare (normluTan * versoAntiorario) >= 0
        // i due vettori "normluTan" e "versoAntiorario" sono già normalizzati
        if(normluTan[0]*versoAntiorario[0]+normluTan[1]*versoAntiorario[1]+normluTan[2]*versoAntiorario[2]> -tau){
            if(inverti != 2){
                inverti = 0;
            }
            // newValue = (versoAntiorario x normluTan) * N
            // (x = prodotto vettoriale, * = prodotto scalare)
            double x = versoAntiorario[1]*normluTan[2] - versoAntiorario[2]*normluTan[1];
            double y = versoAntiorario[2]*normluTan[0] - versoAntiorario[0]*normluTan[2];
            double z = versoAntiorario[0]*normluTan[1] - versoAntiorario[1]*normluTan[0];
            double newValue = x*N[0] + y*N[1] + z*N[2];
            if(newValue > bestValue){
                bestValue = newValue;
                idBestLato = lu;
                idBestNodoEsterno = idNodoEsterno;
                bestluTan = luTan;
            }
        }
    }

    // controllo se il nuovo nodo appartiene già al
    // poligono: se così fosse il poligono risulta
    // chiuso e quindi completo
    for(auto& idN: nodiPoly){
        if(idN == idBestNodoEsterno){
            chiuso = true;
        }
    }
    // altrimenti si indica che il poligono è ancora incompleto
    if(!chiuso){
        nodiPoly.push_back(idBestNodoEsterno);
    }

    // aggiornamento dei valori di nodo e lato corrento
    CurrentEdgeTan = bestluTan;
    CurrentNode = idBestNodoEsterno;
    latiPoly.push_back(idBestLato);
    CurrentEdgeId = idBestLato;
}


void CaricamentoCell2D (PolygonalMesh& mesh, vector<int>& idBordo, vector<int>& idInterno){
    int id2D = -1; // id cell 2D
    Vector3d N = CalcoloNormaleMesh(mesh);
    while(!idBordo.empty()){ // il processo si conclude quando non ci sono più lati nel contenitore idBordo
        id2D++;
        // inizializzazione poligono in costruzione
        int idEdgeIniz = idBordo[0];
        int CurrentEdgeId = idEdgeIniz;
        vector<int> latiPoly = {idEdgeIniz};
        vector<int> nodiPoly = {mesh.Cell1DVertices[idEdgeIniz][0], mesh.Cell1DVertices[idEdgeIniz][1]};
        Vector3d CurrentEdgeTan = mesh.Cell0DCoordinates[nodiPoly[1]] - mesh.Cell0DCoordinates[nodiPoly[0]];
        int CurrentNode = nodiPoly[1];

        // inverti è un intero che può assumere tre valori:
        // 0: condizione di default, nessuna operazione associata
        // 1: il lato è da invertire
        // 2: il lato è già stato invertito, non può succedere
        // che si debba invertire un'altra volta per lo stesso poligono
        int inverti = 0;
        bool chiuso = false;

        while(!chiuso){ // il processo si arresta quando il poligono è completo (chiuso)
            if(inverti == 1){
                // bisogna invertire il senso in cui è orientato il primo lato,
                // il processo riparte da capo
                latiPoly = {idEdgeIniz};
                nodiPoly = {mesh.Cell1DVertices[idEdgeIniz][1], mesh.Cell1DVertices[idEdgeIniz][0]}; //inversione
                CurrentEdgeTan = mesh.Cell0DCoordinates[nodiPoly[1]] - mesh.Cell0DCoordinates[nodiPoly[0]];
                CurrentNode = nodiPoly[1];
                inverti = 2; // una volta invertito il processo non si invertirà una seconda volta
                CurrentEdgeId = idEdgeIniz;
            }
            LatoSuccessivo(CurrentEdgeTan, CurrentNode, CurrentEdgeId, nodiPoly,
                           latiPoly, inverti, chiuso, mesh, N,
                           idBordo,idInterno);
        }

        // caricamento del nuovo poligono nella mesh
        mesh.Cell2DId.insert({id2D, latiPoly.size()});
        mesh.Cell2DVertices.push_back(nodiPoly);
        mesh.Cell2DEdges.push_back(latiPoly);

        // i lati del poligono appena formato vengono rimossi
        // dal contenitore idInterno e aggiunti al contenitore
        // idBordo (se in principio si trovavano in idInterno),
        // oppure vengo rimossi dal contenitore idBordo
        for(auto& l: latiPoly){
            vector<int>::iterator it = find(idInterno.begin(),idInterno.end(), l);
            if(it != idInterno.end()){
                idInterno.erase(it);
                idBordo.push_back(l);
            }
            else{
                idBordo.erase(find(idBordo.begin(), idBordo.end(), l));
            }
        }
    }
    mesh.NumberCell2D = mesh.Cell2DId.size();
}

}
