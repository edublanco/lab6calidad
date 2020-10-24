
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>
#include <cstdlib>
using namespace std;

const long double PI = 3.14159265358979323846;



//=================================================================================================
//====================================----LAB 3----================================================
//=================================================================================================

struct Vertex{

    long double xk =0;
    long double x = 0;
    long double y =0;
    long double x2 =0;
    long double y2 = 0;
    long double xy = 0;
};

class CalculoLab3{

private:

    bool flagN = true;
//---------------------- checa si esta vacio ----------

//.m=3
bool docVacio(ifstream& doc){

    if ( doc.peek() == ifstream::traits_type::eof() )
        {
          return true;
        }
    return false;
}
//------------------------------------------------------



public:
    string S;
    long double contador =0;
    string tempxk = "";
    string tempx = "";
    string tempy = "";
    long double x =0;
    long double y=0;
    long double n = 0;
    bool fcoma = false;
    vector<Vertex> v ;
    bool numeros = false;

    long double B1(long double sumxy, long double N, long double avx, long double avy,  long double sumx2  ){

        long double parteArriba = sumxy -(N*avx*avy);
        long double parteAbajo = sumx2 -(N*(avx*avx));
        long double b1 = parteArriba/parteAbajo;

        return b1;
    }

    long double B0(long double avy,  long double b1, long double avx  ){
        long double b0 = avy -(b1*avx);

        return b0;
    }

    long double R(long double sumx,  long double sumy, long double sumxy, long double sumx2, long double sumy2, long double N  ){

       long double parteArriba = (N *sumxy)-(sumx * sumy);
       long double parteAbajoDer = (N*sumy2)-(sumy*sumy);
       long double parteAbajoIzq =  (N*sumx2)-(sumx*sumx);
       long double parteAbajo = parteAbajoDer * parteAbajoIzq;
       long double parteAbajoRaiz =  sqrt(parteAbajo);
       long double r = parteArriba/parteAbajoRaiz;


        return r;
    }

    long double R2(long double r  ){
        long double r2 = r * r;
        return r2;
    }

    long double Yk(long double b0,  long double b1, long double xk  ){
        long double yk = b0 +(b1*xk);
        return yk;
    }


    void xyk(ifstream& doc,string name){
        contador = 0 ;

        if(docVacio(doc)){
            cout<<"----------------------------------------------\n";
            cout<<"------- documento vacio o inexistente---------\n";
            cout<<"----------------------------------------------\n";
            exit(0);
        }

        while ( getline(doc, S) ) {
                   v.push_back({0, 0 ,0, 0,0,0});
            //getline(doc, S);
            v[contador].x = 0;
            v[contador].y = 0;
            tempx ="";
            tempy ="";
            n=0;

            if((S[n] == ',') || (S[n] == '.') || (S[n] == '0') || (S[n] == '1')|| (S[n] == '2')|| (S[n] == '3') || (S[n] == '4')|| (S[n] == '5')|| (S[n] == '6')|| (S[n] == '7')|| (S[n] == '8')|| (S[n] == '9')){
                numeros = true;

            }else{
                numeros = false;
            }
            if(!numeros){
                cout<<"---------------------------------------\n";
                cout<<"-----programa con input incorrecto-----\n";
                cout<<"---------------------------------------\n";
                exit(0);
            }
            else if(contador ==0){//------------------- sio es la primera linea va aguardar xk
                while( (S[n] == '0') || (S[n] == '1')|| (S[n] == '2')|| (S[n] == '3') || (S[n] == '4')|| (S[n] == '5')|| (S[n] == '6')|| (S[n] == '7')|| (S[n] == '8')|| (S[n] == '9')){
                    tempxk += S[n] ;
                    n++;
                }
            }else  {
                tempxk = "0" ;
                while((S[n] == ',') || (S[n] == '.') || (S[n] == '0') || (S[n] == '1')|| (S[n] == '2')|| (S[n] == '3') || (S[n] == '4')|| (S[n] == '5')|| (S[n] == '6')|| (S[n] == '7')|| (S[n] == '8')|| (S[n] == '9')){

                        if(S[n] == ',') { // -------checa si hay coma
                            fcoma = true;
                            n++;
                        }else if(fcoma) { //--------si hubo coma por el flag de fcoma, lo guarda en tempy y cada loop le va jhcaiendo un append del siguiente numero
                            tempy += S[n];
                            n++;
                        }
                        else{//-------------------- si no ha habido coma lo guarda en la string  tempx y cada loop le va jhcaiendo un append del siguiente numero
                            tempx += S[n] ;
                            n++;
                        }
                }

                fcoma =false;
                x =stof(tempx); // ---- convierte el string x a flotante y lo guarda en el vertice
                y= stof(tempy); // ---- convierte el string y a flotante y lo guarda en el vertice
                v[contador].x += x;
                v[contador].y += y;
                v[contador].x2 += (x*x);
                v[contador].y2 += (y*y);
                v[contador].xy += (x*y);
            }

            v[contador].xk += stof(tempxk);// --- convierte el estring tempxk flotante y lo guarda en el vertice
            contador ++;
        }
    }
};

class Respuesta{

    public:


           void calc( long double b1, long double b0, long double r, long double r2, long double yk, long double xk, long double N, long double sig,long double ran, long double ls, long double li ){

                printf("N  = %.0f \n", (float)N);
                printf("xk = %.0f \n", (float)xk);
                printf("r  = %.5f \n", (float)r);
                printf("r2 = %.5f \n",(float) r2);
                printf("b0 = %.5f \n",(float) b0);
                printf("b1 = %.5f \n",(float) b1);
                printf("yk = %.5f \n", (float)yk);
                printf("ran= %.5f \n", (float)ran);
                printf("sig= %.10f \n", (float)sig);
                printf("LS = %.5f \n", (float)ls);
                printf("LI = %.5f \n", (float)li);

        }
};


class Lab3 {

private:


    vector<Vertex> r;
    CalculoLab3 calculo;
    Respuesta respuesta;
    long double sumx =0;
    long double sumy =0;
    long double sumx2 =0;
    long double sumy2 =0;
    long double sumxy =0;

    long double avy = 0;
    long double N =0;
    long double b1 =0;
    long double b0=0;
    long double r1=0 ;
    long double r2=0 ;
    long double yk =0;
    long double xk =0;

public:
    long double avx =0;

    long double * funLab3(string archivo){
        ifstream doc(archivo);
        int temp = 0;

        static long double r[6];

        calculo.xyk(doc,archivo);
        N = calculo.contador -1;
        xk= calculo.v[0].xk;

        for(int i =0; i < calculo.contador; i++ ){

            sumx += calculo.v[i].x;
            sumx2 += calculo.v[i].x2;
            sumy += calculo.v[i].y;
            sumy2 += calculo.v[i].y2;
            sumxy += calculo.v[i].xy;

        }
        avx = sumx/N;
        avy = sumy/N;
        b1 = calculo.B1(sumxy,N,avx,avy,sumx2);
        b0 = calculo.B0(avy,b1,avx);
        r1 =calculo.R(sumx,sumy,sumxy,sumx2,sumy2,N);
        r2 = calculo.R2(r1);
        yk = calculo.Yk(b0, b1, xk);

        r[0] = b1;
        r[1] = b0;
        r[2] = r1;
        r[3] = r2;
        r[4] = yk;
        r[5] = xk;
        r[6] = N;

        return r;
    }
};


//=================================================================================================
//====================================----END LAB 3----============================================
//=================================================================================================

//=================================================================================================
//====================================----LAB 5----================================================
//=================================================================================================
class CalculoLab5{
private:
//--------------------------------------------------------------------
     long double pAbajo(long double dof){

       long double aux1 = pow((dof*PI),0.5);
       long double aux3 = dof/2.0;
       long double aux2 = tgamma(aux3);

       return (pow((dof*PI),0.5))* aux2;
     }

//--------------------------------------------------------------------
     long double pArriba(long double dof){
        return tgamma((dof+1)/2);
     }
//--------------------------------------------------------------------
    long double pDerecha(long double x,long double dof){
        long double aux1 =(1+((pow(x,2))/dof));
        return pow((1+((pow(x,2))/dof)),-((dof+1)/2));

     }
//--------------------------------------------------------------------

public:
//--------------------------------------------------------------------
    long double simpson(long double x, long double dof, long double num_seg){
        long double p1 =0;
        long double p2 =0;

        long double w = x / num_seg;
        //cout<<"w: "<<w<<"\n";
        long double i =1;
        while(i<= (num_seg - 1 )){
            long double f = ((pArriba(dof)/pAbajo(dof))*pDerecha((w*i),dof));
            p1 += 4.0*f*(w/3);
            i += 2.0;
        }
        long double j =2;

        while(j<= (num_seg - 2 )){
            long double f = ((pArriba(dof)/pAbajo(dof))*pDerecha((w*j),dof));
            p2 +=2.0*f*(w/3);
            j += 2;

        }
        return (((w/3)*((pArriba(dof)/pAbajo(dof))*pDerecha(0, dof))) + p1 + p2 + ((w/3)*((pArriba(dof)/pAbajo(dof))*pDerecha(x, dof))));

    }

};


//=================================================================================================

class P{

    private:

    long double  E = 0.0000001;
    long double num_seg = 10;
    long double p =1;
    long double r1 = 0;
    long double r2 = 0;
    long double oldR = 0;

    CalculoLab5 calculo;


public:


    long double calcP(long double x,long double dof){  // regresa p T(x, DOF)
        if( !(x > 0 || dof > 0) ){
            cout <<"-----datos incorrecto ---------"<< "\n";
            exit(0);
        }

        while(p > E){
            if(r1== 0){
                r1 = calculo.simpson(x,dof, num_seg);
            }

            num_seg = num_seg*2;
            r2 = calculo.simpson(x,dof, num_seg);
            p = r1 - r2;
            r1 = p;
            oldR = r2;

        }
        return oldR;


        }

};

class Lab5{

public:
    long double funLab5(long double pt ,long double dof ){ // regresa x T(P, DOF)

        long double d0 =0; // borre 3 lnieas y agregue 3 y modifique 1
        long double x = 1;

        long double p = 0;
        long double  E = 0.0000001;
        int contador =0;

        P p0;



        d0= x/2;

        if(pt <= 0 || pt > .5){
            cout <<"-----datos incorrecto ---------"<< "\n";
            exit(0);
        }

        p = p0.calcP(x,dof);

        while( abs(p - pt) > E  ){
            P p1;
            if(p < pt){
                x = x + d0;
                p = p1.calcP(x,dof);
                contador++;
            }else if(p > pt ){
                d0= d0 /2;
                x = x - d0;
                p = p1.calcP(x,dof);;
                contador++;
            }


        }
        return x;


    }

};


//=================================================================================================
//====================================----END LAB 5----============================================
//=================================================================================================


class Rango{
private:

    CalculoLab3 calculo;
    Lab5 lab5;

    long double alfa(string archivo, long double n, long double b1, long double b0){
        long double auxDer = 0 ;
        long double auxIzq = 1 / (n-2) ;
        ifstream doc(archivo);
        calculo.xyk(doc,archivo);


        for(int i = 1; i <= n; i ++){
            auxDer +=   (calculo.v[i].y - b0  - (b1*calculo.v[i].x)) * (calculo.v[i].y - b0  - (b1*calculo.v[i].x));
        }

        return sqrt(auxDer * auxIzq);

    }

    long double parteArriba(long double xk, long double xavg){

        return (xk - xavg) * (xk - xavg);
    }

    long double parteAbajo(long double xavg, string archivo, long double n){

        ifstream doc(archivo);
        calculo.xyk(doc,archivo);
        long double auxAbajo = 0;

        for(int i = 1; i <= n; i ++){

            auxAbajo +=   (calculo.v[i].x - xavg) * (calculo.v[i].x - xavg);
        }
        return auxAbajo;
    }



public:
    long double rang(long double n, long double xk , long double xavg,long double b1, long double b0, string archivo ){

        long double der = sqrt(1 + (1/n) + (parteArriba(xk, xavg)/parteAbajo(xavg,archivo,n)));
        long double izq = lab5.funLab5(.35,n - 2 );

        return alfa(archivo,n, b1,b0) * der * izq;
    }

    long double ls(long double r, long double yk){
        return r + yk;
    }

    long double li(long double r, long double yk){
        long double aux;

        aux =yk - r;

        if (aux<0){
            return 0;
        }else{
            return aux;
        }
    }

};




int main(){
    string archivo;
    cin >> archivo ;

    Lab3 lab3;
    P pt;
    Lab5 lab5;
    Rango rango;
    CalculoLab3 calculo;
    Respuesta respuesta;

    long double *p;
    p = lab3.funLab3(archivo);

    long double  N = *(p + 6);
    long double b1 = *(p + 0);
    long double b0= *(p + 1);
    long double r1= *(p + 2) ;
    long double r2= *(p + 3);
    long double yk = *(p + 4);
    long double xk = *(p + 5);
    long double xavg = lab3.avx;
    long double r = rango.rang(N,xk,xavg,b1,b0,archivo);
    long double ls = rango.ls(r, yk);
    long double li = rango.li(r, yk);
    long double sig= 0;
    //long double yavg = lab3.avy;

    //--------------------------------sig-----------------
    long double izq =0;
    long double totalx =0 ;
    long double abajo = 0;

    izq = r1 * (sqrt(N-2));
    abajo = sqrt(1-r2);
    totalx = izq / abajo;
    sig = 1.0 - (2.0* pt.calcP(totalx, N-2));

    //----------------------------- end sig----------------------

    respuesta.calc(b1,b0, r1, r2,yk,xk,N,sig,r, ls,li );

    //---------------------------------------------------------

}



