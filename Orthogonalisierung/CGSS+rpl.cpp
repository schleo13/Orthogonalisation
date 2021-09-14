
#include "Ortho.cpp"
#include <limits>
#include <random>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <chrono>

template<typename T>
array<vector<T>,3> cgs_rpl_step(matrix<T> Q,vector<T> x,T nu){
     srand(time(NULL));
     int rpltol = 1;
T rho;
int n = Q.size1();
int nq = Q.size2();
vector<T> r (nq);
T nux = norm_2(x);
vector<T> y(n);
bool zeronorm;

if(nq==0){
    if(nux == 0){
    for(int i = 0; i < n ; i++)
    y(i) =(((rand()%100))/109.0)-0.5;
    y = y/norm_2(y);   
    rho = 0.0;
    }
    else{
      y = x/nux;
      rho = nux;  
    }
vector<T> rrho(1);
rrho(0)=rho;
array<vector<T>,3> foo {y,r,rrho};
return foo;
}
 
 if(nu < nux)
 nu = nux;
if(nux != 0){
    zeronorm = false;
    y = x/nux;
    nu = nu /nux;
}
else{
    zeronorm = true;
    for(int i = 0; i < n ; i++)
    y(i) =((rand()%100)/109.0)-0.5;
    y = y/norm_2(y);
    nu = 1.0;
}

T nu1 = nu;
T nu2;

while(true){
  
    vector<T> s = prod(trans(Q),y);
    r = r + s;
    y = y - prod(Q,s);
    nu2 = norm_2(y);

    if(nu2 > 0.5*nu1)
    break;

    if(nu2 > rpltol* nu * std::numeric_limits<T>::epsilon()){
        nu1 = nu2;
    } 
    else{
        nu = nu * std::numeric_limits<T>::epsilon();
        nu1 = nu;
       for(int i = 0; i < n ; i++)
        y(i) = ((rand()%100)/109.0)-0.5;
        y = nu * (y/norm_2(y));
        
    }
    // if(boost::math::isnan::isnan(y))
    // rho = FP_NAN;
}
if(!zeronorm){
    rho = norm_2(y);
    y = y/rho;
    rho = rho * nux;
    r = r * nux;
}
else{
    y = y / norm_2(y);
    for(int i = 0; i < nq; i++)
    r(i)= 0.0;
    rho = 0.0;
}
vector<T> rrho(1);
rrho(0)=rho;
array<vector<T>,3> foo {y,r,rrho};
return foo;
}



template<typename T>
array<matrix<T>,2>   cgs_rpl(matrix<T> X){
    //auto start = std::chrono::high_resolution_clock::now();
    int rpltol = 1;
    int m = X.size1();
    int s = X.size2();
    matrix<T> Q(m,s);
     matrix<T> R(s,s);
T tt = 0.0;
    vector<T> w = column(X,0);
    array<vector<T>,3> too = cgs_rpl_step(matrix<T>(m,0),w,tt);
    column(Q,0)= too[0];
    R(0,0) = too[2](0);
    
    for(int k=0; k < s-1; k++){
        range ee(0,Q.size1());
        range k1(0,k);
        matrix<T> Q1 = project(Q,ee,k1);
        vector<T> w= column(X,k+1);
        array<vector<T>,3> too= cgs_rpl_step(Q,w, tt);
        column(Q,k+1)= too[0];
        
        column(R,k+1)= too[1];
        R(k+1,k+1) = too[2](0);
    }
    /*for(int i = 0; i < X.size2();i++){
        vector<T> qq = norm_2(column(X,i))*column(Q,i);
        column(Q,i) = qq;
    }*/
        // auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> float_ms = end - start;
//cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;  
    //cout<<"cgs_rpl_____Q = "<<Q<<"\n"<<endl;
    //cout<<"cgs_rpl_____R = "<<R<<"\n"<<endl;

    array<matrix<T>,2> foo {Q,R};
    
return  foo;
}

void ZeitNormMessen(){ //ntru *2
    int p = 11;
    string ss = "";
    string s1 = "svpc-e1-d";
    string s2 = ".txt";
        vector<double> Tcgs (p); //150=56;
        vector<double> NcgsI (p);
        vector<double> Ncgs (p);
        vector<double> Tmgs (p); //150=56;
        vector<double> NmgsI (p);
        vector<double> Nmgs (p);
        vector<double> Tcgs_rpl (p); //150=56;
        vector<double> Ncgs_rplI (p);
        vector<double> Ncgs_rpl (p);
        vector<double> THouse (p); //150=56;
        vector<double> NHouseI (p);
        vector<double> NHouse (p);
        
        array<matrix<double>,2> foo; //
     int k = 0;
     for(int i = 6; i <= 26; i = i+2 ){   
        matrix<double> X(i,i);
         stringstream sss;
        sss << i;
        string str = sss.str();
        ss = s1 + str + s2;
          X = createMatrixDouble(ss,X);
   
         auto start = std::chrono::high_resolution_clock::now();
        foo = CGS(X);
       
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> float_ms = end - start;
        cout<<"i = "<<i<<"__"<<"Benötigte Zeit : "<<float_ms.count()<<endl;
        Tcgs(k)= float_ms.count();
        NcgsI(k)=  matrixNormInfI(foo[0]);
         Ncgs(k)=  matrixNormInfQR(X,foo[0],foo[1]);

           start = std::chrono::high_resolution_clock::now();
        foo = MGS(X);
         end = std::chrono::high_resolution_clock::now();
         float_ms = end - start;
        cout<<"i = "<<i<<"__"<<"Benötigte Zeit : "<<float_ms.count()<<endl;
        Tmgs(k)= float_ms.count();
        NmgsI(k)=  matrixNormInfI(foo[0]);
         Nmgs(k)= matrixNormInfQR(X,foo[0],foo[1]);

        start = std::chrono::high_resolution_clock::now();
        foo = cgs_rpl(X);
         end = std::chrono::high_resolution_clock::now();
         float_ms = end - start;
        cout<<"i = "<<i<<"__"<<"Benötigte Zeit : "<<float_ms.count()<<endl;
        Tcgs_rpl(k)= float_ms.count();
        Ncgs_rplI(k)= matrixNormInfI(foo[0]);
        Ncgs_rpl(k)=  matrixNormInfQR(X,foo[0],foo[1]);
       

        start = std::chrono::high_resolution_clock::now();
        foo = Householder(X);
         end = std::chrono::high_resolution_clock::now();
         float_ms = end - start;
        cout<<"i = "<<i<<"__"<<"Benötigte Zeit : "<<float_ms.count()<<"\n"<<endl;
        THouse(k)= float_ms.count();
        NHouseI(k)= matrixNormInfI(foo[0]);
         NHouse(k)=  matrixNormInfQR(X,foo[0],foo[1]);

        k++;
     }
     cout<<"svp-e1-d*****TIME***** CGS(1);MGS(2);cgs_rpl(3);Householder(4)"<<endl;
     cout<<""<<Tcgs<<endl;
    cout<<"\n"<<Tmgs<<endl;
     cout<<"\n"<<Tcgs_rpl<<endl;
     //cout<<"\n"<<Tbrpl<<endl;
      cout<<"\n"<<THouse<<"\n"<<endl;
    cout<<"*****INF_NORM_I***** ||I-Q^TQ|| *****"<<endl;
    cout<<""<<NcgsI<<endl;
    cout<<"\n"<<NmgsI<<endl;
     cout<<"\n"<<Ncgs_rplI<<endl;
     // cout<<"\n"<<Nbrpl<<endl;
      cout<<"\n"<<NHouseI<<"\n"<<endl;
        cout<<"*****INF_NORM_QR***** ||X-QR|| *****"<<endl;
    cout<<""<<Ncgs<<endl;
    cout<<"\n"<<Nmgs<<endl;
     cout<<"\n"<<Ncgs_rpl<<endl;
     // cout<<"\n"<<Nbrpl<<endl;
      cout<<"\n"<<NHouse<<"\n"<<endl;
}



int main(){
ZeitNormMessen(); 




    return 0;
}






















