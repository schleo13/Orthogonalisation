#include <iostream>                             // Standardstream-Funktionaliät einbinden
#include <fstream>
#include <string>
#include <sstream>                              // ofstream und ifstream einbinden

#include <boost/numeric/ublas/matrix_vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
using std::cout;
using std::endl;
using namespace boost::numeric::ublas;
using  boost::multiprecision::cpp_dec_float_50;
using  boost::multiprecision::cpp_dec_float_100;
using std::cerr;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ios_base;
using std::array;

template<typename T>
matrix<T> createMatrixDouble1(string txt, matrix<T> M){
    ifstream Quelldatei;
    Quelldatei.open(txt, ios_base::in);
    if(!Quelldatei)
    cerr<< "Konnte nicht geöffnet werden\n" << endl;
    char c;
    int i=0;
    int j=0;  
    string s ="";
    stringstream ss;      
    while(Quelldatei.get(c)){       
        if(c ==' '){
            double n;
             ss >> n;
            M(j,i)= n;
           cout<<" "<<n;
            i++;
             ss.clear();   
        }
        if(c == '\n'){
            double n;
             ss >> n;
            M(j,i)= n;
           cout<<" "<<n;          
             ss.clear();
            j++;
            i=0;
        }
        if(c > 47 && c<58  ){
            //c = c+48;
            ss<<c;
            //cout<<" "<<(int)c-48;
             
        }       
    }
    return M;
}
template<typename T>
matrix<T> createMatrixDouble(string txt, matrix<T> M){
    ifstream Quelldatei;
    Quelldatei.open(txt, ios_base::in);
    if(!Quelldatei)
    cerr<< "Konnte nicht geöffnet werden\n" << endl;
    char c;
    int i=0;
     int j=0;  
    string s ="";
    stringstream ss;      
    while(Quelldatei.get(c)){
        if(c ==' ' ){
            double n;    
             ss >> n;
            M(j,i)= n;
           //cout<<" "<<s;
            i++;
             ss.clear();  
        }  
        if(c == ']' ){
            if(j == M.size1())break;
            double n;
             ss >> n;
            M(j,i)= n;
           //cout<<" "<<s;
            i++;
             ss.clear();  
            j++;
            i=0;
        }
        if(c > 47 && c<58  ){
            //c = c+48;
            ss<<c;
            //cout<<" "<<(int)c-48;
             
        }       
    }
    return M;
}

template<typename T>
matrix<T> createMatrix100(string txt, matrix<T> M){
    ifstream Quelldatei;
    Quelldatei.open(txt, ios_base::in);
    if(!Quelldatei)
    cerr<< "Konnte nicht geöffnet werden\n" << endl;
    char c;
    int i=0;
     int j=0;  
    string s ="";
    stringstream ss;      
    while(Quelldatei.get(c)){
        if(c ==' ' ){    
             ss >> s;
            M(j,i)= cpp_dec_float_100(s);
           //cout<<" "<<s;
            i++;
             ss.clear();  
        }  
        if(c == ']' ){
            if(j == M.size1())break;
             ss >> s;
            M(j,i)= cpp_dec_float_100(s);
           //cout<<" "<<s;
            i++;
             ss.clear();  
            j++;
            i=0;
        }
        if(c > 47 && c<58  ){
            //c = c+48;
            ss<<c;
            //cout<<" "<<(int)c-48;
             
        }       
    }
    return M;
}

template<typename T>
matrix<T> createMatrix50(string txt, matrix<T> M){
    ifstream Quelldatei;
    Quelldatei.open(txt, ios_base::in);
    if(!Quelldatei)
    cerr<< "Konnte nicht geöffnet werden\n" << endl;
    char c;
    int i=0;
    int j=0;  
    string s ="";
    stringstream ss;      
    while(Quelldatei.get(c)){       
        if(c ==' '){
             ss >> s;
            M(j,i)=cpp_dec_float_50(s);
           cout<<" "<<s;
            i++;
             ss.clear();   
        }
        if(c == '\n'){
             ss >> s;
            M(j,i)= cpp_dec_float_50(s);
           cout<<" "<<s;          
             ss.clear();
            j++;
            i=0;
        }
        if(c > 47 && c<58  ){
            //c = c+48;
            ss<<c;
            //cout<<" "<<(int)c-48;
             
        }       
    }
    return M;
}


