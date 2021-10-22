#include <iostream>                            
#include <fstream>
#include <string>
#include <sstream>                             
#include <array>
#define NDEBUG
#include <boost/numeric/ublas/matrix_vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/gmp.hpp>
using std::cout;
using std::endl;
using namespace boost::numeric::ublas;
using  boost::multiprecision::cpp_dec_float_50;
using  boost::multiprecision::cpp_dec_float_100;
using  boost::multiprecision::mpz_int;
using std::cerr;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ios_base;
using std::array;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<150>> cpp_dec_float_150;
typedef boost::multiprecision::number<boost::multiprecision::gmp_float<0>> mpff;

matrix<mpz_int> createMatrix(string txt, matrix<mpz_int> M){
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
        if(c ==' ' || c==',' ){    
             ss >> s;
            M(j,i)= mpz_int(s);
           //cout<<" "<<s;
            i++;
             ss.clear();  
        }  
        if(c=='\n' ){
            if(j == M.size1())break;
             ss >> s;
            M(j,i)= mpz_int(s);
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

matrix<cpp_dec_float_150> createMatrix(string txt, matrix<cpp_dec_float_150> M){
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
        if(c ==' ' || c==',' ){    
             ss >> s;
            M(j,i)= cpp_dec_float_150(s);
           //cout<<" "<<s;
            i++;
             ss.clear();  
        }  
        if(c=='\n' ){
            if(j == M.size1())break;
             ss >> s;
            M(j,i)= cpp_dec_float_150(s);
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

matrix<mpff> createMatrix(string txt, matrix<mpff> M){
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
        if(c ==' ' || c==',' ){    
             ss >> s;
            M(j,i)= mpff(s);
           //cout<<" "<<s;
            i++;
             ss.clear();  
        }  
        if(c=='\n' ){
            if(j == M.size1())break;
             ss >> s;
            M(j,i)= mpff(s);
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




matrix<double> createMatrix(string txt, matrix<double> M){
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
         if(c ==' ' || c==',' ){    
            double n;    
             ss >> n;
            M(j,i)= n;
           //cout<<" "<<s;
            i++;
             ss.clear();  
        }  
        if(c =='\n' ){
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


matrix<double> createMatrixDoubleOctave(string txt, matrix<double> M){
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
         if( c==',' ){    
            double n;    
             ss >> n;
            M(j,i)= n;
           //cout<<" "<<s;
            i++;
             ss.clear();  
        }  
        if(c =='\n' ){
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




matrix<cpp_dec_float_50> createMatrix50(string txt, matrix<cpp_dec_float_50> M){
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
         if(c ==' ' || c==',' ){    
             ss >> s;
            M(j,i)= cpp_dec_float_50(s);
           //cout<<" "<<s;
            i++;
             ss.clear();  
        }  
        if(c == ']' || c=='\n' ){
            if(j == M.size1())break;
             ss >> s;
            M(j,i)= cpp_dec_float_50(s);
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







