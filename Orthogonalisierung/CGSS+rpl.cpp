
#include "Ortho.cpp"
#include <limits>
#include <random>
#include <ctime>
#include <cstdlib>
#include <cmath>

template<typename T>
array<vector<T>,3> cgs_rpl_step(matrix<T> Q,vector<T> x, T nu, int rpltol){
    srand(time(NULL));
T rho;
int n = Q.size1();
int nq = Q.size2();
vector<T> r (nq);
for(int i = 0; i< nq; i++)
r(i)= 0.0;
T nux = norm_2(x);
vector<T> y(n);
bool zeronorm;


if(nq==0){
    if(nux == 0){
    for(int i = 0; i < n ; i++)
    y(i) =(((rand()%100))/100.01)-0.5;
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
T nu = nux;



if(nux != 0){
    zeronorm = false;
    y = x/nux;
    nu = nu /nux;
}
else{
    zeronorm = true;
    for(int i = 0; i < n ; i++)
    y(i) =((rand()%100)/100.01);
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

    if(nu2 > rpltol* nu * std::numeric_limits<T>::epsilon())
    nu1 = nu2;
    else{
        nu = nu * std::numeric_limits<T>::epsilon();
        nu1 = nu;
       for(int i = 0; i < n ; i++)
        y(i) = ((rand()%100)/100.01);
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
array<matrix<T>,2>   cgs_rpl(matrix<T> X,int rpltol){
    int m = X.size1();
    int s = X.size2();
    matrix<T> Q(m,s);
     matrix<T> R(s,s);
    for(int i = 0; i < m; i++){
        for(int j = 0; j < s; j++){
            Q(i,j)= 0.0;
            R(i,j)=0.0;
        }
    }
   
    T tt = 0.0;
      
    vector<T> w = column(X,0);
    array<vector<T>,3> too = cgs_rpl_step(matrix<T>(m,0),w,tt,rpltol);
    column(Q,0)= too[0];
    R(0,0) = too[2](0);

    for(int k=0; k < s-1; k++){
        vector<T> w= column(X,k+1);
        array<vector<T>,3> too= cgs_rpl_step(Q,w,tt,rpltol);
        column(Q,k+1)= too[0];
        column(R,k+1)= too[1];
        R(k+1,k+1) = too[2](0);
    }
    /*for(int i = 0; i < X.size2();i++){
        vector<T> qq = norm_2(column(X,i))*column(Q,i);
        column(Q,i) = qq;
    }*/

    array<matrix<T>,2> foo {Q,R};
    
return  foo;
}



void TestCGSS(){
    matrix<cpp_dec_float_50> MM(3,3);
MM(0,0) =  cpp_dec_float_50("4354556569999999999999");
MM(0,1) =  cpp_dec_float_50("99494494985745999999999999");
MM(0,2) =  cpp_dec_float_50("444444444444444444499999999999");
MM(1,0) =  cpp_dec_float_50("121212121212121212122299999999999");
MM(1,1) =  cpp_dec_float_50("57573837463728485940398299999999999");
MM(1,2) =  cpp_dec_float_50("10000000000000000000000000999999999999");
MM(2,0) =  cpp_dec_float_50("4000404430033030404040404049999999999999");
MM(2,1) =  cpp_dec_float_50("99999999999999999999999999999999999999999999");
MM(2,2) =  cpp_dec_float_50("100000000000000000000000000000011111111111111111");

array<matrix<cpp_dec_float_50>,2> ar = cgs_rpl(MM,1);
matrix<cpp_dec_float_50> Q = ar[0];
matrix<cpp_dec_float_50> R = ar[1];
cout << R<<endl;
cout<<"***********************************************************************************************"<<endl; 
cout<<"***********************************************************************************************"<<endl;
matrix<cpp_dec_float_50> I = prod(trans(Q),Q);
cout << I <<endl;
cout<<"***********************************************************************************************"<<endl; 
cout<<"***********************************************************************************************"<<endl;
matrix<cpp_dec_float_50> X = prod(Q,R);
cout<<X<<endl;
}

void TestCGSSHARD(){
    matrix<cpp_dec_float_100> HARD(60,60);
 HARD = createMatrix100("svp60.txt",HARD);
array< matrix<cpp_dec_float_100>,2> mmm = cgs_rpl(HARD,1);
matrix<cpp_dec_float_100> Q = mmm[0];
matrix<cpp_dec_float_100> R = mmm[1];
 cout<<"***********************************************************************************************"<<endl; 
cout<<"***********************************************************************************************"<<endl;
cout<<"***********************************************************************************************"<<endl;
cout<<prod(trans(Q),Q)<<endl;
cout<<"***********************************************************************************************"<<endl; 
cout<<"***********************************************************************************************"<<endl;
cout<<"***********************************************************************************************"<<endl;
cout<<R<<endl;
cout<<"***********************************************************************************************"<<endl; 
cout<<"***********************************************************************************************"<<endl;
cout<<"***********************************************************************************************"<<endl;
cout<<prod(Q,R)<<endl;
}








