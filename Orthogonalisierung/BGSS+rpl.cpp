#include "CGSS+rpl.cpp"


using  boost::math::round;
template<typename T>
array<matrix<T>,3> bcgs_rpl_step(matrix<T> QQ, matrix<T> X){
   int s = X.size2(); 
   int sk = QQ.size2();
   bool reorth = false;
   vector<T> y(s);
   vector<T> r(s);
   T rho; 
   array<vector<T>,3> aaa;
   vector<T> nu (s);
   for(int k = 0; k <  s; k++)
       nu(k)= norm_2(column(X,k));

   matrix<T> R12 = prod(trans(QQ),X);  
   matrix<T> Y = X - prod(QQ,R12);
   matrix<T> R22(s,s);

    int n = Y.size1(); 
    range ee(0,n);

   for(int k = 0; k < s ; k++){       
       range rr1(0,k);
       range rr2(0,s+1);
      matrix<T> Y11 = project(Y,ee,rr1);
      vector<T> w11 = column(Y,k);
       aaa = cgs_rpl_step(Y11, w11, nu(k));

        column(Y,k) = aaa[0];
        vector<T> rr = column(R22,k);
        rr = project(rr,rr1);
        rr = aaa[1];
        column(R22,k) = r;
        R22(k,k) = aaa[2](0);
       
       if(aaa[2](0) <= 0.5*nu(k))
       reorth = true; 
   }

   if (sk == 0 || !reorth){
    array<matrix<T>,3> foo1 = {Y,R12,R22};
    return foo1;
   }
   
   matrix<T> S12 = prod(trans(QQ),Y);
   
   Y = Y - prod(QQ,S12);
   T tt = 1.0;
   matrix<T>  S22(s,s);

   for(int k = 0; k < s; k++){
     
       range rr1(0,k);
       vector<T> w11(k);
       matrix<T> Y11(n,k);
        Y11 =  project(Y, ee,rr1);
        w11 = column(Y,k);
        aaa = cgs_rpl_step(Y11,w11,tt); 

        if(aaa[2](0) < 0.5){
            range rr3 (0, sk);
            Y11 = project(Y,ee,rr1);
            matrix<T> Y33 (QQ.size1(),QQ.size2()+k);
            for(int i = 0; i < QQ.size2();i++)
            column(Y33,i) = column(QQ,i);
            for(int i = 0; i < k;i++)
            column(Y33,i+QQ.size2()) = column(Y,i);
                tt = 0.0;
            aaa = cgs_rpl_step(Y33,w11,tt);
            column(S12,k) = column(S12,k) + project(aaa[1],rr3);
            range rr4 (sk+1, sk+k);
            r = project(aaa[1],rr4);
        }
        column(Y,k) = aaa[0];
         vector<T> rr = column(S22,k);
        rr = project(rr,rr1);
        rr = r;
        column(S22,k) = r;
        S22(k,k) = aaa[2](0);
   }
         
        

    R12 = R12 + prod(S12,R22);
    R22 = prod(S22,R22);
    array<matrix<T>,3> foo = {Y,R12,R22};
    return foo;

}

template<typename T>
array<matrix<T>,2> bcgs_rpl(matrix<T> XX, int s){
int m = XX.size1();
int n = XX.size2();
int p = n/s;

matrix<T> QQ (m,n);
matrix<T> RR (m,n);

array<matrix<T>,3> aaa;
range kk (0,s);
range ee (0,m);
int sk = s;

matrix<T> Y11 = project(XX,ee,kk);
aaa = bcgs_rpl_step(matrix<T>(m,0), Y11);
project(QQ,ee,kk) = aaa[0];
project(RR,kk,kk) = aaa[2];
int as = 0;
int bs = s;

for(int k = 0; k <p-1; k++){
    as = as +s;
    bs = bs +s;
     range kk (as, bs);
    range rr1(0,sk);
    Y11 = project(QQ,ee,rr1);
    matrix<T> X11 = project(XX,ee,kk);
    aaa = bcgs_rpl_step(Y11, X11);
    project(QQ,ee,kk) = aaa[0];
    
    project(RR,rr1,kk) = aaa[1];
        
    project(RR,kk,kk) = aaa[2];
    
    sk = sk +s;
    
}
array<matrix<T>,2> foo = {QQ,RR};
cout<<"block_rpl_____Q = "<<QQ<<"\n"<<endl;
    cout<<"block_rpl_____R = "<<RR<<"\n"<<endl;
return foo;}





















