#include "CGSS+rpl.cpp"



template<typename T>
array<matrix<T>,3> bcgs_rpl_step(matrix<T> QQ, matrix<T> X){
   
   int s = X.size2(); 
   int sk = QQ.size2();
   bool reorth = false;
   matrix<T> r(s,1);
   T rho; 
   array<matrix<T>,3> aaa;
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
       range rk (k,k+1);
      matrix<T> Y11 = project(Y,ee,rr1);
      matrix<T> w11 = project(Y,ee,rk);

       aaa = cgs_rpl_stepM(Y11, w11, nu(k));
        project(Y,ee,rk) = aaa[0];
       project(R22,rr1,rk)= aaa[1];
        
        R22(k,k) = aaa[2](0,0);
       
       if(aaa[2](0) <= 0.5*nu(k))
       reorth = true; 
      
   }

   if (sk == 0 || !reorth){
    array<matrix<T>,3> foo1 = {Y,R12,R22};
    return foo1;
   }
   
   matrix<T> S12 = prod(trans(QQ),Y);
   Y = Y - prod(QQ,S12);
   matrix<T>  S22(s,s);
    T tt = 1.0;
   for(int k = 0; k < s; k++){
        range rr1(0,k);
        range rk(k,k+1);

        matrix<T> Y11 =  project(Y, ee,rr1);
        matrix<T> w11 = project(Y,ee,rk);
        aaa = cgs_rpl_stepM(Y11,w11,tt);
        matrix<T> r = aaa[1];
       
        
        if(aaa[2](0,0)< 0.5){
            tt = 0.0;
            Y11 = project(Y,ee,rr1);
            matrix<T> Y33 (QQ.size1(),QQ.size2()+k);
            range qa (0,QQ.size1());
            range qe (0,QQ.size2());
            project(Y33,qa,qe) = project(QQ,qa,qe);
            range yk (0,k);
            range qee (QQ.size2(),QQ.size2()+k);
            project(Y33,qa,qee) = project(Y,qa,yk);
           
            aaa = cgs_rpl_stepM(Y33,w11,tt);
 
            range rr3(0,sk);
            range orr(0,1);
            range ss(0,S12.size1()); 
            range rr4(sk,sk+k);

            project(S12,ss,rk) = project(S12,ss,rk) + project(aaa[1],rr3,orr);
            r = project(aaa[1],rr4,orr);
                 
        }
         project(Y,ee,rk) = aaa[0];
        project(S22,rr1,rk) = r;
        S22(k,k) = aaa[2](0,0);
   }      
    R12 = R12 + prod(S12,R22);
    R22 = prod(S22,R22);
    array<matrix<T>,3> foo = {Y,R12,R22};
    return foo;
}


template<typename T>
array<matrix<T>,2> bcgs_rpl(matrix<T> XX, int s){
 // auto start = std::chrono::high_resolution_clock::now();
  
int m = XX.size1();
int n = XX.size2();
int p = n/s;

matrix<T> QQ (m,n);
matrix<T> RR (m,n); 
matrix<T> Y11;
matrix<T> X11;
array<matrix<T>,3> aaa;
range kk (0,s);
range ee (0,m);
int sk = s;

X11 = project(XX,ee,kk);
aaa = bcgs_rpl_step(matrix<T>(m,0), X11);
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
    X11 = project(XX,ee,kk);
    aaa = bcgs_rpl_step(Y11, X11);
    project(QQ,ee,kk) = aaa[0];
    project(RR,rr1,kk) = aaa[1];
    project(RR,kk,kk) = aaa[2];
    sk = sk +s;
      
}
 
array<matrix<T>,2> foo = {QQ,RR};
   // auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> float_ms = end - start;
    //cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;  
    //cout<<"\n"<<matrixNormInfI(foo[0])<<endl;
//cout<<"\n"<<matrixNormInfQR(XX,foo[0],foo[1])<<endl;
return foo;}





















