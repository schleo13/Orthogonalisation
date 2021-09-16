#include "readMatrix.hpp"




//CGS algorithm returns QR
template <typename T>
array<matrix<T>,2> CGS(matrix<T> X){
    //auto start = std::chrono::high_resolution_clock::now();
    matrix<T> Q(X.size1(),X.size2());
     matrix<T> R(X.size1(),X.size1());
    for(int k = 0; k < X.size1();k++){
        vector<T> x = column(X,k);
        for(int j = 0; j < k; j++){
            R(j,k) = inner_prod(column(Q,j),x);
        }
        for(int j = 0; j < k; j++){
            x = x - R(j,k)*column(Q,j);
        }
        R(k,k)=norm_2(x);
        column(Q,k)= x/R(k,k);
    }
     
    array<matrix<T>,2> foo{Q,R};
    
    //auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> float_ms = end - start;
//cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;    
//cout<<"CGS_________Q = "<<Q<<"\n"<<endl;
//cout<<"CGS_________R = "<<R<<"\n"<<endl;
    return foo;
}

//MGS algorithm returns QR
template<typename T>
array<matrix<T>,2> MGS(matrix<T> X){
   //  auto start = std::chrono::high_resolution_clock::now();
    matrix<T> Q (X.size1(),X.size2());
    matrix<T> R(X.size1(),X.size1());
    for(int k = 0; k < X.size2(); k++ ){
        vector<T> x = column(X,k);
        for(int j = 0; j < k; j++){
            R(j,k) = inner_prod(column(Q,j),x);
            x = x - R(j,k)*column(Q,j);
        }
        R(k,k)=norm_2(x);
        column(Q,k)= x/R(k,k);
    }
    
    array<matrix<T>,2> foo = {Q,R};
    //    auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> float_ms = end - start;
    //cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;  
    //cout<<"MGS_____Q = "<<Q<<"\n"<<endl;
    //cout<<"MGS_____R = "<<R<<"\n"<<endl;
    return foo;
}

//Householder algorithm returns QR
template <typename T>
array<matrix<T>,2> Householder(matrix<T> X){
    // auto start = std::chrono::high_resolution_clock::now();
    vector<T> x(X.size1());
    matrix<T> Q(X.size1(), X.size2());
    matrix<T> R = X;
     for(int vc = 0; vc < Q.size1(); vc++)
     Q(vc,vc) = 1.0;
    for(int j =0; j < X.size2(); j++){
        range je (j,X.size2());
        T nu = norm_2(project(column(R,j),je));
        T s = R(j,j);
        if(s < 0){s =-1;}
        else{ s = 1;}
        T u1 = R(j,j)-s*nu;
        vector<T> w = project(column(R,j),je)/u1;
        w(0)=1;
        T tau = -s*u1/nu;
        range ee(0, X.size2());
        
        matrix<T> R1 = project(R,je,ee);
        vector<T> Rr = prod(w,R1);
         matrix<T> Q1 = project(Q,ee,je);
        vector<T> Qr = prod(Q1,w);
         w = tau*w;
        project(R,je,ee)=project(R,je,ee)-outer_prod(w,Rr);
        project(Q,ee,je)=project(Q,ee,je)-outer_prod(Qr,w);
    }
   array<matrix<T>,2> foo = {Q,R};
       // auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> float_ms = end - start;
    //cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;  
    //cout<<"Householder_Q = "<<Q<<"\n"<<endl;
    //cout<<"Householder_R = "<<R<<"\n"<<endl;
    return foo;
}


template<typename T>
matrix<T> gg(vector<T> v,int i,int j){
    matrix<T> R(v.size(),v.size());
    for(int vc = 0; vc <v.size(); vc++){
        R(vc,vc) = 1.0;
    }  
        T nn = sqrt(v(i)*v(i)+v(j)*v(j));       
        T c = v(i)/nn;        
        T s = v(j)/nn;
            R(i,i)=c;
            R(j,j)=c;
            R(i,j)=s;
            R(j,i)=s*(-1);
        
    
    return R;
}
//Givensrotation returns QR
template <typename T>
array<matrix<T>,2> Givens(matrix<T> X){
  //   auto start = std::chrono::high_resolution_clock::now();
    matrix<T> G;
     matrix<T> Q(X.size1(),X.size2());
     matrix<T> R = X;

    for(int vc = 0;vc<X.size1();vc++){
        Q(vc,vc) =1.0;
    }
    for(int i = 0; i <X.size1()-1; i++){
        vector<T> v (X.size1());
      
        for(int j = i+1; j < X.size2();j++){
            v = column(R,i);
            G = gg(v,i,j);
            Q = prod(Q,G);
            R = prod(G,R);
            //cout<<nG<<endl;
        }
    }
     array<matrix<T>,2> foo = {Q,R};
    //auto end = std::chrono::high_resolution_clock::now();
   // std::chrono::duration<double, std::milli> float_ms = end - start;
   // cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;  
    //cout<<"Givens___Q = "<<Q<<"\n"<<endl;
    //cout<<"Givens___R = "<<R<<"\n"<<endl;
    return foo;
}


template<typename T>
T matrixNormInfI(matrix<T> X){
    matrix<T>mone (X.size1(),X.size2());
    for(int vc = 0; vc < X.size2(); vc++)
    mone(vc,vc) = 1.0;
    matrix<T> mm = prod(trans(X),X);
    mm = mone-mm;

   vector<T> v(mm.size1());
    for(int i = 0; i < mm.size1(); i++){
        T c = 0.0;
       for(int j = 0; j< mm.size2();j++){
           if(mm(i,j) < 0) mm(i,j)= mm(i,j)*(-1);
           c = c + mm(i,j);
       }
       v(i)=c;
    }
    T mx = 0;
    for(int i = 0; i < v.size(); i++){
        if(v(i) > mx)
        mx = v(i);
    }
return mx;
}

template<typename T>
T matrixNormInfQR(matrix<T> X, matrix<T> Q, matrix<T>R){
    matrix<T> mm = prod(Q,R);
    mm = X-mm;
  vector<T> v(mm.size1());
    for(int i = 0; i < mm.size1(); i++){
        T c = 0.0;
       for(int j = 0; j< mm.size2();j++){
           if(mm(i,j) < 0) mm(i,j)= mm(i,j)*(-1);
           c = c + mm(i,j);
       }
       v(i)=c;
    }
    T mx = 0;
    for(int i = 0; i < v.size(); i++){
        if(v(i) > mx)
        mx = v(i);
    }
return mx;
}



template<typename T>
T matrixNormMaxI(matrix<T> X){
    matrix<T>mone (X.size1(),X.size2());
    for(int vc = 0; vc < X.size2(); vc++)
    mone(vc,vc) = 1.0;
    matrix<T> mm = prod(trans(X),X);
    mm = mone-mm;
   
    T mx = 0;
    for(int i = 0; i < X.size1(); i++){
        for(int j = 0; j < X.size2(); j++){
            T x = mm(i,j);
        if(x <0){
            x = x*(-1);
        }
        if(x > mx)
        mx = x;
        } 
    }
return mx;
}

template<typename T>
T matrixNormMaxQR(matrix<T> X, matrix<T> Q, matrix<T> R){
    matrix<T> mm = prod(Q,R);
    mm = X-mm;
    T mx = 0;
    for(int i = 0; i < X.size1(); i++){
        for(int j = 0; j < X.size2(); j++){
            T x = mm(i,j);
        if(x <0){
            x = x*(-1);
        }
        if(x > mx)
        mx = x;
        } 
    }
return mx;
}
template<typename T>
matrix<T> randMatrixI(int i, T n ){
    srand(time(NULL));
    matrix<T> X(i,i);
    vector<T> y(i);
    for(int j = 0; j < i; j++)y(j)=rand()*((i*n)*(i*n));
    column(X,0)= y;
    for(int j = 1; j < i; j++)X(j,j)=1.0;
    return X;
}
template<typename T>
matrix<T> randMatrixFull(int i, T n ){
    srand(time(NULL));
    matrix<T> X(i,i);
    
    for(int j = 0; j < i; j++){
        for(int k = 0; k <i ; k++){
            X(j,k)= rand()*((i*n)*(i*n));
        }
    }
    return X;
}



matrix<cpp_dec_float_100> Hilbert(int n, cpp_dec_float_100 a){
    matrix<cpp_dec_float_100> m(n,n);
    for(int i = 0; i < m.size1(); i++){
        for(int j = 0; j <m.size2(); j++){            
            cpp_dec_float_100 h = 1.0/(i+j+1.0);
            m(j,i) = h;
        }
    }
    return m;
}

matrix<double> Hilbert( int n, double a){
    matrix<double> m(n,n);
    for(int i = 0; i < m.size1(); i++){
        for(int j = 0; j <m.size2(); j++){            
            double h = 1.0/(i+j+1.0);
            m(j,i) = h;
        }
    }
    return m;
}

template<typename T> 
matrix<T> clearEps(matrix<T> M){
     cpp_dec_float_100 epp ("0.0000000000001");
for(int i = 0; i < M.size1();i++){
        for(int j = 0; j < M.size2();j++){
            if(M(i,j) > 0 && M(i,j) < epp || M(i,j) < 0 && M(i,j)*(-1) < epp){
                M(i,j) = 0;
            }
        }
    }
    return M;
}







