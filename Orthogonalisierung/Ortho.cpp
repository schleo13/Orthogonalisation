#include "readMatrix.hpp"


//classic CGS algorithm returns Q
template <typename T>
matrix<T> CGS(matrix<T> X){
    matrix<T> Q(X.size1(),X.size2());
    vector<T> q(X.size1());
        q = column(X,0);
        q = q/norm_2(q);
        column(Q,0)= q;
    for(int j = 1; j < X.size1();j++){
        vector<T> x(X.size1());
        x = column(X,j);       
        vector<T> px (X.size1());
        for(int i = 0; i <= j-1; i++ ){
            q = column(Q,i);    
            px =px + (inner_prod(x,q)*q);
        }
        q = x-px;
        q = q/norm_2(q);
        column(Q,j) = q;
    }
    
    return Q;
}

//classic MGS algorithm returns Q
template<typename T>
matrix<T> MGS(matrix<T> X){
    matrix<T> Q = X;
    for(int j = 0; j < Q.size1(); j++){
        vector<T> qn(Q.size1());
        qn = column(Q,j);
        qn = qn/norm_2(qn);
        column(Q,j) = qn;
        vector<T> x(Q.size1());
        for(int i = j+1; i < Q.size1(); i++){           
            x = column(Q,i);
            x = x-((inner_prod(qn,x)/inner_prod(qn,qn))*qn); //inner_prod(qn,qn)
            column(Q,i) = x;           
        }        
    }
    
    return Q;
}
//classic Householder algorithm returns Q
template <typename T>
matrix<T> Householder(matrix<T> X){
    vector<T> x(X.size1());
    matrix<T> Q(X.size1(), X.size2());
     for(int vc = 0; vc < Q.size1(); vc++)
     Q(vc,vc) = 1.0;
    for(int j = 0; j < X.size1(); j++){
        vector<T> px(Q.size1());
        px = column(X,j);
        vector<T> q = px;
        if(px(j) < 0){
            q(j) = q(j)+(-1)*norm_2(px); 
        }
        else{
            q(j) = q(j)/norm_2(px);
        }
        q = q/norm_2(q);
        
        Q=Q-prod(Q,outer_prod(q,2*trans(q)));
    }
    return Q;
}

//unterlayser for Givensrotation
template<typename T>
matrix<T> makeRM(vector<T> v,int i,int j){
    matrix<T> R(v.size(),v.size());
    for(int vc = 0; vc <v.size(); vc++){
        R(vc,vc) = 1.0;
    }  
        T nn = sqrt(v(i)*v(i)+v(j)*v(j));       
        T c = v(j)/nn;        
        T s = v(i)/nn;
            R(i,i)=c;
            R(j,j)=c;
            R(i,j)=s*(-1);
            R(j,i)=s;
        
    
    return R;
}
//classic Givensrotation returns Q
template <typename T>
matrix<T> Givens(matrix<T> X){
    matrix<T> G(X.size1(),X.size2());
     matrix<T> nG(X.size1(),X.size2());
    for(int vc = 0;vc<X.size1();vc++){
        G(vc,vc)= 1.0;
        nG(vc,vc) =1.0;
    }
    for(int j = 0; j <X.size1()-1; j++){
        vector<T> v (X.size1());
        v = column(X,j);
        for(int i = 1+j; i < X.size1();i++){
            G = makeRM(v,i,j);
            nG = prod(nG,G);
            cout<<nG<<endl;
        }
    }
    return nG;
}


template<typename T>
T matrixNormMax(matrix<T> X){
    matrix<T>mone (X.size1(),X.size1());
    for(int vc = 0; vc < X.size1(); vc++)
    mone(vc,vc) = 1.0;

    matrix<T> mm = prod(trans(X),X);
    mm = mone-mm;
   
    T mx = 0;
    for(int i = 0; i < X.size1(); i++){
        for(int j = 0; j < X.size1(); j++){
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

void TestDouble(){
matrix<double> MFF(3,3);
MFF(0,0) = 124559.0;
MFF(0,1) = 345449.0;
MFF(0,2) = 346559.0;
MFF(1,0) = 122769.0;
MFF(1,1) = 243569.0;
MFF(1,2) = 123449.0;
MFF(2,0) = 122229.0;
MFF(2,1) = 999009.0;
MFF(2,2) = 434999.0;

        matrix<double> cgsm = CGS(MFF);
        double c1 = matrixNormMax(cgsm);
        cgsm = prod(trans(cgsm),cgsm);
        cout<<cgsm<<endl;
        matrix<double> mgsm = MGS(MFF);
          double c2 = matrixNormMax(mgsm);
        mgsm = prod(trans(mgsm),mgsm);  
         cout<<mgsm<<endl;
         
        matrix<double> house = Householder(MFF);
        double c3 = matrixNormMax(house);
        house = prod(trans(house),house);
         cout<<house<<endl;
        matrix<double> g = Givens(MFF);
        double c4 = matrixNormMax(g);
        g = prod(trans(g),g);
         cout<<g<<endl;

 cout<<c1 <<" | "<< c2<<" | " <<c3 <<" | "<< c4<< endl;
 cout<<"*************************************************************************************"<<endl;
cout<<"*************************************************************************************"<<endl; 
}


void Test50(){


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

matrix<cpp_dec_float_50> cgsm = CGS(MM);
        cpp_dec_float_50 c1 = matrixNormMax(cgsm);
        cgsm = prod(trans(cgsm),cgsm);
        cout<<cgsm<<endl;
        matrix<cpp_dec_float_50> mgsm = MGS(MM);
          cpp_dec_float_50 c2 = matrixNormMax(mgsm);
        mgsm = prod(trans(mgsm),mgsm);  
         cout<<mgsm<<endl;
         
        matrix<cpp_dec_float_50> house = Householder(MM);
        cpp_dec_float_50 c3 = matrixNormMax(house);
        house = prod(trans(house),house);
         cout<<house<<endl;
        matrix<cpp_dec_float_50> g = Givens(MM);
        cpp_dec_float_50 c4 = matrixNormMax(g);
        g = prod(trans(g),g);
         cout<<g<<endl;

 cout<<c1 <<" | "<< c2<<" | " <<c3 <<" | "<< c4<< endl; 

 
cout<<"*************************************************************************************"<<endl;
cout<<"*************************************************************************************"<<endl;

}


template<typename T>
matrix<T> clearEps(matrix<T> m){
    for(int i = 0; i < m.size1(); i++){
        for(int j = 0; j <m.size2(); j++){
            if(m(j,i) < 0.000000000000001){
                m(j,i) = 0;
            }
        }
    }
    return m;
}

matrix<cpp_dec_float_100> Hilbert100(unsigned int n){
    matrix<cpp_dec_float_100> m(n,n);
    for(int i = 0; i < m.size1(); i++){
        for(int j = 0; j <m.size2(); j++){            
            cpp_dec_float_100 h = 1.0/(i+j+1.0);
            m(j,i) = h;
        }
    }
    return m;
}
matrix<cpp_dec_float_100> ones100(unsigned int n){
    matrix<cpp_dec_float_100> m(n,n);
    for(int i = 0; i < m.size1(); i++){      
            m(i,i) = cpp_dec_float_100("1");
    }
    return m;
}
matrix<double> HilbertDouble(unsigned int n){
    matrix<double> m(n,n);
    for(int i = 0; i < m.size1(); i++){
        for(int j = 0; j <m.size2(); j++){            
            double h = 1.0/(i+j+1.0);
            m(j,i) = h;
        }
    }
    return m;
}
matrix<double> onesDouble(unsigned int n){
    matrix<double> m(n,n);
    for(int i = 0; i < m.size1(); i++){      
            m(i,i) = 1.0;
    }
    return m;
}


void TestALL(){
    for(int i = 10; i<= 50; i= i+2){
        matrix<cpp_dec_float_100> mm = Hilbert100(i);
        matrix<cpp_dec_float_100> cgsm = CGS(mm);
        matrix<cpp_dec_float_100> mgsm = MGS(mm);
        matrix<cpp_dec_float_100> house = Householder(mm);
        matrix<cpp_dec_float_100> g = Givens(mm);
        cpp_dec_float_100 c1 = matrixNormMax(cgsm);
         cpp_dec_float_100 c2 = matrixNormMax(mgsm);
          cpp_dec_float_100 c3 = matrixNormMax(house);
           cpp_dec_float_100 c4 = matrixNormMax(g);
        cout<<"size: "<<i<<" ||| " <<c1 <<" | "<< c2<<" | " <<c3 <<" | "<< c4<< endl; 
    }
}

void tt(){
matrix<cpp_dec_float_100> mm (100,100); 
 mm = createMatrix100("ss100.txt",mm);
     matrix<cpp_dec_float_100> mgsm = Givens(mm);
cpp_dec_float_100 c2 = matrixNormMax(mgsm);
cout<<mgsm<<endl;   
cout<<"*************************************************************************************"<<endl;
cout<<"*************************************************************************************"<<endl;
mgsm = prod(trans(mgsm),mgsm);
mgsm = clearEps(mgsm);
 cout<<mgsm<<endl;   
cout<<"*************************************************************************************"<<endl;
cout<<"*************************************************************************************"<<endl;
cout<<"*************************************************************************************"<<endl;
cout<<"*************************************************************************************"<<endl;
    cout<<c2<<endl;
}






