#include "BGSS+rpl.cpp"
#include <boost/math/special_functions/round.hpp>

template<typename T>
matrix<T> LLLR(matrix<T> X){
    array<matrix<T>,2> aaa;
    double delta = 0.75;
    aaa = cgs_rpl(X,1);
    matrix<T> R = aaa[1];
    matrix<T> Q = aaa[0];
    int k = 1;
    cout<<R<<endl;
    while(k < X.size2()){
        for(int j = k-1; j >= 0; j--){
            T roo =R(j,k)/R(j,j);
            if(roo > 0.5 || roo*(-1) > 0.5){
                cout<<"j = "<<j<<"       roo = "<<roo<<endl;
                roo = boost::math::round(roo);
                column(X,k) = column(X,k) - roo*column(X,j);
                 aaa = cgs_rpl(X,1);
        R = aaa[1];
            }    
        }
        T left = delta*R(k-1,k-1)*R(k-1,k-1);
        T right = R(k,k-1)*R(k,k-1)+R(k,k)*R(k,k);
        cout<<left<<" < "<<right<<endl;
        if(left < right){
            k++;
        }
        else{
            cout << "\nSwap: k = "<<k<< endl;
            vector<T> vv = column(X,k-1);
            column(X,k-1) = column(X,k);
            column(X,k) = vv;
             aaa = cgs_rpl(X,1);
             R = aaa[1];
            if(k > 1){
                k--;
            }
            else{
                k = 1;
            }
        }
    }
    cout<<"----------------------------------------------------------"<<endl;
    cout<<X<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    return X;  
}

template<typename T>
matrix<T> LLLQ(matrix<T> X){
    double delta = 0.75;
    matrix<T> Q = MGS(X);
    int k = 1;

    while(k < X.size2()){
        for(int j = k-1; j >= 0; j--){
            T roo = inner_prod(column(X,k),column(Q,j))/inner_prod(column(Q,j),column(Q,j));
            if(roo > 0.5 || roo*(-1) > 0.5){
                 cout<<"j = "<<j<<"       roo = "<<roo<<endl;
                 roo = round(roo);
            column(X,k) = column(X,k) - (roo*column(X,j));
            Q = MGS(X); 
            }
        }
        T mu = inner_prod(column(X,k),column(Q,k-1))/inner_prod(column(Q,k-1),column(Q,k-1));  
        T left = inner_prod(column(Q,k),column(Q,k));
        T right = (delta-mu*mu)*inner_prod(column(Q,k-1),column(Q,k-1));
        cout<<left<<" > "<<right<<endl;
        if(left >= right){
            k++;
        }
        else{
            cout << "\nSwap: k = "<<k<< endl;
            vector<T> vv = column(X,k-1);
            column(X,k-1) = column(X,k);
            column(X,k) = vv;  
            Q = MGS(X);
            if(k-1 == 0){
                k = 1;
            }
            else{
                k--;
            }
        }
    }
    cout<<"----------------------------------------------------------"<<endl;
    cout<<X<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    return X;  
}

void MatrizenTest(){
matrix<double> V (3,3);

V(0,0)= 1.0;V(0,1)=-1.0;V(0,2)=3.0;
V(1,0)= 1.0;V(1,1)=0.0;V(1,2)=5.0;
V(2,0)= 1.0;V(2,1)=2.0;V(2,2)=6.0;

matrix<double> B (3,3);

B(0,0)= 15.0;B(0,1)=46.0;B(0,2)=32.0;
B(1,0)= 23.0;B(1,1)=15.0;B(1,2)=1.0;
B(2,0)= 11.0;B(2,1)=3.0;B(2,2)=1.0;

matrix<double> N (2,2);

N(0,0)= 201.0;N(0,1)=1648.0;
N(1,0)= 37.0;N(1,1)=297.0;

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

matrix<double> M(6,6);
M(0,0) = 35051.0; M(0,1) = 0;   M(0,2) = 0;   M(0,3) = 0;   M(0,4) = 0;   M(0,5) = 0;
M(1,0) = 1669.0 ; M(1,1) = 1.0; M(1,2) = 0;   M(1,3) = 0;   M(1,4) = 0;   M(1,5) = 0;
M(2,0) = 4828.0 ; M(2,1) = 0;   M(2,2) = 1.0; M(2,3) = 0;   M(2,4) = 0;   M(2,5) = 0;
M(3,0) = 33763.0; M(3,1) = 0;   M(3,2) = 0;   M(3,3) = 1.0; M(3,4) = 0;   M(3,5) = 0;
M(4,0) = 4875.0 ; M(4,1) = 0;   M(4,2) = 0;   M(4,3) = 0;   M(4,4) = 1.0; M(4,5) = 0;
M(5,0) = 32724.0; M(5,1) = 0;   M(5,2) = 0;   M(5,3) = 0;   M(5,4) = 0;   M(5,5) = 1.0; 
 
matrix<double> R(5,5);
R(0,0) = 1000.0; R(0,1) = 0.0; R(0,2) = 0.0; R(0,3) = 0.0; R(0,4) = 0.0; 
R(1,0) = 2222.0 ; R(1,1) = 1.0; R(1,2) = 0.0; R(1,3) = 0.0; R(1,4) = 0.0; 
R(2,0) = 3333.0 ; R(2,1) = 0.0; R(2,2) = 1.0; R(2,3) = 0.0; R(2,4) = 0.0; 
R(3,0) = 4000.0; R(3,1) = 0.0; R(3,2) = 0.0; R(3,3) = 1.0; R(3,4) = 0.0; 
R(4,0) = 5555.0 ; R(4,1) = 0.0; R(4,2) = 0.0; R(4,3) = 0.0; R(4,4) = 1.0;
}







