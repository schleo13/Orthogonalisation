#include "BGSS+rpl.cpp"
#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/pow.hpp>

template<typename T>
matrix<T> LLL2(matrix<T> X){
    array<matrix<T>,2> aaa;
    double delta = 3.0/4.0;
    aaa = cgs_rpl(X,1);
    matrix<T> R = aaa[1];
    int k = 1;
    while(k < X.size2()){
         
        for(int j = k-1; j >= 0; j--){
            //T roo = round(inner_prod(column(X,k),column(M,j))/inner_prod(column(M,j),column(M,j)));
            T roo = boost::math::round(R(j,k)/R(j,j));
            if(roo != 0){
                cout<<"j = "<<j<<"       roo = "<<roo<<endl;
            column(X,k) = column(X,k) - roo*column(X,j);
            aaa = cgs_rpl(X,1);
            R = aaa[1]; 
            }
              
        }
        //T mu = inner_prod(column(X,k),column(M,k-1))/inner_prod(column(M,k-1),column(M,k-1));  
       // T mu = R(k-1,k)/R(k-1,k-1); cout<<" mu = "<<mu<<" ; ";
        //T left = inner_prod(column(M,k),column(M,k));
        T left = delta*R(k-1,k-1)*R(k-1,k-1);
        //T right = ((3.0/4.0)-mu*mu)*inner_prod(column(M,k-1),column(M,k-1));
        //T right = (delta-mu*mu)*R(k-1,k-1)*R(k-1,k-1);
        T right = R(k,k-1)*R(k,k-1)+R(k,k)*R(k,k);
        cout<<left<<" > "<<right<<endl;
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


int main(){
    
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




cout<<"***********************************************************************************************"<<endl;

 cout<<"***********************************************************************************************"<<endl; 
cout<<"***********************************************************************************************"<<endl;
cout<<"***********************************************************************************************"<<endl;


matrix<cpp_dec_float_100> HARD(60,60);
 //HARD = createMatrix100("svp60.txt",HARD);
 matrix<cpp_dec_float_100> AA = LLL2(V);
matrix<cpp_dec_float_100> BB = LLL2(N);

matrix<cpp_dec_float_100> CC = LLL2(MM);


cout<<"***********************************************************************************************"<<endl; 
cout<<"***********************************************************************************************"<<endl;
cout<<"***********************************************************************************************"<<endl;



cout<<"***********************************************************************************************"<<endl; 
cout<<"***********************************************************************************************"<<endl;
cout<<"***********************************************************************************************"<<endl;
//matrix<cpp_dec_float_100> HH = LLL2(HARD);



    return 0;
}
