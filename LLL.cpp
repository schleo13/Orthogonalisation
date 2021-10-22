#include "BGSS+rpl.cpp"
#include <boost/math/special_functions/round.hpp>

/*
Berechnet durch den naiven LLL in Blockschreibweise die reduzierte Matrix.
Diese Methode wird durch die Matrix R definiert.
*/
template<typename T>
matrix<T> LLLBlock(matrix<T> X,double delta, int bsize){
    auto start = std::chrono::high_resolution_clock::now();
   
    array<matrix<T>,2> aaa;
    aaa = bcgs_rpl(X,bsize);
    matrix<T> R = aaa[1];
    int k = 1;
    int swap = 0;
     range ee(0,X.size1());
     int kk = bsize;
     int blocksize = bsize;
     matrix<T> M;
    
    while(k < X.size2()){
     for(int j = k-1; j >= 0; j--){
                T roo =R(j,k)/R(j,j);
            if(roo >= 0.5 || roo*(-1) >= 0.5 ){
                roo = boost::math::round(roo);
                column(X,k) = column(X,k) - (roo*column(X,j));
                column(R,k) = column(R,k) - (roo*column(R,j));
            }
     }    
        T left = delta*R(k-1,k-1)*R(k-1,k-1);
       T right = (R(k,k-1)*R(k,k-1))+(R(k,k)*R(k,k));
      
        if(left <= right){
            k++;
            if(k == kk-1 && kk+blocksize <= X.size2()){
                 kk = kk +blocksize;
            }
        }
        else{
         
           swap++;
           if(k > blocksize && k == kk-blocksize-1){
               kk = kk-blocksize;
           }
            column(X,k).swap(column(X,k-1));
             range rk (0,kk);
             M = project(X,ee,rk);
             aaa = bcgs_rpl(M,bsize);
             project(R,ee,rk) = aaa[1];
            
            if(k > 1){
                k--;
            }
            else{
                k = 1;
            }   
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_ms = end - start;
    cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;    
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"Anzal swaps = "<<swap<<"      NormMAX = "<<matrixNormMax(X)<<"   NormMIN = "<<matrixNormMin(X)<<endl;
     cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
    return X;  
}



/*
Berechnet die Reduziert mit naiven LLL definiert durch R.
*/
template<typename T>
matrix<T> LLLR(matrix<T> X,double delta){
     auto start = std::chrono::high_resolution_clock::now();
    array<matrix<T>,2> aaa;
    aaa = MGS(X);
    matrix<T> R = aaa[1];
    int k = 1;
    int swap = 0;
   
    while(k < X.size2()){
       
     for(int j = k-1; j >= 0; j--){
                T roo =R(j,k)/R(j,j);
            if(roo >= 0.5 || roo*(-1) >= 0.5 ){
                roo = boost::math::round(roo);
                column(X,k) = column(X,k) - (roo*column(X,j));
                column(R,k) = column(R,k) - (roo*column(R,j));
                
            }
     }       
       T left = delta*R(k-1,k-1)*R(k-1,k-1);
       T right = (R(k,k-1)*R(k,k-1))+(R(k,k)*R(k,k));
      
        if(left <= right){
            k++;
        }
        else{
           swap++;
            column(X,k).swap(column(X,k-1));
            aaa = MGS(X);
            R = aaa[1];
            if(k > 1){
                k--;
            }
            else{
                k = 1;
            }   
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_ms = end - start;
    cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;    
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"Anzal swaps = "<<swap<<"          NormMAX = "<<matrixNormMax(X)<<"   NormMIN = "<<matrixNormMin(X)<<endl;
     cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
    return X;  
}


/*
Berechnet die Reduziert mit naiven LLL mit festen Datentypen.
*/
matrix<mpz_int> LLLRMPZ(matrix<mpz_int> X,double delta){
     auto start = std::chrono::high_resolution_clock::now();
    array<matrix<mpff>,2> aaa;
    matrix<mpff> A = X;
    aaa = cgs_rpl(A);
    matrix<mpff,column_major> R = aaa[1];
    int k = 1;
    int swap = 0;
    while(k < X.size2()){
     for(int j = k-1; j >= 0; j--){
                mpz_int roo (boost::math::round(R(j,k)/R(j,j)));
            if(roo >= 1|| roo*(-1) >= 1 ){ 
                column(X,k) = column(X,k) - (roo*column(X,j));
                column(R,k) = column(R,k) - (roo*column(R,j));
                
            }
     }       
       mpff left = delta*R(k-1,k-1)*R(k-1,k-1);
       mpff right = (R(k,k-1)*R(k,k-1))+(R(k,k)*R(k,k));
      
        if(left <= right){
            k++;
            
        }
        else{
           swap++;
            column(X,k).swap(column(X,k-1));
            A = X;
            aaa = cgs_rpl(A);
            R = aaa[1];
            if(k > 1){
                k--;
            }
            else{
                k = 1;
            }   
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_ms = end - start;
    cout<<"Benötigte Zeit : "<<float_ms.count()<<endl;    
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"Anzal swaps = "<<swap<<"          NormMAX = "<<matrixNormMax(X)<<"   NormMIN = "<<matrixNormMin(X)<<endl;
     cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
    return X;  
}





/*
Berechnet die Reduziert mit naiven LLL definiert durch Q.
Orthogonalisierungsverfahren darf kein Orthonormalverfahren sein.
*/
template<typename T>
matrix<T> LLLQ(matrix<T> X){
    double delta = 0.75;
    array<matrix<T>,2> aaa;
    aaa = MGS(X);
    matrix<T> Q = aaa[0];
    int k = 1;
    int swap = 0;
    while(k < X.size2()){
        for(int j = k-1; j >= 0; j--){
            T roo = inner_prod(column(X,k),column(Q,j))/inner_prod(column(Q,j),column(Q,j));
            if(roo >=0.5 || roo*(-1) >= 0.5){
                // cout<<"j = "<<j<<"       roo = "<<roo<<endl;
                  roo = boost::math::round(roo);
                  //cout<<roo<<endl;
            column(X,k) = column(X,k) - roo*column(X,j);
            aaa = MGS(X);
            Q = aaa[0]; 
            }
        }
        T mu = inner_prod(column(X,k),column(Q,k-1))/inner_prod(column(Q,k-1),column(Q,k-1));  
        T left = inner_prod(column(Q,k),column(Q,k));
        T right = (delta-mu*mu)*inner_prod(column(Q,k-1),column(Q,k-1));
        //cout<<left<<" > "<<right<<endl;
        if(left > right){
            k++;
        }
        else{
           // cout << "\nSwap: k = "<<k<< endl;
           swap++;
           
            vector<T> vv = column(X,k-1);
            column(X,k-1) = column(X,k);
            column(X,k) = vv;
            aaa = MGS(X);
            Q = aaa[0]; 
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
    //cout<<matrixNormMax(X)<<endl;
    cout<<swap<<endl;
    return X;  
}


void LLLTest1(){
int   dim = 40;
string s = "svpchallengedim40seed0.txt";
mpff::default_precision(130);
matrix<mpff> A (dim,dim);
A = createMatrix(s,A);
A = trans(A);
LLLR(A,0.75);
LLLBlock(A,0.75,5);
LLLBlock(A,0.75,10);
LLLBlock(A,0.75,20);
}

void LLLTest2(){
int   dim = 40;
string s = "svpchallengedim40seed0.txt";
//mpff::default_precision(130);
matrix<mpff> A (dim,dim);
A = createMatrix(s,A);
A = trans(A);

for(double i = 0.5; i <= 0.76; i = i + 0.05)
LLLBlock(A,i,5);

}

void LLLTest3(){
int   dim = 10;
string s = "matrix100.txt";
matrix<mpff> A (dim,dim);
A = createMatrix(s,A);
A = trans(A);


LLLR(A,0.75);
LLLBlock(A,0.75,5);


}

void LLLTest4(){
int   dim = 40;
string s = "svpchallengedim40seed0.txt";
mpff::default_precision(130);
matrix<mpff> A (dim,dim);
A = createMatrix(s,A);
A = trans(A);

LLLR(A,0.75);
LLLBlock(A,0.75,5);


}








