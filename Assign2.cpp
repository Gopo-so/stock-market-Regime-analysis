#include<iostream>
#include<sstream>
#include<vector>
#include<cmath>
using namespace std;
 long double pi[2];
 long double A[2][2];
 long double mu[2];
 long double xi[10031][2][2];
 long double gamma_arr[10031][2];
 long double sigma_sq[2];
 long double obs[10031];
 long double alpha[10031][2];
 long double beta_arr[10031][2];

long double density(long double x, int state) {//pdf
    long double coefficient = 1 / sqrt(2 * M_PI * sigma_sq[state]);
    return coefficient * exp(-0.5 * pow((x - mu[state]), 2) / sigma_sq[state]);
}

void HMM(int steps){
    for(int i=0;i<500;i++){//converging by running for 300 times
    //alpha base case
    alpha[0][0] = pi[0]*density(obs[0],0);
    alpha[0][1] = pi[1]*density(obs[0],1);
    
    //alpha gen case
    for(int t=0;t<(steps-1);t++){//t+1 we are calculating
        
        for(int i=0;i<2;i++){
            long double sum1=0.0;
            for(int j=0;j<2;j++){
            sum1+=alpha[t][j]*A[j][i];
            }
            alpha[t+1][i] = density(obs[t+1],i)*sum1;
        }
    }
    //beta_arr base case
    beta_arr[steps-1][0] = 1;
    beta_arr[steps-1][1] = 1;
    //beta_arr gen
   for(int t=steps-2;t>=0;t--){
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){
                beta_arr[t][i] += beta_arr[t+1][j]*A[i][j]*density(obs[t+1],j);
        }
    }
   }
   //gamma_arr 
    for(int t=0;t<steps;t++){
        
        for(int i=0;i<2;i++){
            long double sum2 =0.0;
            for(int j=0;j<2;j++){
                sum2+=alpha[t][j]*beta_arr[t][j];
        }
        gamma_arr[t][i] = alpha[t][i]*beta_arr[t][i]/sum2;
    }

    }
    //xeta
    for(int t=0;t<steps-1;t++){
    long double x_sum =0.0;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            x_sum+=alpha[t][i]*A[i][j]*beta_arr[t+1][j]*density(obs[t+1],j);
        }
    }
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            xi[t][i][j] = alpha[t][i]*A[i][j]*beta_arr[t+1][j]*density(obs[t+1],j)/x_sum;
        }
    }

    }
    //updating parameters
    //pi HMMs by gamma_arr values of first input
    pi[0] = gamma_arr[0][0];
    pi[1] = gamma_arr[0][1];
    //transition matrix,mu,sigma
    long double c_xi[2][2]= {0.0};
    long double c_xid[2]= {0.0};
    long double c_gamma_arr[2]= {0.0};
    long double c_mu[2]= {0.0};
    long double c_sigma[2] = {0.0};
    for(int t=0;t<steps;t++){
        for(int i=0;i<2;i++){
            if(t<steps-1){
                for(int j=0;j<2;j++){
                c_xi[i][j] +=xi[t][i][j];
                } 
                c_xid[i] +=gamma_arr[t][i];
            }
            c_gamma_arr[i] += gamma_arr[t][i];
            c_mu[i]+= obs[t]*gamma_arr[t][i];
            c_sigma[i]+=pow((obs[t] - mu[i]),2)*gamma_arr[t][i];
        }
    }
    for(int i=0;i<2;i++){
        mu[i] = c_mu[i]/c_gamma_arr[i];
        sigma_sq[i] = c_sigma[i]/c_gamma_arr[i];
        for(int j=0;j<2;j++){
            A[i][j] = c_xi[i][j]/c_xid[i];
        }

    }


    }
  
}
void check(int t){
  
    if(mu[0]>mu[1]){
        if(gamma_arr[t][0]>0.5){
            cout<<"Bull"<<endl;
        }
        else{
            cout<<"Bear"<<endl;
        }
    }else{
        if(gamma_arr[t][1]>0.5){
         cout<<"Bull"<<endl;
        }
        else{
            cout<<"Bear"<<endl;
        }
    }
}

int main(){
    mu[0]=mu[1]=0;
    sigma_sq[0]=sigma_sq[1] =1;

   cin>>A[0][0]>>A[0][1]>>A[1][0]>>A[1][1];
   cin>>pi[0]>>pi[1];
   int steps;
   cin>>steps;
   
   for(int i=0;i<steps;i++){
    cin>>obs[i];
   }HMM(steps);
   
    
   for(int i=0;i<steps;i++){
    check(i);
   
   }
   return 0;
}