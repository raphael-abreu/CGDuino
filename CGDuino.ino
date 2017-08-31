// see online at https://circuits.io/circuits/3034264-the-unnamed-circuit/
// Raphael Abreu, Out 2016

#define precision 0.0001
#define decimal_places 5

void setup() {
  Serial.begin(9600);      // open the serial port at 9600 bps:  
}

float vectors_dot_prod(const float *x, const float *y, int n){
  float res = 0;
  int i;
  for (i = 0; i < n; i++){
      res += x[i] * y[i];
  }
  return res;
}


float pattern_dot_prod(int i, const float *y,const float *pattern, int patternSize, int n){
  int k =0,j=0;
  float res =0;
   for(j = i,k=floor(patternSize/2);j<n && k<patternSize;j++,k++){
  		res += y[j]*pattern[k];
  }
   for(j = i-1,k=floor(patternSize/2)-1;j>=0 && k>=0;j--,k--){
  		res += y[j]*pattern[k];
  }
  return res;
}



float getMax(float* array, int size){
 int max = array[0];
 for (int i=1; i<size; i++){
   if (max<array[i]){
     max = array[i];
   }
 }
     return max;

}  

float absoluto(float x){
if (x>=0)
    return x;
if (x <0)
    return x*-1;
}

float eval(float x,float y){
float temp = x-y;
if (absoluto(temp)>precision)
    return 1;
else
    return 0;
} 
 

void loop() {  
  
/*

#This is a MATRIX FREE conjugate gradient method, wich means that you don't have to compute a matrix X. instead you have function or, a pattern that is repeated on the diagonal (N diagonal sparse matrices).
#As an example, this method can be used to solve the discrete poisson matrix for heat propagation: https://en.wikipedia.org/wiki/Discrete_Poisson_equation
#Have fun!

*/

float pattern[7]{-1,0,-1,4,-1,0,-1};
float b[9] = {2, 1, 5, 2, 1, 5, 2, 1, 5};
float x[9] = {0,0,0,0,0,0,0,0,0};
float tes[6] = {5,2,1,3,4,6};
int patternSize=7;
  
int count = 1,i,j,k,n=9;
float result,passo,B,rMax,r_oldMax;

float Ab[n];
float Ax[n];
float Ap[n];
float r[n];
float r_old[n];
float p[n];

 Serial.print("\n");
Serial.print("\n");
Serial.print("\n");

  

// fake_MatVec
for (i = 0; i < n; i++){
     Ax[i] = pattern_dot_prod(i, x,  pattern, patternSize, n);
} 

  
// r = b - Ax;
for (i = 0; i < n; i++){
    r[i] = b[i]-Ax[i]; //residual
    p[i] = r[i]; // in the first iteraction, p is the residual
}

do{ 
Serial.print("*******");Serial.print(count);Serial.print("*******\n"); 
  

 
  // fake_MatVec
for (i = 0; i < n; i++){
     Ap[i] = pattern_dot_prod(i, p,  pattern, patternSize, n);
} 




 // direction = (rT[i]* r[i])/(pT[i])* Ap[i]
passo = (vectors_dot_prod(r,r,n) / vectors_dot_prod(p,Ap,n));



// update search direction
for (i = 0; i < n; i++){
    x[i] = x[i] + (passo*p[i]);
}


//compute the new residual
for (i = 0; i < n; i++){
   r_old[i] = r[i]; // saving for use below
   r[i] = r[i]-(passo*Ap[i]);
  
}


//B = rT*r / old_rT*old_r
B = (vectors_dot_prod(r,r,n) / vectors_dot_prod(r_old,r_old,n));
Serial.print("Beta:"); 
Serial.print(B,decimal_places); 
Serial.print("\n");


//p = r + b*p
for (i = 0; i < n; i++){
p[i] = r[i] + B*p[i];
}
for(i=0;i<n;i++){
Serial.print("X[");
Serial.print(i);
Serial.print("]:");
Serial.print(x[i],decimal_places);
Serial.print("\n");

}
Serial.print("-------------\n");
 rMax = getMax(r,n);
 r_oldMax = getMax(r_old,n);
  count++; 
  }while(B>precision);
  
  
  delay(20000);
  Serial.println(""); 
}

