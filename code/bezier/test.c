//g++ -o test test.c -g  -Wall -W -ansi -Dcimg_use_vt100 -Dcimg_display=0 -fno-tree-pre -lm
//sfdp graph2.dot -Tpdf -o a.pdf -Goverlap=prism -Goverlap_scaling=2 -Gsep=+40 -Gsplines
// neato graph2.dot -Tpdf -o a.pdf

//ffmpeg -i xxx.avi -f image2 -vf scale=iw/8:-1,fps=fps=0.05 img%03d.jpg

  //montage img* -tile 1x -geometry +0+0 test.jpg

// ls graphx1.dot | while read a; do sfdp $a -Tpdf -o $a.pdf -Goverlap=prism1000 -Goverlap_scaling=2 -Gsep=+10; open $a.pdf; done

#include "CImg.h"
#include <vector>
#include <algorithm>
#include "fheap.c"
#include <stdlib.h>
#include <string.h>
#define NOMEMAX 100
#define S 16
#define NB 256/S

#define INFTY 100000
// quante linee di statistica sulla immagine originale
#define divisioni 100

using namespace cimg_library;
using namespace std;

CImg<unsigned char> stro1;
CImg<unsigned char> stro2;
CImg<unsigned char> src;

struct edge{
  int x,y,z,flag;
  float w;
  vector<int> nedge; // codice edge successivo
  vector<float> wedge; // peso edge successivo
  
};
struct char3{
  unsigned char x,y,z;
};

long A[NB][NB][NB];
long B[NB][NB][NB];
long C[NB][NB][NB];
long D[NB][NB][NB];
float E[NB][NB][NB];
vector<edge> alist[NB][NB][NB];
int gnode[NB][NB][NB];
vector<char3> gchar3;

float pos[NB][NB][NB][2];

int code=1;


CImg<double> stat2;


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);


#define PI 3.1415


int rec(int x, int y, int z){

  if (C[x][y][z]>0)
    return C[x][y][z];

  if (B[x][y][z]==0)
    return 0;

  long best=0;
  int nx,ny,nz;
  for (int dx=-1;dx<=1;dx++)
    for (int dy=-1;dy<=1;dy++)
      for (int dz=-1;dz<=1;dz++)
	if (x+dx>=0 && y+dy>=0 && z+dz>=0 &&
	    x+dx<NB && y+dy<NB && z+dz<NB){
	  if (B[x+dx][y+dy][z+dz]>best){
	    best=B[x+dx][y+dy][z+dz];
	    nx=x+dx;
	    ny=y+dy;
	    nz=z+dz;
	  }		  
	}
  if (best>0 && nx==x && ny==y && nz==z){ // sono il max
    if (C[x][y][z]==0){
      C[x][y][z]=code++;
      return C[x][y][z];
    }
    else
      return C[x][y][z];
  }
  //printf("%d %d %d -> %d %d %d\n",x,y,z,nx,ny,nz);
  //non massimo
  int co=-1;
  if (C[nx][ny][nz]>0){ // gia' assegnato
    co=C[nx][ny][nz];
  }
  else
    co=rec(nx,ny,nz);
  C[x][y][z]=co;
  return co;
}


vector<int> label;
vector<float> w;  //peso di ogni nodo

vector<int> code_map; /// per ogni codice dice il nodo massimo

int rec_label(int x, int y, int z){

  if (C[x][y][z]>=0)
    return C[x][y][z];

  if (B[x][y][z]==0)
    return 0;

  int st1=gnode[x][y][z]; // non ottimizzato (posso gia lavorare con xyz)
  int dbg=(st1==3026)||(st1==3027);
  if (dbg)
  printf("%d %d %d: node %d\n",x,y,z,st1);
  char3 no=gchar3[st1];

  if (dbg)
  printf("%d %d %d: node %d check %d %d %d\n",x,y,z,st1,no.x,no.y,no.z);
  int found=-1;
  char3 no1;
  no1.x=-1;
  
  for (int ct1=0;ct1<alist[no.x][no.y][no.z].size();ct1++){ // per ogni arco
    int neigh=gnode[alist[no.x][no.y][no.z][ct1].x][alist[no.x][no.y][no.z][ct1].y][alist[no.x][no.y][no.z][ct1].z];
    
  if (dbg)
    printf("w %f -> %f arco %f (%d -> %d), %d %d %d\n",w[st1],w[neigh],alist[no.x][no.y][no.z][ct1].w,st1,neigh,no1.x,no1.y,no1.z);
    if (w[neigh]>w[st1]){ // solo monotono
      // scelgo miglior candidato
      if (found==-1){
	found=neigh;//ct1
	no1=gchar3[neigh];
      }
      else{
	//	if (alist[no.x][no.y][no.z][ct1].w>alist[no.x][no.y][no.z][found].w){
	if (w[neigh]>w[found]){
	  found=neigh;//ct1;
	  no1=gchar3[neigh];
	}
      }
    }
  }
  if (dbg)
  if (found>=0)
  printf("found %d (%d %d %d), code %d, code next %d\n",found,no1.x,no1.y,no1.z,code,C[no1.x][no1.y][no1.z]);
  else
    printf("not found\n");

  if (found<0){ // sono il max
    if (C[x][y][z]==-1){
      printf("Max %d: code %d\n",st1,code);
      code_map.push_back(st1);
      C[x][y][z]=code++;
      return C[x][y][z];
    }
    else
      return C[x][y][z];
  }
  //printf("%d %d %d -> %d %d %d\n",x,y,z,nx,ny,nz);
  //non massimo
  int co=-1;
  if (C[no1.x][no1.y][no1.z]>=0){ // gia' assegnato    
    co=C[no1.x][no1.y][no1.z];
  }
  else
    co=rec_label(no1.x,no1.y,no1.z);
  C[x][y][z]=co;
  return co;
}

struct float4{
  float x,y,z,w;
};
  struct st{
    long m1,M1;
    long m2,M2;
    float best; // best val di interscambio
    float col;
    int idx;
    float acc;
    float w;

    float arco;

    int node,ed;
    
    int x1,y1,z1;
    int x2,y2,z2;
  };
struct int3{
  int x,y,z;
};
struct int2{
  int x,y;
};
struct float2{
  float x,y;
};

void rgb2hsl(float r, float g, float b, float& h, float& s, float& l){
  float maxf = fmax(fmax(r, g), b);
  float minf = fmin(fmin(r, g), b);
  l = (maxf + minf) / 2;

    if(maxf == minf){
        h = s = 0; // achromatic
    }else{
      float d = maxf - minf;
        s = (l > 0.5f) ? d / (2.0f - maxf - minf) : d / (maxf + minf);

        if (r > g && r > b)
            h = (g - b) / d + (g < b ? 6.0f : 0.0f);

        else if (g > b)
            h = (b - r) / d + 2.0f;

        else
            h = (r - g) / d + 4.0f;

        h /= 6.0f;
    }
}

float norm(float* a){
 return pow(a[0]*a[0]+a[1]*a[1]+a[2]*a[2],0.5f);
}
void normal(float* a){
  float norm=pow(a[0]*a[0]+a[1]*a[1]+a[2]*a[2],0.5f);
  if (norm>0){
    a[0]/=norm;
    a[1]/=norm;
    a[2]/=norm;
  }
}


vector<int> node;
vector<float> noder;
vector<float2> nodex;
vector< vector <float > > adj;
vector< vector <int > > ord;
vector<vector<int> > ordn; //vicini
vector<vector<int> > ordp; //permutazione corrente
vector<int> ciclo;
vector< vector <int2> > stack;

vector< vector<int > > dep; // per ogni code, metto nodi che sono sotto
vector< vector<int > > depi; // per ogni code, metto nodi che sono sopra

vector<float4> col;
float max1=0;

float sizen(int i){
  float temp=pow(col[i].w/max1,0.5f);
  return 10*temp+(1-temp)*0.5;
}

void rec_a(int i){
  //assegna seguenti
  int dbg=1;
  for (int j=0;j<code;j++){
    if (adj[i][j]>0){
      if (noder[j]==1000 || noder[j]<noder[i]+sizen(i)){
	if (dbg) printf("noderA %d (%d) %f\n",i,j,noder[i]);
	noder[j]=noder[i]+sizen(i);
	if (dbg) printf("noderB %d %f\n",i,noder[i]);
	rec_a(j);
      }
      if (noder[j]!=1000 && noder[j]>noder[i]+sizen(i)){ // provo a ridurre
	int highest=-1000;
	for (int j1=0;j1<code;j1++){
	  if (adj[j1][j]>0) //calcolo minimo entranti
	    if (noder[j1]!=1000)
	      if (highest<noder[j1]+sizen(j1))
		highest=noder[j1]+sizen(j1);
	}
	if (highest<noder[j]){
	  //	  printf("lift %d -> %f\n",j,highest);
	  if (dbg) printf("noderC %d (%d) %f\n",j,i,noder[j]);
	  noder[j]=highest;
	  if (dbg) printf("noderD %d %f\n",j,noder[j]);
	  rec_a(j);
	}
      }
    }
    if (adj[j][i]>0){
      if (noder[j]==1000 || noder[j]>noder[i]-sizen(j)){
	if (dbg) printf("noderE %d (%d) %f\n",j,i,noder[j]);
	noder[j]=noder[i]-sizen(j);
	if (dbg) printf("noderF %d %f\n",j,noder[j]);
	rec_a(j);
      }
    }
  }
}

void rec_x(int i){
  //posizione x
  int us=0;
  for (int j=0;j<code;j++)
    if (adj[i][j]>0)
      us++;

  for (int j=0;j<code;j++){
    if (adj[i][j]>0){
      if (nodex[j].x==-1){
	noder[j]=noder[i]+sizen(i);
	rec_a(j);
      }
      if (noder[j]!=1000 && noder[j]>noder[i]+sizen(i)){ // provo a ridurre
	int highest=-1000;
	for (int j1=0;j1<code;j1++){
	  if (adj[j1][j]>0) //calcolo minimo entranti
	    if (noder[j1]!=1000)
	      if (highest<noder[j1]+sizen(j1))
		highest=noder[j1]+sizen(j1);
	}
	if (highest<noder[j]){
	  //	  printf("lift %d -> %f\n",j,highest);
	  noder[j]=highest;
	  rec_a(j);
	}
      }
    }
    if (adj[j][i]>0){
      if (noder[j]==1000 || noder[j]>noder[i]-sizen(j)){
	noder[j]=noder[i]-sizen(j);
	rec_a(j);
      }
    }
  }
}

int rec_v(int i){
  //  printf("visit %d: %d\n",i,node[i]);
  if (node[i]==0) {
    ciclo.push_back(i);
    return -1;
  }
  if (node[i]==-1){
    node[i]=0;
    int res=0;
    for (int k=0;k<code;k++)
      if (adj[i][k]>0){
	if (rec_v(k)==-1){
	  res=-1;
	  k=code; //kill
	}
      }
    if (res==-1)
      ciclo.push_back(i);
    node[i]=1;	
    return res;
  }
  return 0;
}

vector< vector<int> > lista1;
vector< vector<int> > lista2;

inline float compare(int i,int k){
  float cti=0;
  int i1,i2,k1,k2;
  i1=lista1[i][0];
  i2=lista1[i][lista1[i].size()-1];
  k1=lista1[k][0];
  k2=lista1[k][lista1[k].size()-1];
  cti=-fabs(
    (noder[i2]+noder[i1])/2
    -
    (noder[k2]+noder[k1])/2);
  return cti;

  for (int i1=0;i1<lista1[i].size();i1++)
    for (int k1=0;k1<lista1[k].size();k1++)
      if (lista1[i][i1]==lista1[k][k1])
	cti+=1;//*fabs(lista1w[i]-lista1w[k])/(lista1w[i]+lista1w[k])*fabs(noder[lista1[i][i1]]-noder[lista1[k][k1]]);//*(col[lista1[i][i1]].w);
  return cti;
}

int transitiva(int nextp){
  int ok=1;
  int mod=1;
  while (mod){
    mod=0;
    for (int i=0;i<code;i++)
      for (int j=0;j<code;j++)
	if (abs(ord[i][j])==1)
	  for (int k=0;k<code;k++){
	    if (ord[i][j]==ord[j][k] &&
		ord[i][k]==0 && ord[k][i]==0){
	      ord[i][k]=ord[i][j];
	      ord[k][i]=-ord[i][j];
	      //printf("addT %d %d: %d\n",i,k,ord[i][j]);
	      int2 t;
	      t.x=i;
	      t.y=k;
	      stack[nextp].push_back(t);
	      mod=1;
	    }
	    if (ord[i][j]==ord[j][k] &&
		(ord[i][k]==-ord[i][j] || ord[k][i]==ord[i][j])){
	      printf("FAILT-- %d %d\n",i,k);
	      ok=0;
	    }

	  }
  }
  return ok;
}

int test_verticale(int& r1, int& r2){
  int ok=1;

  // per ogni coppia vertici
  // per ogni arco uscente
  for (int i=0;i<code;i++)
    for (int j=0;j<code;j++)
      if (i!=j && 
	  noder[i]<noder[j] && 
	  abs(ord[i][j])!=2 && abs(ord[i][j])!=0){ // conosco la precedenza
	for (int i1=0;i1<code;i1++){ // con lista ad si velocizza
	  if (ord[i][i1]==2)
	    for (int j1=0;j1<code;j1++){ // con lista ad si velocizza
	      if (ord[j][j1]==2)
		if (abs(ord[i1][j1])!=2 && abs(ord[i1][j1])!=0){//conosco la precedenza

		  //check
		  if (noder[i1]>noder[j] && // c'e' possibile intersezione y
		      ord[i][j]!=ord[i1][j1]
		      ){
		    int a=i1;
		    int b=j1;
		    // da togliere
		    int inters=0;
		    int best=-1;
		    float bestv=1000;
		    for (int p1=0;p1<dep[a].size();p1++){
		      int found=0;
		      for (int p2=0;p2<dep[b].size();p2++){
			if (dep[a][p1]==dep[b][p2]){
			  found=1;
			  if (bestv>noder[dep[a][p1]]){
			    bestv=noder[dep[a][p1]];
			    best=dep[a][p1];
			    inters=1;
			  }
			    
			}
		      }
		    }
		    printf("highest int %d, %d\n",inters,best);
		    if (0&& inters==0){
		      a=i;   // interseco da sopra (piu' impreciso)
		      b=j;
		      // da togliere
		      best=-1;
		      bestv=1000;
		      for (int p1=0;p1<dep[a].size();p1++){
			int found=0;
			for (int p2=0;p2<dep[b].size();p2++){
			  if (dep[a][p1]==dep[b][p2]){
			    found=1;
			    if (bestv>noder[dep[a][p1]]){
			      bestv=noder[dep[a][p1]];
			      best=dep[a][p1];
			      inters=1;
			    }			    
			  }
			}
		      }
		      printf("highest int %d, %d\n",inters,best);
		    }


		      if (inters==0){
			printf("indip\n");
		      }
		      else{
			printf("FAILC-- %d %d %d %d (%f %f %f %f, %d %d)\n",i,j,i1,j1,noder[i],noder[j],noder[i1],noder[j1],ord[i][j],ord[i1][j1]);

			int c=best;
		    
			vector<int> l1;

			inters=0;
			for (int p1=0;p1<dep[a].size();p1++){
			  int found=0;
			  for (int p2=0;p2<depi[c].size();p2++){
			    if (dep[a][p1]==depi[c][p2]){
			      found=1;
			      p2=depi[c].size();
			    }
			  }
			  if (found) {
			    l1.push_back(dep[a][p1]);
			    printf("%d\n",dep[a][p1]);
			  }
			}

			for (int p1=0;p1<dep[b].size();p1++){
			  int found=0;
			  for (int p2=0;p2<depi[c].size();p2++){
			    if (dep[b][p1]==depi[c][p2]){
			      found=1;
			      p2=depi[c].size();
			    }
			  }
			  if (found) {
			    l1.push_back(dep[b][p1]);
			    printf("%d\n",dep[b][p1]);
			  }
			}
			printf("hi\n");
			int bx=-1,by=-1;
			float bv=10e10;
			for (int i=0;i<l1.size();i++)
			  for (int j=0;j<l1.size();j++)
			    if (adj[l1[i]][l1[j]]>0){
			      printf("sel %d %d %f\n",l1[i],l1[j],adj[l1[i]][l1[j]]);
			      if (bv>adj[l1[i]][l1[j]]){
				bx=l1[i];
				by=l1[j];
				bv=adj[l1[i]][l1[j]];
			      }
			    }
			printf("remove %d %d %f\n",bx,by,bv);		    
			r1=bx;r2=by;
			ok=0;				
		      }
		  }
		}
	    }
	}
      }

  return ok;
}

float ell_d[3]; // autovalori
float** ell_v;  // autovettori
float ell_c[3]; // baricentro

void jacobi(float **a, int n, float d[], float **v, int *nrot)
// Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
// output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
// v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
// a. nrot returns the number of Jacobi rotations that were required.
{
 int j,iq,ip,i;
 float tresh,theta,tau,t,sm,s,h,g,c;
 vector<float> b;
 vector<float> z;
 b.resize(n+1);
 z.resize(n+1);
 for (ip=1;ip<=n;ip++) 
  {// Initialize to the identity matrix.
   for (iq=1;iq<=n;iq++) 
    v[ip][iq]=0.0;
    
   v[ip][ip]=1.0;
  }
  
for (ip=1;ip<=n;ip++) 
{// Initialize b and d to the diagonal of a. 
  b[ip]=d[ip]=a[ip][ip];
  z[ip]=0.0; //This vector will accumulate terms of the form tapq as in equation (11.1.14).
}
*nrot=0;
for (i=1;i<=50;i++) {
 sm=0.0;
 for (ip=1;ip<=n-1;ip++) { //Sum off-diagonal elements.
 for (iq=ip+1;iq<=n;iq++)
 sm += fabs(a[ip][iq]);
 }

if (sm == 0.0) { //The normal return, which relies on quadratic convergence to machine underflow.
  return;
}

if (i < 4)
  tresh=0.2*sm/(n*n); //...on the first three sweeps.
else
  tresh=0.0; //...thereafter.
for (ip=1;ip<=n-1;ip++) {
 for (iq=ip+1;iq<=n;iq++) {
 g=100.0*fabs(a[ip][iq]); //After four sweeps, skip the rotation if the off-diagonal element is small.
 if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
     && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
  a[ip][iq]=0.0;
 else if (fabs(a[ip][iq]) > tresh) {
   h=d[iq]-d[ip];
   if ((float)(fabs(h)+g) == (float)fabs(h))
    t=(a[ip][iq])/h; //t = 1/(2.)
   else {
    theta=0.5*h/(a[ip][iq]); //Equation (11.1.10).
    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
    if (theta < 0.0) t = -t;
   }
 c=1.0/sqrt(1+t*t);
 s=t*c;
 tau=s/(1.0+c);
 h=t*a[ip][iq];
 z[ip] -= h;
 z[iq] += h;
 d[ip] -= h;
 d[iq] += h;
 a[ip][iq]=0.0;
 for (j=1;j<=ip-1;j++) { //Case of rotations 1 = j < p.
   ROTATE(a,j,ip,j,iq)
 }
 for (j=ip+1;j<=iq-1;j++) { //Case of rotations p < j < q.
   ROTATE(a,ip,j,j,iq)
 }
 for (j=iq+1;j<=n;j++) { //Case of rotations q < j = n.
   ROTATE(a,ip,j,iq,j)
 }
 for (j=1;j<=n;j++) {
   ROTATE(v,j,ip,j,iq)
 } 
 ++(*nrot);
}
}
}

for (ip=1;ip<=n;ip++) {
 b[ip] += z[ip];
 d[ip]=b[ip]; //Update d with the sum of tapq,
 z[ip]=0.0; //and reinitialize z.
}
}
printf("Too many iterations in routine jacobi");
}









void stats(){
  printf("welcome\n");
  int sizey=400;
  CImg<int> stat1(1024,NB*NB,1,3,255);



  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++)
	E[x][y][z]=0;

  vector<int3> colp;
  for (int i=0;i<code;i++){
    float4 temp;
    int3 temp1;
    temp1.x=-1;
    col.push_back(temp);
    colp.push_back(temp1);
    col[i].x=0;
    col[i].y=0;
    col[i].z=0;
    col[i].w=0;
  }
  float max2=0;
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
	if (C[x][y][z]>0){
	  if (max2<B[x][y][z]){
	    max2=B[x][y][z];
	  }
	  if (colp[C[x][y][z]].x==-1 ||
	      B[colp[C[x][y][z]].x][colp[C[x][y][z]].y][colp[C[x][y][z]].z]<B[x][y][z]){
	    colp[C[x][y][z]].x=x;
	    colp[C[x][y][z]].y=y;
	    colp[C[x][y][z]].z=z;
	  }
	  col[C[x][y][z]].x+=B[x][y][z]*x*S;
	  col[C[x][y][z]].y+=B[x][y][z]*y*S;
	  col[C[x][y][z]].z+=B[x][y][z]*z*S;
	  col[C[x][y][z]].w+=B[x][y][z];	  
	}
      }

  for (int i=0;i<code;i++){
    col[i].x/=col[i].w;
    col[i].y/=col[i].w;
    col[i].z/=col[i].w;
    if (max1<col[i].w)
      max1=col[i].w;

    if (col[i].x>255) col[i].x=255;
    if (col[i].y>255) col[i].y=255;
    if (col[i].z>255) col[i].z=255;
  }
  

  printf("ok\n");
  //analizzo adiacenze col
  vector < vector< vector<int3> > > neicol;
  vector < vector<st> > nei;
  nei.resize(code);
  neicol.resize(code);
  for (int i=0;i<code;i++){
    nei[i].resize(code);
    neicol[i].resize(code);    
    for (int j=0;j<code;j++){
      neicol[i][j].clear();
      nei[i][j].m1=1e10;
      nei[i][j].M1=0;
      nei[i][j].m2=1e10;
      nei[i][j].M2=0;
      nei[i][j].best=0;
      nei[i][j].acc=0;
      nei[i][j].w=0;
    }
  }
  
  vector<long> listmax;
  vector<int3> listmaxpt;
  listmax.resize(code);
  listmaxpt.resize(code);

  vector < vector<int3> > listbest1; // punti partenza per ogni coppia faccia (interno)

  vector<int3> listbest; // punti partenza per ogni faccia (interno)
  vector<long> listbestv; // quanto riesco a fare scavallando li (uso peso della faccia adiacente)
  listbest.resize(code);
  listbest1.resize(code);
  listbestv.resize(code);
  for (int i=0;i<code;i++){
    listbestv[i]=0;
  listbest1[i].resize(code);
  }


  printf("ok\n");

    for (int x=0;x<NB;x++)
      for (int y=0;y<NB;y++)
	for (int z=0;z<NB;z++)
	  if (C[x][y][z]>0 
	      ){
	    int posx=x;
	    int posy=y+z*NB;	 

	    int c=0;
	    if (z<NB/2)
	      c=255;
	      //	    printf("%d %d %d, %d --> %d %d\n",x,y,z,C[x][y][z],posx,posy);
	    stat1(posx,posy,0,0)=col[C[x][y][z]].x*0.3+255*0.7;
	    stat1(posx,posy,0,1)=col[C[x][y][z]].y*0.3+255*0.7;
	    stat1(posx,posy,0,2)=col[C[x][y][z]].z*0.3+255*0.7;
	    
	    /*	    if (B[x][y][z]>0.1*col[C[x][y][z]].w){
	      stat1(posx,posy,0,0)=(int)((float)x*S*B[x][y][z]/max1*300);
	      stat1(posx,posy,0,1)=(int)((float)y*S*B[x][y][z]/max1*300);
	      stat1(posx,posy,0,2)=(int)((float)z*S*B[x][y][z]/max1*300);
	    }
	    */

	}

  printf("ok\n");

    /////////analisi monotona
      for (int x=0;x<NB;x++)
	for (int y=0;y<NB;y++)
	  for (int z=0;z<NB;z++)
	    if (B[x][y][z]>0){
	      int i=C[x][y][z];
	      // sono sul bordo con altra zona
	      int startok=0;
	      int c1=-1;
	      for (int dx=-1;dx<=1;dx++)
		for (int dy=-1;dy<=1;dy++)
		  for (int dz=-1;dz<=1;dz++)
		    //if (fabs(dx)+fabs(dy)+fabs(dz)<=1)
		      {
		      if (x+dx>=0 && y+dy>=0 && z+dz>=0 &&
			  x+dx<NB && y+dy<NB && z+dz<NB) 
			{
			  if (C[x+dx][y+dy][z+dz]!=i && C[x+dx][y+dy][z+dz]!=0){
			    startok=1;
			    c1=C[x+dx][y+dy][z+dz];
			    //printf("start %d %d %d, n %d %d %d, %d %d\n",x,y,z,x+dx,y+dy,z+dz,i,c1);
			  }
			}
		      else{
			// per ora non considero bordo
			startok=1; // bordo
			c1=0;
			//		      printf("starB %d %d %d (%d) n %d %d %d (%d)\n",
			//     x,y,z,i,x+dx,y+dy,z+dz,C[x+dx][y+dy][z+dz]);
		      }
		    }

	      // secondo punto non deve essere a contatto con l'altra faccia (non navigo lungo taglio, ma entro la mia regione)
	      //quale sarebbe il successivo
	      if(0)
	      if (startok) {
		long best=0;
		int nx=0,ny=0,nz=0;
		for (int dx1=-1;dx1<=1;dx1++)
		  for (int dy1=-1;dy1<=1;dy1++)
		    for (int dz1=-1;dz1<=1;dz1++)
		      if (x+dx1>=0 && y+dy1>=0 && z+dy1>=0 &&
			  x+dx1<NB && y+dy1<NB && z+dy1<NB)
			if (B[x+dx1][y+dy1][z+dz1]>best){
			  best=B[x+dx1][y+dy1][z+dz1];
			  nx=x+dx1;
			  ny=y+dy1;
			  nz=z+dz1;
			}
		if (nx!=x || ny!=y || nz!=z){ // prosegue
		  // controlla che nessun vicino abbia codice c1
		  for (int dx1=-1;dx1<=1;dx1++)
		    for (int dy1=-1;dy1<=1;dy1++)
		      for (int dz1=-1;dz1<=1;dz1++)
			if (fabs(dx1)+fabs(dy1)+fabs(dz1)<=1)
			  if (x+dx1>=0 && y+dy1>=0 && z+dy1>=0 &&
			      x+dx1<NB && y+dy1<NB && z+dy1<NB){
			    if (C[x+dx1][y+dy1][z+dz1]==c1){ // non mi sono allontanato!!
			      startok=0;
			      //printf("non stacca\n");
			    }
			  }
			  else{
			    if (c1==0) // se era bordo e sono ancora su bordo
			      startok=0;
			  }
		  
		}
	      }

	      if (startok){
		int ct1=0;
		int3 p,last;
		p.x=x;
		p.y=y;
		p.z=z;
		last.x=-1;
		float score=0;
		float ang=0;
		int mod=1;

		while (mod){
		  mod=0;
		  //next
		  long best=0;
		  int nx=0,ny=0,nz=0;
		  for (int dx=-1;dx<=1;dx++)
		    for (int dy=-1;dy<=1;dy++)
		      for (int dz=-1;dz<=1;dz++)
			if (p.x+dx>=0 && p.y+dy>=0 && p.z+dz>=0 &&
			    p.x+dx<NB && p.y+dy<NB && p.z+dz<NB)
			  //			  if (fabs(dx)+fabs(dy)+fabs(dz)<=2)
			    {
			      if (B[p.x+dx][p.y+dy][p.z+dz]>best){
				//		      printf("%d %d, %d %d %d\n",B[p.x+dx][p.y+dy][p.z+dz],best,dx,dy,dz);
				best=B[p.x+dx][p.y+dy][p.z+dz];
				nx=p.x+dx;
				ny=p.y+dy;
				nz=p.z+dz;
			      }		  
			    }
		  if (nx!=p.x || ny!=p.y || nz!=p.z){
		    mod=1;
		    score+=B[p.x][p.y][p.z];
		    if (last.x!=-1){
		      float v11,v12,v13;
		      float v21,v22,v23;
		      v11=p.x-last.x;
		      v12=p.y-last.y;
		      v13=p.z-last.z;
		      v21=nx-p.x;
		      v22=ny-p.y;
		      v23=nz-p.z;
		      float norm=pow(v11*v11+v12*v12+v13*v13,0.5f);
		      v11/=norm;v12/=norm;v13/=norm;
		      norm=pow(v21*v21+v22*v22+v23*v23,0.5f);
		      v21/=norm;v22/=norm;v23/=norm;
		      float t=v11*v21+v12*v22+v13*v23;
		      t=pow((t+1)/2,2.0f);
		      ang+=t;
		      //printf("%f %f %f - %f %f %f ",v11,v12,v13,v21,v22,v23);
		    }
		    //printf("%d %d %d -> %d %d %d: %f %f (%d)\n",p.x,p.y,p.z,nx,ny,nz,score,ang,last.x);	      
		    last.x=p.x;
		    last.y=p.y;
		    last.z=p.z;
		    p.x=nx;
		    p.y=ny;
		    p.z=nz;
		    ct1++;
		  }
		}
		//		printf("%d %d %d: n %d, score %f %f, ang %f --> index %f\n",x,y,z,ct1,score,score/ct1,ang/ct1,score*ang/ct1);
		E[x][y][z]=score/ct1/B[p.x][p.y][p.z]*ang; // normalizzo per intensita' max
		// va bene per selezione relativa tra coppie di zone
	      }
	    }

      vector<int3> pair1;
      vector<int3> pair2;
      vector<float> pairv;
      float maxpair=0;
      //superclassifica coppie
      for (int x=0;x<NB;x++)
	for (int y=0;y<NB;y++)
	  for (int z=0;z<NB;z++)
	    if (E[x][y][z]>0){
	      int max=1;
	      int s=1;

	      int3 bestn;
	      float bestnv=0;

	      for (int dx=-s;dx<=s;dx++)
		for (int dy=-s;dy<=s;dy++)
		  for (int dz=-s;dz<=s;dz++)
		    if (x+dx>=0 && y+dy>=0 && z+dz>=0 &&
			x+dx<NB && y+dy<NB && z+dz<NB)
		      if (C[x+dx][y+dy][z+dz]!=C[x][y][z] &&
			  C[x+dx][y+dy][z+dz]!=0 &&
			  E[x+dx][y+dy][z+dz]>0) {
			if (C[x][y][z]==5)
 			printf("Valuto %d %d\n",C[x][y][z],C[x+dx][y+dy][z+dz]);
			// devo risalire con entrambi e calcolare angoli
			int x1=x+dx;
			int y1=y+dy;
			int z1=z+dz;
		      
			long best=0;
			int nx=0,ny=0,nz=0;
			for (int dx1=-1;dx1<=1;dx1++)
			  for (int dy1=-1;dy1<=1;dy1++)
			    for (int dz1=-1;dz1<=1;dz1++)
			      if (x+dx1>=0 && y+dy1>=0 && z+dz1>=0 &&
				  x+dx1<NB && y+dy1<NB && z+dz1<NB){
				  if (B[x+dx1][y+dy1][z+dz1]>best){
				    //		      printf("%d %d, %d %d %d\n",B[p.x+dx][p.y+dy][p.z+dz],best,dx,dy,dz);
				    best=B[x+dx1][y+dy1][z+dz1];
				    nx=x+dx1;
				    ny=y+dy1;
				    nz=z+dz1;
				  }		  
				}
			best=0;
			int nx1=0,ny1=0,nz1=0;
			for (int dx1=-1;dx1<=1;dx1++)
			  for (int dy1=-1;dy1<=1;dy1++)
			    for (int dz1=-1;dz1<=1;dz1++)
			      if (x1+dx1>=0 && y1+dy1>=0 && z1+dz1>=0 &&
				  x1+dx1<NB && y1+dy1<NB && z1+dz1<NB){
				  if (B[x1+dx1][y1+dy1][z1+dz1]>best){
				    //		      printf("%d %d, %d %d %d\n",B[p.x+dx][p.y+dy][p.z+dz],best,dx,dy,dz);
				    best=B[x1+dx1][y1+dy1][z1+dz1];
				    nx1=x1+dx1;
				    ny1=y1+dy1;
				    nz1=z1+dz1;
				  }		  
				}
		      
			// angoli nx x x1 nx1
			float v11,v12,v13;
			float v21,v22,v23;
			float v31,v32,v33;
			v11=x1-nx1;
			v12=y1-ny1;
			v13=z1-nz1;
			v21=x-x1;
			v22=y-y1;
			v23=z-z1;
			v31=nx-x;
			v32=ny-y;
			v33=nz-z;
			float norm=pow(v11*v11+v12*v12+v13*v13,0.5f);
			v11/=norm;v12/=norm;v13/=norm;
			norm=pow(v21*v21+v22*v22+v23*v23,0.5f);
			v21/=norm;v22/=norm;v23/=norm;
			norm=pow(v31*v31+v32*v32+v33*v33,0.5f);
			v31/=norm;v32/=norm;v33/=norm;
			float t1=v11*v21+v12*v22+v13*v23;
			float t2=v31*v21+v32*v22+v33*v23;
			t1=pow((t1+1)/2,2.0f);
			t2=pow((t2+1)/2,2.0f);
			float corr=(t1*t2);
			//			printf("%d %d %d . %d %d %d -> %f * %f = %f\n",x,y,z,x+dx,y+dy,z+dz,E[x][y][z]+E[x+dx][y+dy][z+dz],corr,(E[x][y][z]+E[x+dx][y+dy][z+dz])*corr);
			// semplicemente sommo, per confronto relativo tra zone
			float newv=(E[x][y][z]+E[x+dx][y+dy][z+dz])*corr*(B[x][y][z]+B[x+dx][y+dy][z+dz]);
			if (newv>bestnv){
			  bestnv=newv;
			  bestn.x=x1;
			  bestn.y=y1;
			  bestn.z=z1;
			}
		      }

	      if (bestnv>0){
		if (C[x][y][z]==5)
		  printf("ok\n");
		int3 temp;
		temp.x=x;
		temp.y=y;
		temp.z=z;
		pair1.push_back(temp);
		temp.x=bestn.x;
		temp.y=bestn.y;
		temp.z=bestn.z;
		pair2.push_back(temp);
		pairv.push_back(bestnv);
		if (maxpair<bestnv)
		  maxpair=bestnv;
	      }

	    }

      printf("Quanti pair %d\n",pairv.size());
      if(0)
      for (int i=0;i<pairv.size();i++){
	printf("%d %d %d . %d %d %d -> %f\n",
	       pair1[i].x,pair1[i].y,pair1[i].z,
	       pair2[i].x,pair2[i].y,pair2[i].z,
	       pairv[i]);
      }

      // prendo il migliore tra due codici (perdo caso raro? di due doppie dita che si toccano)
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++){
	  nei[i][j].best=0;	  
	  nei[i][j].idx=-1;	  
	}

      for (int i=0;i<pairv.size();i++){
	//trovo coppia
	int c1,c2;
	c1=C[pair1[i].x][pair1[i].y][pair1[i].z];
	c2=C[pair2[i].x][pair2[i].y][pair2[i].z];
	if (c1>c2){
	  c2=C[pair1[i].x][pair1[i].y][pair1[i].z];
	  c1=C[pair2[i].x][pair2[i].y][pair2[i].z];
	  int3 temp;
	  temp=pair1[i];
	  pair1[i]=pair2[i];
	  pair2[i]=temp;
	}
	if (nei[c1][c2].best<pairv[i]){
	  //printf("canc %d %d, %d %f\n",c1,c2,nei[c1][c2].idx,pairv[nei[c1][c2].idx]);
	  nei[c1][c2].best=pairv[i];
	  if (nei[c1][c2].idx>=0) // cancello
	    pairv[nei[c1][c2].idx]=0;
	  nei[c1][c2].idx=i;
	}
      }



      //sort pair
      if(0)
      for (int j=0;j<pairv.size();j++){
	  for (int i=0;i<pairv.size()-1;i++){
	  if (pairv[i]<pairv[i+1]){
	    float t=pairv[i];pairv[i]=pairv[i+1];pairv[i+1]=t;
	    int3 p=pair1[i];pair1[i]=pair1[i+1];pair1[i+1]=p;
	    p=pair2[i];pair2[i]=pair2[i+1];pair2[i+1]=p;
	  }
	}
      }


      //rescoring: ora non mi interessa piu' il miglior percorso tra due zone
      //ma voglio valutare i percorsi tra loro
      /*
      maxpair=0;
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (nei[i][j].best>0){
	    int idx;
	    idx=nei[i][j].idx;
	    float v=
	      B[pair1[idx].x][pair1[idx].y][pair1[idx].z]+
	      B[pair2[idx].x][pair2[idx].y][pair2[idx].z];
	    nei[i][j].best=v;
	    pairv[idx]=v;
	    if (maxpair<v)
	      maxpair=v;
      	  }
      */

  FILE* gr1=fopen("graph1.dot","w+");
  fprintf(gr1,"digraph{\nnode [shape = circle, fixedsize=true, style = filled, fillcolor = palegreen];\n");
  for (int i=0;i<code;i++){
    if (i==0)
    fprintf(gr1,"%d [fillcolor =\".0 .0 .0\"];\n",i);
    else{
      float temp=pow(col[i].w/max1,0.5f);

      temp=10*temp+(1-temp)*0.5;
      fprintf(gr1,"%d [width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",i,
	      temp,
	      (int)col[i].x,(int)col[i].y,(int)col[i].z);
    }
  }

    ///seleziona massimi locali su coppie
    int ct2=0;
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (nei[i][j].best>0){
	    int pos=0;	  
	    int idx=nei[i][j].idx;
	    {
	      int x=pair1[idx].x;
	      int y=pair1[idx].y;
	      int z=pair1[idx].z;
	      int x1=pair2[idx].x;
	      int y1=pair2[idx].y;
	      int z1=pair2[idx].z;
	      float st1a=0;
	      float st1b=0;
	      float st2a=0;
	      float st2b=0;
	      float st1ac=0;
	      float st1bc=0;
	      float st2ac=0;
	      float st2bc=0;
	      float st1;
	      float st2;
	      float mint=1e10;
	      {
		vector<int3> temp;
		for (int k=0;k<2;k++) {
		  vector<int3> tempp;
		  int mod=1;
		  int3 p;
		  if (k==0){
		    p.x=x;
		    p.y=y;
		    p.z=z;
		  }
		  else{
		    p.x=x1;
		    p.y=y1;
		    p.z=z1;
		  }
		  while (mod){
		    mod=0;
		    int3 temp1;
		    temp1.x=p.x;
		    temp1.y=p.y;
		    temp1.z=p.z;
		    temp.push_back(temp1);	      
		    tempp.push_back(temp1);	      
		    //next
		    long best=0;
		    int nx,ny,nz;
		    for (int dx=-1;dx<=1;dx++)
		      for (int dy=-1;dy<=1;dy++)
			for (int dz=-1;dz<=1;dz++)
			  if (p.x+dx>=0 && p.y+dy>=0 && p.z+dz>=0 &&
			      p.x+dx<NB && p.y+dy<NB && p.z+dz<NB){
			    if (B[p.x+dx][p.y+dy][p.z+dz]>best){
			      best=B[p.x+dx][p.y+dy][p.z+dz];
			      nx=p.x+dx;
			      ny=p.y+dy;
			      nz=p.z+dz;
			    }		  
			  }
		    if (nx!=p.x || ny!=p.y || nz!=p.z){
		      mod=1;
		      p.x=nx;
		      p.y=ny;
		      p.z=nz;
		    }
		  }


		  // estraggo da tempp le stime gradiente
		  for (int i1=0;i1<tempp.size();i1++){
		    int la=tempp.size()-1;
		    float maxx=B[tempp[la].x][tempp[la].y][tempp[la].z];
		    float newv=B[tempp[i1].x][tempp[i1].y][tempp[i1].z]/maxx;
		    if (B[tempp[i1].x][tempp[i1].y][tempp[i1].z]<mint)
		      mint=B[tempp[i1].x][tempp[i1].y][tempp[i1].z];
		    if (k==0){
		      if (i1<tempp.size()/2)
			{st1a+=newv;st1ac++;}
		      else
			{st1b+=newv;st1bc++;}
		    }
		    else{
		    if (i1<tempp.size()/2)
		      {st2a+=newv;st2ac++;}
		    else
		      {st2b+=newv;st2bc++;}
		    }
		  }
		  //reverse se primo
		  if (k==0){
		    for (int i1=0;i1<temp.size()/2;i1++){
		      int3 t;
		      t.x=temp[i1].x;
		      t.y=temp[i1].y;
		      t.z=temp[i1].z;
		      temp[i1].x=temp[temp.size()-1-i1].x;
		      temp[i1].y=temp[temp.size()-1-i1].y;
		      temp[i1].z=temp[temp.size()-1-i1].z;
		      temp[temp.size()-1-i1].x=t.x;
		      temp[temp.size()-1-i1].y=t.y;
		      temp[temp.size()-1-i1].z=t.z;
		    }		  
		  }


		}
		st1a/=st1ac;
		st1b/=st1bc;
		st2a/=st2ac;
		st2b/=st2bc;
		st1=fabs(st1a-st1b);
		st2=fabs(st2a-st2b);
		float st=(st1+st2)/2;
		printf("%f %f, %f %f %f %f\n",st1,st2,st1a,st1b,st2a,st2b);
		printf("%f %f %f:\n",st,mint,(1-st)*mint);
		for (int i1=0;i1<temp.size();i1++){
		  printf("%d -> ", 
			 //temp[i1].x,temp[i1].y,temp[i1].z,
			 B[temp[i1].x][temp[i1].y][temp[i1].z]);
		  int posx=temp[i1].x;
		  int posy=temp[i1].y+temp[i1].z*NB;	 		  
		  stat1(posx,posy,0,0)=temp[i1].x*S;
		  stat1(posx,posy,0,1)=temp[i1].y*S;
		  stat1(posx,posy,0,2)=temp[i1].z*S;		  
		  stat1(i1,ct2,0,0)=temp[i1].x*S;
		  stat1(i1,ct2,0,1)=temp[i1].y*S;
		  stat1(i1,ct2,0,2)=temp[i1].z*S;
		}	
		//st=1-((1-st)*(1-dc));

		// copia in array
		for (int i1=0;i1<temp.size();i1++){
		  neicol[i][j].push_back(temp[i1]);
		  neicol[j][i].push_back(temp[temp.size()-1-i1]);
		}

		nei[i][j].col=0;
		for (int i1=0;i1<temp.size()-1;i1++){
		  float delta=pow(
		    pow((float)temp[i1].x-temp[i1+1].x,2.0f)+
		    pow((float)temp[i1].x-temp[i1+1].x,2.0f)+
		    pow((float)temp[i1].x-temp[i1+1].x,2.0f),0.5f);
		  nei[i][j].col+=delta;
		}	


		printf("\n%d-%d, %d: %d %d %d: %f (%f)\n",i,j,idx,x,y,z,nei[i][j].col,st);
		//riciclo best per mettere il valore
		nei[i][j].best=(1-st)*mint;

		ct2++;
		printf("\n");


	      }
	    }
	  }

      /*
      // calcolo colore
    ///seleziona massimi locali su coppie
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (nei[i][j].best>0){
	    int pos=0;	  
	    int idx=nei[i][j].idx;
	    int x=pair1[idx].x;
	    int y=pair1[idx].y;
	    int z=pair1[idx].z;
	    int x1=pair2[idx].x;
	    int y1=pair2[idx].y;
	    int z1=pair2[idx].z;
	    //distanza colore

	    float m1=(col[C[x][y][z]].x+col[C[x][y][z]].y+col[C[x][y][z]].z)/3/255;
	    float m2=(col[C[x][y][z]].x+col[C[x][y][z]].y+col[C[x][y][z]].z)/3/255;

	    float h1,s1,l1;
	    float h2,s2,l2;
	    rgb2hsl(col[C[x][y][z]].x/255.0,col[C[x][y][z]].y/255.0,col[C[x][y][z]].z/255.0,h1,s1,l1);
	    rgb2hsl(col[C[x1][y1][z1]].x/255.0,col[C[x1][y1][z1]].y/255.0,col[C[x1][y1][z1]].z/255.0,h2,s2,l2);

	    float deltah=fabs(h1-h2);
	    if (deltah>0.5)
	      deltah=fabs(1+h1-h2);

	    float mins=fmin(s1,s2);
	    float penalita=3*(1-mins)+1*(mins);
	    deltah*=penalita;
	    float dc=pow(pow(deltah,2.0f)+
			 pow(l1-l2,2.0f)+
			 pow(s1-s2,2.0f),0.2f);
	    //	    if (dc>1) dc=1;
			  
	    printf("%d %d -> %f pen %f (%f %f %f  %f %f %f)\n",i,j,dc,
		   penalita,
		   h1,s1,l1,
		   h2,s2,l2);
	    	    nei[i][j].col=dc;
	  }
*/

      float rangemin=1e10;
      float rangemax=0;
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (nei[i][j].best>0){
	    if (rangemin>nei[i][j].best)
	      rangemin=nei[i][j].best;
	    if (rangemax<nei[i][j].best)
	      rangemax=nei[i][j].best;
	  }      


      float rangeminc=1e10;
      float rangemaxc=0;
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++){
	    if (rangeminc>nei[i][j].col)
	      rangeminc=nei[i][j].col;
	    if (rangemaxc<nei[i][j].col)
	      rangemaxc=nei[i][j].col;
	  }      

      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  nei[i][j].arco=0;

      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (nei[i][j].best>0){
	    float stscaled=1-(nei[i][j].best-rangemin)/(rangemax-rangemin); // inverto, perche' ho proporzionale a minimo sella
	    float colscaled=(nei[i][j].col-rangeminc)/(rangemaxc-rangeminc)+0.001;
	    //stscaled=colscaled;
	    //colscaled=0;
	    //stscaled 0..1--> diventa velocita' con cui si attraversa arco
	    stscaled=stscaled*0.01+(1-stscaled)*1;
	    float time=colscaled/stscaled;
	    //if (time>1) time=1;
	    nei[i][j].arco=time;
	  }
      float rangemint=1e10;
      float rangemaxt=0;
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (nei[i][j].best>0){
	    if (rangemint>nei[i][j].arco)
	      rangemint=nei[i][j].arco;
	    if (rangemaxt<nei[i][j].arco)
	      rangemaxt=nei[i][j].arco;
	  }      

      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (nei[i][j].best>0){
	    nei[i][j].arco=(nei[i][j].arco-rangemint)/(rangemaxt-rangemint);
	    printf("%d -> %d: %f\n",i,j,nei[i][j].arco);
	  }
      /*      //prune
      if (0)
  for (int i=0;i<code;i++)
    for (int j=i+1;j<code;j++)
      if (nei[i][j].arco>0){ // per ogni arco
	  // se i >=3 e j>=3 archi e ij e' arco piu' debole --> tolgo

	float ref=nei[i][j].arco;
	float corr=col[i].w+col[j].w;
	ref*=corr;

	int cti=0;
	int ctj=0;

	for (int k=i+1;k<code;k++)
	  if (nei[i][k].arco>0){
	    float corr1=col[i].w+col[k].w;
	    if(corr1*nei[i][k].arco>=ref)
	    cti++;
	  }
	for (int k=1;k<i;k++)
	  if (nei[k][i].arco>0){
	    float corr1=col[k].w+col[i].w;
	    if (corr1*nei[k][i].arco>=ref)
	    cti++;
	  }
	for (int k=j+1;k<code;k++)
	  if (nei[j][k].arco>0){
	    float corr1=col[j].w+col[k].w;
	    if (corr1*nei[j][k].arco>=ref)
	    ctj++;
	  }
	for (int k=1;k<j;k++)
	  if (nei[k][j].arco>0){
	    float corr1=col[k].w+col[j].w;
	    if(corr1*nei[k][j].arco>=ref)
	    ctj++;
	}

	if (cti>=3 && ctj>=3){
	  printf("kkill %d %d\n",i,j);
	  nei[i][j].arco=0;
	}

	}


  if(0)
  for (int i=0;i<code;i++)
    for (int j=i+1;j<code;j++)
      for (int k=j+1;k<code;k++)
      if (nei[i][j].arco>0)
      if (nei[j][k].arco>0)
	if (nei[i][k].arco>0){

	  if (nei[i][j].arco < nei[i][k].arco && nei[i][j].arco < nei[j][k].arco){
	    nei[i][j].arco=0;
	    printf("kill %d %d (%d)\n",i,j,k);
	  }

	  if (nei[i][k].arco < nei[i][j].arco && nei[i][k].arco < nei[j][k].arco){
	    nei[i][k].arco=0;
	    printf("kill %d %d (%d)\n",i,k,j);
	  }

	  if (nei[j][k].arco < nei[i][j].arco && nei[j][k].arco < nei[i][k].arco){
	    nei[j][k].arco=0;
	    printf("kill %d %d (%d)\n",j,k,i);
	  }
	}
      */


    ///scrive grafo (solo per display nascondo archi brutti
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (nei[i][j].arco>0){
	    int pos=0;	  
	    int idx=nei[i][j].idx;
	    int x=pair1[idx].x;
	    int y=pair1[idx].y;
	    int z=pair1[idx].z;
	    int x1=pair2[idx].x;
	    int y1=pair2[idx].y;
	    int z1=pair2[idx].z;
	    //scrifo gr
	    float corr;
	    if (col[i].w<col[j].w)
	      corr=col[j].w/col[i].w;
	    else
	      corr=col[i].w/col[j].w;
	    corr=pow(corr,0.2f);
	    float arco=nei[i][j].arco*corr;
	    //if (arco<1)
	      {
	      printf("arco %d %d %f %f -> %f\n",i,j,nei[i][j].arco,corr,arco);
	      fprintf(gr1,"%d -> %d [dir=none penwidth=10 color = \"",C[x][y][z],C[x1][y1][z1]);
	      arco=pow(arco,0.5f);
	      float r=1,g=0,b=0;
	      if (arco<0.33) {r=arco*3; g=0;b=0;}
	      if (arco>=0.33 && arco<0.66) {r=1;g=(arco-0.33)*3;b=0;}
	      if (arco>0.66) {r=1;g=1;b=(arco-0.66)*3;}
	      if (arco>1) {r=0;g=1;}
	      fprintf(gr1,"#%02x%02x%02x",(int)(255*r),(int)(255*g),(int)(255*b));
	      fprintf(gr1,"\"];\n");
	    }
	  }
      fprintf(gr1,"}\n");
      fclose(gr1);


      // visita greedy
      vector<int> vis;
      int ctp=0;
      for (int starti=1;starti<code;starti++){
	int start=starti;
	printf("greedy ------------- %d\n",start);

	// seleziona migliori 2 da cui uscire
	int bi=-1,bi1=-1;
	float bs=1e10;
	for (int i=0;i<code;i++)
	  if (nei[start][i].arco+nei[i][start].arco>0)
	    for (int i1=0;i1<code;i1++)
	      if (i!=i1 && nei[start][i1].arco+nei[i1][start].arco>0){
		  float v1[3];
		  float v2[3];
		  float p1[3];
		  float p2[3];
		  float p3[3];
		  p1[0]=colp[i].x;		
		  p1[1]=colp[i].y;		
		  p1[2]=colp[i].z;		
		  p2[0]=colp[start].x;		
		  p2[1]=colp[start].y;		
		  p2[2]=colp[start].z;		
		  p3[0]=colp[i1].x;		
		  p3[1]=colp[i1].y;		
		  p3[2]=colp[i1].z;
		  v1[0]=p1[0]-p2[0];
		  v1[1]=p1[1]-p2[1];
		  v1[2]=p1[2]-p2[2];
		  v2[0]=p2[0]-p3[0];
		  v2[1]=p2[1]-p3[1];
		  v2[2]=p2[2]-p3[2];
		  normal(v1);
		  normal(v2);
		  float pen=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
		float corr;
		if (col[i].w<col[start].w)
		  corr=col[start].w/col[i].w;
		else
		  corr=col[i].w/col[start].w;
		corr=pow(corr,0.2f);
		float move1=(nei[start][i].arco+nei[i][start].arco)*corr*(2-pen);
		if (col[i1].w<col[start].w)
		  corr=col[start].w/col[i1].w;
		else
		  corr=col[i1].w/col[start].w;
		corr=pow(corr,0.2f);
		float move2=(nei[start][i1].arco+nei[i1][start].arco)*corr*(2-pen);
		float m=move1+move2;
		/*
		printf("a %f %f, m %f %f\n",nei[start][i].arco+nei[i][start].arco,
		       nei[start][i1].arco+nei[i1][start].arco,
		       move1,move2);
		*/
		if (m<bs && m<=1){ // non parto da salti grandi
		  bs=m;
		  bi=i;
		  bi1=i1;
		}
	      }
	printf("esco con %d %d, %f\n",bi,bi1,bs);

	vector<int> lista;

	for (int k=0;k<2;k++){
	vector<int> temp;
	  if (k==0)
	    start=bi;
	  else
	    start=bi1;
	  int last=starti;
	  vis.clear();
	  for (int i=0;i<code;i++)
	    vis.push_back(0);
	
	  while (start>0 && vis[start]==0){
	    temp.push_back(start);
	    printf("%d - ",start);
	    vis[start]=1;
	    //cerco neig
	    int n=-1;
	    float bn=1e10;
	    float penn=0;
	    for (int i=0;i<code;i++)
	      if (vis[i]==0){
		if (nei[start][i].arco+nei[i][start].arco>0){

		  float pen=1;
		  if (last>=0){ //angolo buono
		    float v1[3];
		    float v2[3];
		    float p1[3];
		    float p2[3];
		    float p3[3];
		    p1[0]=colp[i].x;		
		    p1[1]=colp[i].y;		
		    p1[2]=colp[i].z;		
		    p2[0]=colp[start].x;		
		    p2[1]=colp[start].y;		
		    p2[2]=colp[start].z;		
		    p3[0]=colp[last].x;		
		    p3[1]=colp[last].y;		
		    p3[2]=colp[last].z;
		    v1[0]=p1[0]-p2[0];
		    v1[1]=p1[1]-p2[1];
		    v1[2]=p1[2]-p2[2];
		    v2[0]=p2[0]-p3[0];
		    v2[1]=p2[1]-p3[1];
		    v2[2]=p2[2]-p3[2];
		    normal(v1);
		    normal(v2);
		    float ps=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
		    pen=ps;
		  }

		  float corr;
		  if (col[i].w<col[start].w)
		    corr=col[start].w/col[i].w;
		  else
		    corr=col[i].w/col[start].w;
		  corr=pow(corr,0.2f);
		  float move=(nei[start][i].arco+nei[i][start].arco)*corr*pow(2-pen,2.0f);
		  //printf("test %d arco %f move %f,  corr %f,  pen %f\n",i,nei[start][i].arco+nei[i][start].arco,move,corr,pen);
		  if (bn>move 
		      && pen>=0.5
		      && move<=0.3
		      ){
		    bn=move;
		    n=i;
		    penn=pen;
		  }	      
		}
	      }
	    printf("val %f pen %f\n",bn,penn);
	    last=start;
	    start=n;
	  }
	  printf("\n");
	if (k==0){
	  //inverte
	  for (int i=0;i<temp.size();i++)
	    lista.push_back(temp[temp.size()-1-i]);
	  lista.push_back(starti);
	}else{
	  for (int i=0;i<temp.size();i++)
	    lista.push_back(temp[i]);
	}

	}



	if (lista.size()>1){	  
	  lista1.resize(ctp+1);
	  for (int i=0;i<lista.size();i++)
	    lista1[ctp].push_back(lista[i]);
	  ctp++;
	}



      }


      /// analisi percorsi

      // spezzo loop
      if (0)
      for (int j=0;j<lista1.size();j++){
	int common=-1;
	for (int i=0;i<lista1[j].size();i++)
	  if (lista1[j][i]==lista1[j][lista1[j].size()-1-i])
	    common=i;
	  else
	    i=lista1[j].size();//break
	if (common>=0){
	  //cerca piu' dist
	  float bv=0;
	  int b=-1;
	  for (int i=common+1;i<lista1[j].size()-common-1;i++){
	    float v=
	      pow(col[lista1[j][i]].x-col[lista1[j][0]].x,2)+
	      pow(col[lista1[j][i]].y-col[lista1[j][0]].y,2)+
	      pow(col[lista1[j][i]].z-col[lista1[j][0]].z,2);
	    if (v>bv){
	      bv=v;
	      b=i;
	    }
	  }
	  
	  printf("loop %d: %d, split %d\n",j,common,b);

	  // aggiungo coda
	  int ls=lista1.size();
	  lista1.resize(ls+1);
	  for (int i=lista1[j].size()-1;i>b;i--){
	    lista1[ls].push_back(lista1[j][i]);
	  }
	  lista1[j].resize(b+1);

	} 
      }


      int mod=1;
      while(mod) {
	mod=0;
      for (int j=0;j<lista1.size();j++){
	printf("%d, %d : ",j,lista1[j].size());
	for (int i=0;i<lista1[j].size();i++)
	  printf("%d - ",lista1[j][i]);
	printf("\n");
      }

      //togli doppioni
      for (int j=0;j<lista1.size();j++){
	for (int j1=0;j1<lista1.size();j1++)
	  if (j!=j1 && lista1[j].size()>0 && 
	      lista1[j].size()==lista1[j1].size()){
	    int match=0;
	    for (int i1=0; i1<lista1[j].size();i1++){
	      if (lista1[j][i1]==lista1[j1][i1])
		match++;
	    }
	    if (match>0 && match==lista1[j].size()){ // rimuovo doppione
	      lista1[j1].clear();
	      printf("rimosso %d\n",j1);
	    }
	  }
      }
      // inverso
      for (int j=0;j<lista1.size();j++){
	for (int j1=0;j1<lista1.size();j1++)
	  if (j!=j1 && lista1[j].size()>0 && 
	      lista1[j].size()==lista1[j1].size()){
	    int match=0;
	    for (int i1=0; i1<lista1[j].size();i1++){
	      if (lista1[j][i1]==lista1[j1][lista1[j].size()-1-i1])
		match++;
	    }
	    if (match>0 && match==lista1[j].size()){ // rimuovo doppione
	      lista1[j1].clear();
	      printf("rimosso INV %d\n",j1);
	    }
	  }
      }

      for (int j=0;j<lista1.size();j++){
	for (int j1=0;j1<lista1.size();j1++)
	  if (j!=j1){
	    for (int i2=0;i2<lista1[j1].size();i2++){
	      int match=0;
	      for (int i1=0; i1<lista1[j].size() && i2+i1<lista1[j1].size();i1++){
		if (lista1[j][lista1[j].size()-1-i1]==lista1[j1][i2+i1])
		  match++;
	      }
	      if (match>0 && match==lista1[j].size()){ // rimuovo doppione
		lista1[j].clear();
		printf("sotto sequenza INV %d %d (%d)\n",j,j1,match);		  
	      }
	    }
	  }
      }


      for (int j=0;j<lista1.size();j++){
	for (int j1=0;j1<lista1.size();j1++)
	  if (j!=j1){
	    for (int i2=0;i2<lista1[j1].size();i2++){
	      int match=0;
	      for (int i1=0; i1<lista1[j].size() && i2+i1<lista1[j1].size();i1++){
		if (lista1[j][i1]==lista1[j1][i2+i1])
		  match++;
	      }
	      if (match>0 && match==lista1[j].size()){ // rimuovo doppione
		lista1[j].clear();
		printf("sotto sequenza %d %d (%d)\n",j,j1,match);		  
	      }
	    }
	  }
      }

      for (int j=0;j<lista1.size();j++){
	if (lista1[j].size()==0){
	  lista1.erase(lista1.begin()+j);
	  j--;
	}
      }
      }
      
      for (int j=0;j<lista1.size();j++){
	printf("%d: ",j);
	for (int i=0;i<lista1[j].size();i++)
	  printf("%d - ",lista1[j][i]);
	printf("\n");
      }

      node.resize(code);
      noder.resize(code);
      nodex.resize(code);

      adj.resize(code);
      for (int i=0;i<code;i++){
	adj[i].resize(code);
	for (int i1=0;i1<code;i1++)
	  adj[i][i1]=0;
      }
      ord.resize(code);
      for (int i=0;i<code;i++){
	ord[i].resize(code);
	for (int i1=0;i1<code;i1++)
	  ord[i][i1]=0;
      }

      for (int j=0;j<lista1.size();j++){

	float w=0;
	for (int k=0;k<lista1[j].size();k++)
	  w+=col[lista1[j][k]].w;

	for (int k=0;k<lista1[j].size()-1;k++){
	  int a=lista1[j][k];
	  int b=lista1[j][k+1];
	  adj[a][b]+=w;//sizen(a)+sizen(b);
	}
      }


      for (int i=0;i<code;i++)
	for (int i1=0;i1<code;i1++)
	  if (nei[i][i1].arco>0)
	  adj[i][i1]*=nei[i][i1].arco;






      int maxc=0;
      //tolgo cicli (arco piu' debole...)
      int restart1=1;
      while (restart1){
	restart1=0;
	for (int i=0;i<code;i++){
	  node[i]=-1;
	}
	for (int i=0;i<code;i++){
	  if (node[i]==-1){
	    ciclo.clear();
	    int ok=rec_v(i);
	    if (ok==-1){
	      float arc=10e10;
	      int idx=-1;
	      for (int k=0;k<ciclo.size()-1;k++){
		printf("%d ",ciclo[k]);
		if (idx==-1 || arc>adj[ciclo[k+1]][ciclo[k]]){
		  arc=adj[ciclo[k+1]][ciclo[k]];
		  idx=k;
		}		
		if (k>0 && ciclo[k+1]==ciclo[0])
		  k=ciclo.size();
	      }
	      printf("remove %d %d\n",ciclo[idx+1],ciclo[idx]);
	      adj[ciclo[idx+1]][ciclo[idx]]=0;

	      //spezzo lista che contiene

	      for (int j=0;j<lista1.size();j++){
		for (int k=0;k<lista1[j].size()-1;k++)
		  if (lista1[j][k]==ciclo[idx+1] &&
		      lista1[j][k+1]==ciclo[idx]){
		    //spezzo
		    lista1.resize(lista1.size()+1);
		    for (int k1=k+1;k1<lista1[j].size();k1++)
		      lista1[lista1.size()-1].push_back(lista1[j][k1]);
		    lista1[j].resize(k+1);
		    if (lista1[lista1.size()-1].size()<=1)
		      lista1.resize(lista1.size()-1); // ritolgo
		  }
	      }

	      restart1=1;
	      i=code;
	    }
	  }
	}
      }

      //tolgo corti
      for (int j=0;j<lista1.size();j++){
	if (lista1[j].size()<=1){
	  lista1.erase(lista1.begin()+j);
	  j--;
	}
      }
      
      printf("dopo spezzamento\n");
      for (int j=0;j<lista1.size();j++){

	float w=0;
	for (int k=0;k<lista1[j].size();k++)
	  w+=col[lista1[j][k]].w;

	printf("%d (%f): ",j,w);
	for (int i=0;i<lista1[j].size();i++)
	  printf("%d - ",lista1[j][i]);
	printf("\n");
      }


      //clean graph
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  for (int k=1;k<code;k++)
	    if (adj[i][j]>0)
	      if (adj[j][k]>0)
		if (adj[i][k]>0){
		  if (adj[i][j] <= adj[i][k] && adj[i][j] <= adj[j][k]){
		    adj[i][j]=0;
		    printf("kill %d %d (%d)\n",i,j,k);
		  }

		  if (adj[i][k] <= adj[i][j] && adj[i][k] <= adj[j][k]){
		    adj[i][k]=0;
		    printf("kill %d %d (%d)\n",i,k,j);
		  }

		  if (adj[j][k] <= adj[i][j] && adj[j][k] <= adj[i][k]){
		    adj[j][k]=0;
		    printf("kill %d %d (%d)\n",j,k,i);
		  }

		}











      FILE* gr2=fopen("graph2.dot","w+");
      fprintf(gr2,"digraph{\nsize=\"20.00,10.00\"\nnode [fontsize = 40, shape = circle, fixedsize=true, style = filled, fillcolor = palegreen];\n");
      for (int i=1;i<code;i++){

	  // se singleton salto
	  int uno=0;
	  for (int j=1;j<code;j++)
	    if (adj[i][j]>0 || adj[j][i]>0)
	      uno=1;
	  if (uno){
	    if (i==0)
	      fprintf(gr2,"%d [label =\"%d\", fillcolor =\".0 .0 .0\"];\n",i,node[i]);
	    else{
	      fprintf(gr2,"%d [label =\"%d\", width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",i,i,
		      sizen(i),
		      (int)col[i].x,(int)col[i].y,(int)col[i].z);
	    }
	  }
      }

      for (int i=0;i<code;i++){
	for (int i1=0;i1<code;i1++)
	  if (adj[i][i1]>0)
	    fprintf(gr2,"%d -> %d; \n",i,i1);
      }
      fprintf(gr2,"}\n");
      fclose(gr2);


      //disponi
      float maxv=0;
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
	if (maxv<B[x][y][z]) maxv=B[x][y][z];
      }
  printf("m %f\n",maxv);
      for (int x=0;x<NB;x++)
	for (int y=0;y<NB;y++)
	  for (int z=0;z<NB;z++)
      { //skippo l'ultimo a parte ultimo stretch
	    float h1,s1,l1;

	    //B[neicol[ni][nj][j].x][neicol[ni][nj][j].y][neicol[ni][nj][j].z]+
	    rgb2hsl(x*S/255.0,y*S/255.0,z*S/255.0,h1,s1,l1);
	    //	    rgb2hsl(col[i].x/255.0,col[i].y/255.0,col[i].z/255.0,h1,s1,l1);
	    //printf("%f %f %f -> %f %f %f\n",col[i].x/255.0,col[i].y/255.0,col[i].z/255.0,h1,s1,l1);

	  float quanto=(float)B[x][y][z]/maxv;
	  //if (z==0&&quanto>0) printf("%f\n",quanto);
	  quanto*=50;
	  if (quanto>1) quanto=1;

	    float px,py;
	    px=512+512*l1*sin(h1*2*3.1415);
	    py=512+512*l1*cos(h1*2*3.1415);

	    int lato=10;//(int)(sizen(i)*20);
	    for (int sidex=0;sidex<lato;sidex++)
	      for (int sidey=0;sidey<lato;sidey++){
		if ((int)px+sidex<1024 &&
		    (int)px+sidex>0 &&
		    (int)py+sidey<1024 &&
		    (int)py+sidey<1024>0){

		  stat1((int)px+sidex,(int)py+sidey,0,0)=(int)(x*S*quanto+(stat1((int)px+sidex,(int)py+sidey,0,0))*(1-quanto));
		  stat1((int)px+sidex,(int)py+sidey,0,1)=(int)(y*S*quanto+(stat1((int)px+sidex,(int)py+sidey,0,1))*(1-quanto));
		  stat1((int)px+sidex,(int)py+sidey,0,2)=(int)(z*S*quanto+(stat1((int)px+sidex,(int)py+sidey,0,2))*(1-quanto));

		}
	      }
	  }
      
      if (0)
      for (int i=0;i<code;i++){
	int ok=0;
	for (int i1=0;i1<code;i1++)
	  if (adj[i][i1]+adj[i1][i]>0)

	    for (int j=0;j<neicol[i][i1].size();j++){ //skippo l'ultimo a parte ultimo stretch
	    float h1,s1,l1;

	    //B[neicol[ni][nj][j].x][neicol[ni][nj][j].y][neicol[ni][nj][j].z]+
	    rgb2hsl(neicol[i][i1][j].x/255.0,neicol[i][i1][j].y/255.0,neicol[i][i1][j].z/255.0,h1,s1,l1);
	    //	    rgb2hsl(col[i].x/255.0,col[i].y/255.0,col[i].z/255.0,h1,s1,l1);
	    printf("%f %f %f -> %f %f %f\n",col[i].x/255.0,col[i].y/255.0,col[i].z/255.0,h1,s1,l1);

	    float px,py;
	    px=512+512*s1*sin(h1*2*3.1415);
	    py=512+512*s1*cos(h1*2*3.1415);

	    int lato=6;//(int)(sizen(i)*20);
	    for (int sidex=0;sidex<lato;sidex++)
	      for (int sidey=0;sidey<lato;sidey++){
		if ((int)px+sidex<1024 &&
		    (int)px+sidex>0 &&
		    (int)py+sidey<1024 &&
		    (int)py+sidey<1024>0){
		  stat1((int)px+sidex,(int)py+sidey,0,0)=col[i].x;
		  stat1((int)px+sidex,(int)py+sidey,0,1)=col[i].y;
		  stat1((int)px+sidex,(int)py+sidey,0,2)=col[i].z;
		  
		}
		  

	      }

	  }
      }







      int base=0;

      int loop=0;
	while (loop<20){
	loop++;

      float rangeminadj=1e20;
      float rangemaxadj=0;
      for (int i=1;i<code;i++)
	for (int j=1;j<code;j++)
	  if (adj[i][j]>0){
	    if (rangeminadj>adj[i][j])
	      rangeminadj=adj[i][j];
	    if (rangemaxadj<adj[i][j])
	      rangemaxadj=adj[i][j];
	  }      

      /////////
      vector<int> start;
      vector<int> end;
      vector<float> val;
      vector <int> prev;
      val.resize(code);
      prev.resize(code);
      for (int i=0;i<code;i++){
	int s=0;
	int u=0;
	for (int i1=0;i1<code;i1++){
	  if (adj[i1][i]>0 || adj[i][i1]>0)
	    s++;
	  if (adj[i][i1]>0)
	    u=1;
	}
	if (s>0)
	  start.push_back(i);
	int es=0;
	int eu=1;
	for (int i1=0;i1<code;i1++){
	  if (adj[i1][i]>0)
	    es=1;
	  if (adj[i][i1]>0)
	    eu=0;
	}

	if (es&&eu)
	  end.push_back(i);
      }

      printf("start\n");
      printf("s %d e %d\n",start.size(),end.size());
      vector<int> best;
      float bestv=0;
      for (int i=0;i<start.size();i++){
	//printf("%d\n",start[i]);

	fheap_t* H=fh_alloc(code);
	for (long i1=0;i1<code;i1++){
	  val[i1] =INFTY;
	  prev[i1]=-1;
	}

	val[start[i]]= 0;
	fh_insert(H, start[i], val[start[i]]);

	while (H->n>0){
	  long u=-1;  
	  u=fh_delete_min(H);
	  //printf("--- %d\n",u);
	  for (int i1=1;i1<code;i1++){
	    float v=adj[u][i1]+adj[i1][u];

	    //correzione per colore
	    float co=(col[i1].x-col[u].x+
		      col[i1].y-col[u].y+
		      col[i1].z-col[u].z)/3/255.0/2+0.5;
	    co=1*(1-co)+0.9*(co);
	    v*=co;

	    if (v>0){
	      float w=log(2-(v)/rangemaxadj)/log(2);
	      //printf("%d %d %f %f %f\n",u,i1,adj[u][i1],(adj[u][i1]*0.999)/rangemaxadj,w);
	      float newval=val[u]+w;
	      if (val[i1]==INFTY){
		fh_insert(H, i1, newval);
		val[i1] = newval;
		prev[i1]=u;		
	      }
	      else{
	      if (val[i1]>newval){
		fh_insert(H, i1, newval);
		val[i1] = newval;
		prev[i1]=u;		
		fh_decrease_key(H, i1, newval);											
	      }
	      }
	  }
	}
	}

	//risultati
	for (int i1=1;i1<code;i1++){
	  int u=i1;//end[i1];
	  int ct=1;
	  float len=0;
	  while (u!=-1 && prev[u]!=-1){
	    //printf("%d ",u);
	    len+=sizen(u);//col[u].w/max1;
	    u=prev[u];
	    ct++;
	  }
	  if (0&&val[i1]<INFTY)
	    printf(": end %d: val %f, ct %d (%f, %f),  %f %f\n",i1,val[i1],ct,len,val[i1]/len,pow(2,val[i1]),pow(2,val[i1])/ct);
	  float v=len;//pow(2,val[end[i1]])/ct
	  if (bestv< v){
	    bestv=v;
	    best.clear();
	    int u=i1;//end[i1];
	    int ct=0;
	    while (u!=-1){
	      best.insert(best.begin(),u);//push_back(u);
	      u=prev[u];
	    }
	  }
	}

      }

      if (bestv>0&&best.size()>1){
      printf("best (%f):\n",bestv);
      for (int i=0;i<best.size();i++)
	printf("%d ",best[i]);
      printf("\n");


      vector<int3> temp;
      vector<float> tempw;
      float totw=0;
      for (int b=0;b<best.size()-1;b++)
	totw+=sizen(best[b])+sizen(best[b+1]);

      for (int b=0;b<best.size()-1;b++){
	int ni=best[b];
	int nj=best[b+1];
	float w=0;
	float local=0;
	for (int j=0;j<neicol[ni][nj].size()-1;j++){ //skippo l'ultimo a parte ultimo stretch
	  float corr=(float)(j+0.5)/(neicol[ni][nj].size());
	  corr=pow(1+sin(corr*3.1415),2);
	  local+=corr*pow(
		     B[neicol[ni][nj][j].x][neicol[ni][nj][j].y][neicol[ni][nj][j].z]+
		     B[neicol[ni][nj][j+1].x][neicol[ni][nj][j+1].y][neicol[ni][nj][j+1].z]
		     ,0.5);
	}
	float gw=(sizen(ni)+sizen(nj))/totw;
	for (int j=0;j<neicol[ni][nj].size()-1;j++){ //skippo l'ultimo a parte ultimo stretch
	  float corr=(float)(j+0.5)/(neicol[ni][nj].size());
	  corr=pow(1+sin(corr*3.1415),2);
	  temp.push_back(neicol[ni][nj][j]);
	  float lw=corr*pow(
		       B[neicol[ni][nj][j].x][neicol[ni][nj][j].y][neicol[ni][nj][j].z]+
		       B[neicol[ni][nj][j+1].x][neicol[ni][nj][j+1].y][neicol[ni][nj][j+1].z],0.5)/local;

	  tempw.push_back(lw*gw);
	}
	if (b==best.size()-2)
	  temp.push_back(neicol[ni][nj][neicol[ni][nj].size()-1]);
      }




      int passo=100*sqrt(bestv);
      int side=(passo*stro1.height())/stro1.width();

      stro2.resize(passo,(passo*stro1.height())/stro1.width(),-100,-100,5);

      printf("%d x %d\n",passo,side);
	float ct=0;
	//printf("draw %d %d, %f, %d\n",temp.size(),tempw.size(),bestv,base);
	//printf("%d %d\n",side,passo);

	for (int px=base;px<base+side;px++)
	  if (px<1024)
	    for (int k=0;k<passo && k<1024;k++){
	      float alpha=stro2(k,px-base,0,0)/255.0;
	      float grigio=stro2(k,px-base,0,1);
	      ct=k*1.0/passo;

	      ct+=(alpha-0.5)*(grigio/255.0-0.5)*0.2;

	      //anticipo un po'
	      ct=ct*1.3;

	      if (ct<0) ct=0;
	      if (ct>0.999) ct=0.999;
	      //printf("k %d ct %f del %f\n",k,ct,(grigio/255.0-0.5)*0.5);
	      // cerco indice
	      int r1=0;
	      float r1f=tempw[r1++];
	      float r1fl=0;//last
	      //printf("%d %f %f\n",r1,r1f,r1fl);
	      while (r1f<ct) {r1fl=r1f;r1f+=tempw[r1++];}
	  
	      //printf("%f: %f %f %d (%d)\n",ct,r1fl,r1f,r1,temp.size());
	      //stat1(px,pos,0,0)=(1-se)*temp[k].x*S+se*temp[k+1].x*S;
	      float se=(ct-r1fl)/(r1f-r1fl);
	      float r,g,b;
	      r=((1-se)*temp[r1-1].x*S+se*temp[r1].x*S);
	      g=((1-se)*temp[r1-1].y*S+se*temp[r1].y*S);
	      b=((1-se)*temp[r1-1].z*S+se*temp[r1].z*S);
	      //printf("ok\n");

	      float r2,g1,b1;
	      r2=alpha*r+(1-alpha)*grigio; // grigio
	      g1=alpha*g+(1-alpha)*grigio; // grigio
	      b1=alpha*b+(1-alpha)*grigio; // grigio

	      stat1(k,px,0,0)=(int)r2;
	      stat1(k,px,0,1)=(int)g1;
	      stat1(k,px,0,2)=(int)b1;
	      //printf("%d %d %f\n",px,k,r);

	    }
	base+=side;

	int isolo=1;
	for (int i=0;i<best.size()-(1-isolo);i++){
	if (isolo)
	  for (int j=1;j<code;j++){ //isolo il nodo
	  adj[best[i]][j]=0;
	  adj[j][best[i]]=0;
	}
	else{
	adj[best[i]][best[i+1]]=0;
	adj[best[i+1]][best[i]]=0;
	}
      }
      }
	}

	/*
      int a=fork();
      if (a==0){
	printf("eseguo\n");
	execlp("./go1", "go1", NULL);
	exit(0);
      }else{
	printf("aspetto\n");
	int retv;
	wait(&retv);
      }


      vector <float2> p;
      p.resize(code);
	
      FILE* fi=fopen("n.txt","r");
      while (!feof(fi)){
	int i;
	float2 t;
	fscanf(fi,"%d %f %f\n",&i,&t.x,&t.y);
	//printf("%d %f %f\n",i,t.x,t.y);
	p[i]=t;
      }
      fclose(fi);
*/

	/*
      if (0){
      FILE* fo=fopen("output.ps","w+");

      fprintf(fo,"%%!PS-Adobe-3.0 EPSF-3.0\n%%%%Creator: me\n%%%%Title: test\n%%%%CreationDate: today\n%%%%DocumentData: Clean7Bit\n%%%%Origin: 0 0\n%%%%BoundingBox: 0 0 60 60\n%%%%LanguageLevel: 2\n%%%%Pages: 1\n%%%%Page: 1 1\n");
      for (int ni=1;ni<code;ni++)
	  for (int nj=ni+1;nj<code;nj++)
	    if (adj[ni][nj]+ adj[nj][ni]>0){
	int ctj=0;
	vector<int3> temp;
	vector<float> tempw;
	float w=0;
	float local=0;
	for (int j=0;j<neicol[ni][nj].size()-1;j++) //skippo l'ultimo a parte ultimo stretch
	  local+=pow(
		     B[neicol[ni][nj][j].x][neicol[ni][nj][j].y][neicol[ni][nj][j].z]+
		     B[neicol[ni][nj][j+1].x][neicol[ni][nj][j+1].y][neicol[ni][nj][j+1].z]
		     ,0.5);
	for (int j=0;j<neicol[ni][nj].size()-1;j++){ //skippo l'ultimo a parte ultimo stretch
	  temp.push_back(neicol[ni][nj][j]);
	  float lw=pow(
		       B[neicol[ni][nj][j].x][neicol[ni][nj][j].y][neicol[ni][nj][j].z]+
		       B[neicol[ni][nj][j+1].x][neicol[ni][nj][j+1].y][neicol[ni][nj][j+1].z],0.5)/local;
	  tempw.push_back(lw);
	}
	temp.push_back(neicol[ni][nj][neicol[ni][nj].size()-1]);

	float passo=100;
	int side=10;
	float ct=0;

	for (int k=0;k<passo;k++){
	  ct=k*1.0/passo;
	  // cerco indice
	  int r1=0;
	  float r1f=tempw[r1++];
	  float r1fl=0;//last
	  while (r1f<ct) {r1fl=r1f;r1f+=tempw[r1++];}
	  printf("%f: %f %f %d (%d)\n",ct,r1fl,r1f,r1,temp.size());
	  //stat1(px,pos,0,0)=(1-se)*temp[k].x*S+se*temp[k+1].x*S;
	  float se=(ct-r1fl)/(r1f-r1fl);
	  float r,g,b;
	  r=((1-se)*temp[r1-1].x*S+se*temp[r1].x*S)/255.0;
	  g=((1-se)*temp[r1-1].y*S+se*temp[r1].y*S)/255.0;
	  b=((1-se)*temp[r1-1].z*S+se*temp[r1].z*S)/255.0;
	  float x1,x2,y1,y2;
	  float ct1=ct+1.0/passo;
	  float d1=(1-ct)*dims+(ct)*dime;
	  float d2=(1-ct1)*dims+(ct1)*dime;
	  x1=(1-ct)*p1.x+(ct)*p2.x;
	  y1=(1-ct)*p1.y+(ct)*p2.y;
	  x2=(1-ct1)*p1.x+(ct1)*p2.x;
	  y2=(1-ct1)*p1.y+(ct1)*p2.y;
	  fprintf(fo,"newpath\n0 setlinewidth\n%3.3f %3.3f moveto\n%3.3f %3.3f lineto\n%3.3f %3.3f lineto\n%3.3f %3.3f lineto\nclosepath\ngsave\n%1.3f %1.3f %1.3f setrgbcolor\nfill\ngrestore\n%1.3f %1.3f %1.3f setrgbcolor\n1 setlinewidth\nstroke\n",
		  x1-no.x*d1,y1-no.y*d1,
		  x2-no.x*d2,y2-no.y*d2,
		  x2+no.x*d2,y2+no.y*d2,
		  x1+no.x*d1,y1+no.y*d1,
		  r,g,b,r,g,b);

	  printf("%f %f %f\n",
		 (1-se)*temp[r1-1].x*S+se*temp[r1].x*S,
		 (1-se)*temp[r1-1].y*S+se*temp[r1].y*S,
		 (1-se)*temp[r1-1].z*S+se*temp[r1].z*S);

	}		
	fprintf(fo,"%%%%EOF\n");

*/

	/*

newpath
100 100 moveto
0 100 rlineto
100 0 rlineto
0 -100 rlineto
-100 0 rlineto
closepath
gsave
0.5 1 0.5 setrgbcolor
fill
grestore
1 0 0 setrgbcolor
4 setlinewidth
stroke
showpage
*/



	/*
      }
      
      fclose(fo);
	    }
*/

      lista2.resize(lista1.size());
      for (int i=0;i<lista1.size();i++){
	for (int j=0;j<lista1[i].size();j++)
	  lista2[i].push_back(lista1[i][j]);
      }

      //display
      if(0)
      for (int j1=0;j1<lista2.size();j1++){
	int ctj=0;
	vector<int3> temp;
	vector<float> tempw;
	float w=0;
	for (int i=0;i<lista2[j1].size()-1;i++){
	  int ni=lista2[j1][i],nj=lista2[j1][i+1];
	  //printf("%d -> %d: %d\n",ni,nj,neicol[ni][nj].size());
	  int last=1;
	  if (i==lista2[j1].size()-2)
	    last=0;
	  float local=0;
	  for (int j=0;j<neicol[ni][nj].size()-last;j++) //skippo l'ultimo a parte ultimo stretch
	    local+=B[neicol[ni][nj][j].x][neicol[ni][nj][j].y][neicol[ni][nj][j].z];
	  float warco=sizen(ni)+sizen(nj);

	  for (int j=0;j<neicol[ni][nj].size()-last;j++){ //skippo l'ultimo a parte ultimo stretch
	    temp.push_back(neicol[ni][nj][j]);
	    float lw=B[neicol[ni][nj][j].x][neicol[ni][nj][j].y][neicol[ni][nj][j].z]/local;
	    tempw.push_back(lw*warco);
	  }
	  w+=warco;
	}
	for (int k=0;k<temp.size();k++)
	  printf("%1.4f ",(float)tempw[k]/w*100);
	printf("\n");
	  
	float ct=0*sizey;
	int side=10;
	for (int k=0;k<temp.size()-1;k++){
	  float newct=ct+(float)tempw[k]/w*sizey;
	  for (int pos=(int)ceil(ct);pos<=(int)floor(newct);pos++){
	    for (int px=side*j1;px<side*(j1+1);px++){
	      float se=(pos-ct)/(floor(newct)+1-ceil(ct));
	      //printf("%d %f %f: %f: %d %d -> %f\n",pos,ct,newct,se,temp[k].x,temp[k+1].x,(1-se)*temp[k].x*S+se*temp[k+1].x*S);
	      stat1(px,pos,0,0)=(1-se)*temp[k].x*S+se*temp[k+1].x*S;
	      stat1(px,pos,0,1)=(1-se)*temp[k].y*S+se*temp[k+1].y*S;
	      stat1(px,pos,0,2)=(1-se)*temp[k].z*S+se*temp[k+1].z*S;
	    }
	  }
	  ct=newct;
	}		
      }










    ///seleziona massimi locali su bordo 
    if (0)
    for (int i=1;i<code;i++){
      printf("--------- %d\n",i);
      int pos=0;	  
      for (int x=0;x<NB;x++)
	for (int y=0;y<NB;y++)
	  for (int z=0;z<NB;z++)
	    if (C[x][y][z]==i && E[x][y][z]>0){
	      int bordo=0;
	      int max=1;
	      int s=1;
	      for (int dx=-s;dx<=s;dx++)
		for (int dy=-s;dy<=s;dy++)
		  for (int dz=-s;dz<=s;dz++)
		    if (x+dx>=0 && y+dy>=0 && z+dz>=0 &&
			x+dx<NB && y+dy<NB && z+dz<NB) {
		      if (C[x+dx][y+dy][z+dz]==i && 
			  E[x+dx][y+dy][z+dz]>E[x][y][z])
			max=0;
			}
		    else
		      bordo=1;
	      if (max && bordo){
		vector<int3> temp;
		printf("%d: %d %d %d: %f\n",i,x,y,z,E[x][y][z]);
		int mod=1;
		int3 p;
		p.x=x;
		p.y=y;
		p.z=z;
		while (mod){
		  mod=0;
		  int3 temp1;
		  temp1.x=p.x;
		  temp1.y=p.y;
		  temp1.z=p.z;
		  temp.push_back(temp1);	      
		  //next
		  long best=0;
		  int nx,ny,nz;
		  for (int dx=-1;dx<=1;dx++)
		    for (int dy=-1;dy<=1;dy++)
		      for (int dz=-1;dz<=1;dz++)
			if (p.x+dx>=0 && p.y+dy>=0 && p.z+dz>=0 &&
			    p.x+dx<NB && p.y+dy<NB && p.z+dz<NB){
			  if (B[p.x+dx][p.y+dy][p.z+dz]>best){
			    best=B[p.x+dx][p.y+dy][p.z+dz];
			    nx=p.x+dx;
			    ny=p.y+dy;
			    nz=p.z+dz;
			  }		  
			}
		  if (nx!=p.x || ny!=p.y || nz!=p.z){
		    mod=1;
		    p.x=nx;
		    p.y=ny;
		    p.z=nz;
		  }
		}
		for (int i1=0;i1<temp.size();i1++){
		  printf("%d %d %d (%d) -> ", temp[i1].x,temp[i1].y,temp[i1].z,B[temp[i1].x][temp[i1].y][temp[i1].z]);
		  int posx=temp[i1].x;
		  int posy=temp[i1].y+temp[i1].z*NB;	 		  
		  stat1(posx,posy,0,0)=temp[i1].x*S;
		  stat1(posx,posy,0,1)=temp[i1].y*S;
		  stat1(posx,posy,0,2)=temp[i1].z*S;		  
		  stat1(i1,ct2,0,0)=temp[i1].x*S;
		  stat1(i1,ct2,0,1)=temp[i1].y*S;
		  stat1(i1,ct2,0,2)=temp[i1].z*S;		  
		  
		}
		ct2++;
		printf("\n");

	      }
	    }
    }


    stat1.save("stats1.jpg");
    
    int ct=0;
    for (int i=1;i<code;i++)
    for (int j=1;j<code;j++)
      if(nei[i][j].M1>0 && nei[i][j].M2>0)
	ct++;

    CImg<int> stat2(100,ct,1,3,255);

    ct=0;
    // disegna best iterativo
    for (int i=1;i<code;i++)
      for (int j=1;j<code;j++)
	if(nei[i][j].M1>0 && nei[i][j].M2>0) {
	  vector<int3> temp;
	  int pos=0;
	  for (int k=0;k<2;k++){

	    int3 p;
	    if (k==0){
	      p.x=nei[i][j].x1;
	      p.y=nei[i][j].y1;
	      p.z=nei[i][j].z1;
	    }
	    else{
	      p.x=nei[i][j].x2;
	      p.y=nei[i][j].y2;
	      p.z=nei[i][j].z2;
	    }
	
	    int mod=1;
	    while (mod){
	      mod=0;
	      int3 temp1;
	      temp1.x=p.x*S;
	      temp1.y=p.y*S;
	      temp1.z=p.z*S;
	      temp.push_back(temp1);	      
	      //next
	      long best=0;
	      int nx,ny,nz;
	      for (int dx=-1;dx<=1;dx++)
		for (int dy=-1;dy<=1;dy++)
		  for (int dz=-1;dz<=1;dz++)
		    if (p.x+dx>=0 && p.y+dy>=0 && p.z+dz>=0 &&
			p.x+dx<NB && p.y+dy<NB && p.z+dz<NB){
		      if (B[p.x+dx][p.y+dy][p.z+dz]>best){
			best=B[p.x+dx][p.y+dy][p.z+dz];
			nx=p.x+dx;
			ny=p.y+dy;
			nz=p.z+dz;
		      }		  
		    }
	      if (nx!=p.x || ny!=p.y || nz!=p.z){
		mod=1;
		p.x=nx;
		p.y=ny;
		p.z=nz;
	      }
	    }
	    if (k==0) // dopo primo, inverto

	      for (int i1=0;i1<temp.size()/2;i1++){
		int3 t1;
		t1=temp[i1];
		temp[i1]=temp[temp.size()-i1-1];
		temp[temp.size()-i1-1]=t1;
	      }
	  }
	  for (int i1=0;i1<temp.size();i1++){
	    stat2(i1,ct,0,0)=temp[i1].x;
	    stat2(i1,ct,0,1)=temp[i1].y;
	    stat2(i1,ct,0,2)=temp[i1].z;
	  }

	  ct++;
	}

    stat2.save("stats2.jpg");






}

void linea(int p1x,int p1y,int p2x, int p2y){
  printf("linea %d %d %d %d\n",p1x,p1y,p2x,p2y);
  //  for (int x=0;x<sx;x+=1){
  //  for (int y=0;y<sy;y+=1){
  int steps=abs(p1x-p2x);
  if (steps<abs(p1y-p2y))
    steps=abs(p1y-p2y);

  for (int s=0;s<steps;s++) {
    float s1=(float)s/steps;
    float s2=(float)(s+1)/steps;
    int x=(int)round(s1*p2x+(1-s1)*p1x);
    int y=(int)round(s1*p2y+(1-s1)*p1y);
    int x1=(int)round(s2*p2x+(1-s2)*p1x);
    int y1=(int)round(s2*p2y+(1-s2)*p1y);
    //    printf("%d %d %f, %d %d - %d %d\n",s,steps,(s1*p1y+(1-s1)*p2y),x,y,x1,y1);
    unsigned char r,g,b;
    unsigned char r1,g1,b1;
    r=src(x,y,0,0);
    g=src(x,y,0,1);
    b=src(x,y,0,2);
      
      r1=src(x1,y1,0,0);
      g1=src(x1,y1,0,1);
      b1=src(x1,y1,0,2);

      /// linea di colore
      int stepsc=abs(r-r1);
      if (stepsc<abs(g-g1))
	stepsc=abs(g-g1);
      if (stepsc<abs(b-b1))
	stepsc=abs(b-b1);

      if (stepsc==0){ // self loop!
    int r2=r;
    int g2=g;
    int b2=b;
      A[r2/S][g2/S][b2/S]++;
    //    printf("%d %d %d -> %d %d %d\n",r2,g2,b2,r3,g3,b3);
	int found=-1;
	for (int ct=0;ct<alist[r2/S][g2/S][b2/S].size();ct++)
	  if (alist[r2/S][g2/S][b2/S][ct].x==r2/S &&
	      alist[r2/S][g2/S][b2/S][ct].y==g2/S &&
	      alist[r2/S][g2/S][b2/S][ct].z==b2/S){
	    found=ct;
	    ct=alist[r2/S][g2/S][b2/S].size();
	  }
	if (found>=0){
	  alist[r2/S][g2/S][b2/S][found].w+=1;//+abs(r/S-r1/S)+abs(g/S-g1/S)+abs(b/S-b1/S);
	}
	else{//aggiungi
	  edge t;
	  t.x=r2/S;
	  t.y=g2/S;
	  t.z=b2/S;
	  t.w=1.0;
	  alist[r2/S][g2/S][b2/S].push_back(t);
	}
      }
      
      //printf("\n");
  for (int sc=0;sc<stepsc;sc++) {
    float sc1=(float)sc/stepsc;
    float sc2=(float)(sc+1)/stepsc;
    int r2=(int)round(sc1*r1+(1-sc1)*r);
    int g2=(int)round(sc1*g1+(1-sc1)*g);
    int b2=(int)round(sc1*b1+(1-sc1)*b);
    int r3=(int)round(sc2*r1+(1-sc2)*r);
    int g3=(int)round(sc2*g1+(1-sc2)*g);
    int b3=(int)round(sc2*b1+(1-sc2)*b);
    
    A[r3/S][g3/S][b3/S]++;
    if (sc==0)
      A[r2/S][g2/S][b2/S]++;
    
    //    printf("%d %d %d -> %d %d %d\n",r2,g2,b2,r3,g3,b3);
	int found=-1;
	for (int ct=0;ct<alist[r2/S][g2/S][b2/S].size();ct++)
	  if (alist[r2/S][g2/S][b2/S][ct].x==r3/S &&
	      alist[r2/S][g2/S][b2/S][ct].y==g3/S &&
	      alist[r2/S][g2/S][b2/S][ct].z==b3/S){
	    found=ct;
	    ct=alist[r2/S][g2/S][b2/S].size();
	  }
	if (found>=0){
	  alist[r2/S][g2/S][b2/S][found].w+=1.0/(stepsc);//+abs(r/S-r1/S)+abs(g/S-g1/S)+abs(b/S-b1/S);
	}
	else{//aggiungi
	  edge t;
	  t.x=r3/S;
	  t.y=g3/S;
	  t.z=b3/S;
	  t.w=1.0/(stepsc);
	  alist[r2/S][g2/S][b2/S].push_back(t);
	}

	found=-1;
	for (int ct=0;ct<alist[r3/S][g3/S][b3/S].size();ct++)
	  if (alist[r3/S][g3/S][b3/S][ct].x==r2/S &&
	      alist[r3/S][g3/S][b3/S][ct].y==g2/S &&
	      alist[r3/S][g3/S][b3/S][ct].z==b2/S){
	    found=ct;
	    ct=alist[r3/S][g3/S][b3/S].size();
	  }
	if (found>=0){
	  alist[r3/S][g3/S][b3/S][found].w+=1.0/(stepsc);//+abs(r/S-r1/S)+abs(g/S-g1/S)+abs(b/S-b1/S);
	}
	else{//aggiungi
	  edge t;
	  t.x=r2/S;
	  t.y=g2/S;
	  t.z=b2/S;
	  t.w=1.0/(stepsc);
	  alist[r3/S][g3/S][b3/S].push_back(t);
	}
    }
  }
}



// r,g,b values are from 0 to 255
// h = [0,1], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)

void RGBtoHSV( float r, float g, float b, float *h, float *s, float *v )
{
	float min, max, delta;
	r/=255;
	g/=255;
	b/=255;
	min = r;
	if (min>g) min=g;
	if (min>b) min=b;
	max = r;
	if (max<g) max=g;
	if (max<b) max=b;
	
	*v = max;				// v

	delta = max - min;

	if( max != 0 )
		*s = delta / max;		// s
	else {
		// r = g = b = 0		// s = 0, v is undefined
		*s = 0;
		*h = 0;
		return;
	}

	if (delta==0) *h=0;
	else{
	if( r == max )
		*h = ( g - b ) / delta;		// between yellow & magenta
	else if( g == max )
		*h = 2 + ( b - r ) / delta;	// between cyan & yellow
	else
		*h = 4 + ( r - g ) / delta;	// between magenta & cyan

	*h *= 60;				// degrees
	if( *h < 0 )
		*h += 360;
	*h/=360;
	}
}

void palettei(int i, int&r, int&g, int& b)
{
 int a=(i%15)+1;
 r=rand()%128+128;;//255*(1*(a&4)/4);
 g=rand()%128+128;;//255*(1*((a&2)/2+(a&8)/8)/2);
 b=128*(i%2)+128;;//255*(1*(a&1)/1);
}

int main(int argc, char *argv[]) {

  CImg<int> stat11(1024,1024,1,3,255);

  CImg<unsigned char> stro("stroke1.jpg");
  stro1=CImg<unsigned char>(stro);
  stro2=CImg<unsigned char>(stro);
  {
    int width = stro.width();
    int height = stro.height();
    unsigned char r,g,b;
	unsigned char gr1 = 0;
	unsigned char gr2 = 0;
	int max=0;
	/* Convert RGB image to grayscale image */
	for(int i=0;i<width;i++){
	  for(int j=0;j<height;j++){
 
	    //Return a pointer to a located pixel value. 
	    r = stro(i,j,0,0); // First channel RED
	    g = stro(i,j,0,1); // Second channel GREEN
	    b = stro(i,j,0,2); // Third channel BLUE

	    gr1 = round(0.333*((double)r) + 0.333*((double)g) + 0.333*((double)b));
 
	    //PAL and NTSC
	    //Y = 0.299*R + 0.587*G + 0.114*B 
	    int t=abs(gr1-r);
	    stro1(i,j,0,0) = t;
	    stro1(i,j,0,1) = gr1;
	    stro1(i,j,0,2) = gr1;
	    if (max<t) max=t;
	  } 
	}
	max=(int)(max*0.6);
	printf("max %d\n",max);
	for(int i=0;i<width;i++){
	  for(int j=0;j<height;j++){
	    int t=(int)((float)stro1(i,j,0,0)/max*255);
	    if (t>255) t=255;
	    stro1(i,j,0,0) = t;
	    stro1(i,j,0,1) = stro1(i,j,0,1);
	    stro1(i,j,0,2) = stro1(i,j,0,2);
	  }
	}

  }

  stro2=stro1;
  stro1.save("stroke2.jpg");

  src=CImg<unsigned char>(argv[1]);

  int sx = src.width();
  int sy = src.height();
  src.blur(1);
  printf("size %dx%d\n",sx,sy);
  if (0)
  for (int x=0;x<sx;x++){
    for (int y=0;y<sy;y++){
      printf("%3d %3d %3d.",
	     src(x,y,0,0),
	     src(x,y,0,1),
	     src(x,y,0,2)
	     );
    }
    printf("\n");
  }
  

  for (int x=0;x<NB;x++)
  for (int y=0;y<NB;y++)
    for (int z=0;z<NB;z++){
    A[x][y][z]=0;
    B[x][y][z]=0;
    C[x][y][z]=0;
    pos[x][y][z][0]=0;
    pos[x][y][z][1]=0;
    alist[x][y][z].clear();    
    gnode[x][y][z]=-1;
    }

  srand(0);

  //linea(0,0,0,sy-1);

  int th=divisioni;//sx;
  if (th>sx) th=sx;
  int step=1;//+sx/50;
    for (int q=0;q<th;q+=step){
    int p1x,p1y,p2x,p2y;
    p1x=rand()%sx;
    p1y=rand()%sy;
    p2x=rand()%sx;
    p2y=rand()%sy;
    p1x=sx*q/th;
    p1y=0;
    p2x=sx*q/th;
    p2y=sy-1;
    linea(p1x,p1y,p2x,p2y);
  }
    th=divisioni;//sy;
  if (th>sy) th=sy;
  step=1;
  for (int q=0;q<th;q+=step){
    int p1x,p1y,p2x,p2y;
    p1x=0;
    p1y=sy*q/th;
    p2x=sx-1;
    p2y=sy*q/th;
    linea(p1x,p1y,p2x,p2y);
  }

    
  //normalizza per distanza rgb arco
  if(0)
  for (int x=0;x<NB;x++)
  for (int y=0;y<NB;y++)
    for (int z=0;z<NB;z++){
      if (alist[x][y][z].size()>0){
	for (int ct=0;ct<alist[x][y][z].size();ct++){
	  float norm=
	    pow(alist[x][y][z][ct].x-x,2.0f)+
	    pow(alist[x][y][z][ct].y-y,2.0f)+
	    pow(alist[x][y][z][ct].z-z,2.0f);
	  //alist[x][y][z][ct].w=1;
	  if (norm>0)
	    alist[x][y][z][ct].w/=pow(norm,1);
	}
      }
    }



  if (0)
  for (int x=0;x<NB;x++)
  for (int y=0;y<NB;y++)
    for (int z=0;z<NB;z++){
      if (alist[x][y][z].size()>0){
	printf("%d %d %d: ",x,y,z);
	for (int ct=0;ct<alist[x][y][z].size();ct++)
	  printf("%d %d %d %3.0f, ",alist[x][y][z][ct].x,alist[x][y][z][ct].y,alist[x][y][z][ct].z,alist[x][y][z][ct].w);
	printf("\n");
      }
      }

  /////////prova sp
  int nnode=0;
  //int gnode[NB][NB][NB];
  //vector<char3> gchar3;

  if (0)
    for (int x=0;x<NB;x++)
  for (int y=0;y<NB;y++)
    for (int z=0;z<NB;z++){
      if (alist[x][y][z].size()>0){
	printf("%d %d %d: %d\n",x,y,z,alist[x][y][z].size());
	for (int ct=0;ct<alist[x][y][z].size();ct++){
	  printf("  %d %d %d %f\n",alist[x][y][z][ct].x,alist[x][y][z][ct].y,alist[x][y][z][ct].z,alist[x][y][z][ct].w);
	}
      }
    }

      
       for (int x=0;x<NB;x++)
  for (int y=0;y<NB;y++)
    for (int z=0;z<NB;z++){
      if (alist[x][y][z].size()>0){
	gnode[x][y][z]=nnode++;
	char3 t;
	t.x=x;
	t.y=y;
	t.z=z;
	gchar3.push_back(t);
      }
      else
	gnode[x][y][z]=-1;
    }

  vector<float> val;
  vector <int> prev;
  val.resize(nnode);
  prev.resize(nnode);
  w.resize(nnode);
  int base=0;

  float stro_normal=0;

  int loop=0;
  int lastnode=-1;


    float avgw=0;
    int avgct=0;

  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
	int idx=gnode[x][y][z];
	if (idx>=0)
	for (int ct=0;ct<alist[x][y][z].size();ct++)
	  if (alist[x][y][z][ct].w>0)
	  {	 
	  avgw+=alist[x][y][z][ct].w;
	  avgct++;
	}
      }
  avgw/=avgct;

  //trim
    if(0)
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
	int idx=gnode[x][y][z];
	if (idx>=0)
	for (int ct=0;ct<alist[x][y][z].size();ct++)
	  if (alist[x][y][z][ct].w>avgw*40){
	    alist[x][y][z][ct].w=avgw*40;
	  }
      }


    // tolgo archi deboli (nel vicinato di un nodo)
    if (1)
    for (int x=0;x<NB;x++)
      for (int y=0;y<NB;y++)
	for (int z=0;z<NB;z++){
	  int idx=gnode[x][y][z];
	  float maxw=0;
	  for (int ct=0;ct<alist[x][y][z].size();ct++){
	    float w=alist[x][y][z][ct].w;
	    int nex=gnode[alist[x][y][z][ct].x][alist[x][y][z][ct].y][alist[x][y][z][ct].z];
	    if (idx!=nex && maxw<w) maxw=w;
	  }

	  for (int ct=alist[x][y][z].size()-1;ct>=0;ct--)
	    if (alist[x][y][z][ct].w<=0.001*maxw){ // tolgo arco debole (controlla anche simmetrico, anche se non era debole)
	      // vicino
	      int nex=gnode[alist[x][y][z][ct].x][alist[x][y][z][ct].y][alist[x][y][z][ct].z];
	      char3 x1=gchar3[nex];	      
	      alist[x][y][z].erase(alist[x][y][z].begin()+ct);
	      //cerco e rimuovo arco vicino
	      for (int ct1=alist[x1.x][x1.y][x1.z].size()-1;ct1>=0;ct1--){
		if (idx==gnode[alist[x1.x][x1.y][x1.z][ct1].x][alist[x1.x][x1.y][x1.z][ct1].y][alist[x1.x][x1.y][x1.z][ct1].z])
		    alist[x1.x][x1.y][x1.z].erase(alist[x1.x][x1.y][x1.z].begin()+ct1);		    
	      }	      
	    }	  
	}


    

    
  float rangemaxg=0;
  float rangemaxgn=0;
  for (int x=0;x<NB;x++)
  for (int y=0;y<NB;y++)
    for (int z=0;z<NB;z++){
      int idx=gnode[x][y][z];
      float wt=0;
      if (idx>=0){
      for (int ct=0;ct<alist[x][y][z].size();ct++){
	if ( idx==gnode[alist[x][y][z][ct].x][alist[x][y][z][ct].y][alist[x][y][z][ct].z]) // selfloop
	  wt+=alist[x][y][z][ct].w;
	if ( idx!=gnode[alist[x][y][z][ct].x][alist[x][y][z][ct].y][alist[x][y][z][ct].z] // no selfloop
	     &&
	     rangemaxg<alist[x][y][z][ct].w){
	  rangemaxg=alist[x][y][z][ct].w;
	}
      }
      if (rangemaxgn<wt)
	rangemaxgn=wt;
      pos[x][y][z][0]+=(float)x/sx;
      pos[x][y][z][1]+=(float)y/sy;
      w[idx]=wt;
      }
    }
  printf("avg %f max %f\n",avgw,rangemaxg);


  /* // test di assorbimento grafo
    // cerco nodo minimo, cerco nodo vicino con arco a peso minimo e merge
  int lowest=-1;
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
      int idx=gnode[x][y][z];
      if (alist[x][y][z].size()>0 && w[idx]>=0 && (lowest<0 ||w[idx]<w[lowest]))
	lowest=idx;
      }
  //
  printf("low %d %f\n",lowest,w[lowest]);
  {
    float wt=0;
    int idx=lowest;
    char3 n=gchar3[idx];
    int ctmin=-1;
    if (idx>=0){
      for (int ct=0;ct<alist[n.x][n.y][n.z].size();ct++){
	if ( idx!=gnode[alist[n.x][n.y][n.z][ct].x][alist[n.x][n.y][n.z][ct].y][alist[n.x][n.y][n.z][ct].z] // no selfloop
	     &&
	     (ctmin<0 || alist[n.x][n.y][n.z][ctmin].w>alist[n.x][n.y][n.z][ct].w)){
	  ctmin=ct;
	}
      }
    }
    printf("ctmin %d\n",ctmin);
  int idx1=gnode[alist[n.x][n.y][n.z][ctmin].x][alist[n.x][n.y][n.z][ctmin].y][alist[n.x][n.y][n.z][ctmin].z];
  char3 n1=gchar3[idx1];
  printf("idx1 %d\n",idx1);
  // scarico tutti archi uscenti idx su idx1
  for (int ct=0;ct<alist[n.x][n.y][n.z].size();ct++){
    int dest=gnode[alist[n.x][n.y][n.z][ct].x][alist[n.x][n.y][n.z][ct].y][alist[n.x][n.y][n.z][ct].z];
    char3 n2=gchar3[dest];
    //cerca arco idx1-->dest
    printf("dest %d\n",dest);
    int found=-1;
    for (int ct1=0;ct1<alist[n1.x][n1.y][n1.z].size();ct1++){
      if (dest==gnode[alist[n1.x][n1.y][n1.z][ct1].x][alist[n1.x][n1.y][n1.z][ct1].y][alist[n1.x][n1.y][n1.z][ct1].z])
	found=ct1;
    }
    printf("found %d\n",found);
    if (found>=0) {
      printf("val %f\n",alist[n1.x][n1.y][n1.z][found].w);
      alist[n1.x][n1.y][n1.z][found].w+=alist[n.x][n.y][n.z][ctmin].w;
      printf("val %f\n",alist[n1.x][n1.y][n1.z][found].w);
      //cerca simmetrico (c'e') dest-->idx1
      int found=-1;
      for (int ct2=0;ct2<alist[n2.x][n2.y][n2.z].size();ct2++){
	if (idx1==gnode[alist[n2.x][n2.y][n2.z][ct2].x][alist[n2.x][n2.y][n2.z][ct2].y][alist[n2.x][n2.y][n2.z][ct2].z])
	  found=ct2;
      }
      printf("found1 %d\n",found);
      if (found>=0){
      printf("val %f\n",alist[n2.x][n2.y][n2.z][found].w);
	alist[n2.x][n2.y][n2.z][found].w+=alist[n.x][n.y][n.z][ctmin].w;
      printf("val %f\n",alist[n2.x][n2.y][n2.z][found].w);
      }
      
      //cancella arco e simmetrico idx1->idx
      printf("val %f\n",alist[n.x][n.y][n.z][ctmin].w);
      alist[n.x][n.y][n.z][ctmin].w=0;
      printf("val %f\n",alist[n.x][n.y][n.z][ctmin].w);
      found=-1;
      for (int ct2=0;ct2<alist[n1.x][n1.y][n1.z].size();ct2++){
	if (idx==gnode[alist[n1.x][n1.y][n1.z][ct2].x][alist[n1.x][n1.y][n1.z][ct2].y][alist[n1.x][n1.y][n1.z][ct2].z])
	  found=ct2;
      }
      if (found>=0){
      printf("val %f\n",alist[n1.x][n1.y][n1.z][found].w);
	alist[n1.x][n1.y][n1.z][found].w=0;
      printf("val %f\n",alist[n1.x][n1.y][n1.z][found].w);
      }
    }
    else{ // aggiungi arco idx1->dest, dest->idx1
      printf("add\n");
	  edge t;
	  t.x=n2.x;
	  t.y=n2.y;
	  t.z=n2.z;
	  t.w=alist[n.x][n.y][n.z][ct].w;
	  alist[n1.x][n1.y][n1.z].push_back(t);
	  t.x=n1.x;
	  t.y=n1.y;
	  t.z=n1.z;
	  alist[n2.x][n2.y][n2.z].push_back(t);      
    }
    // cancella nodo
    //gnode[n.x][n.y][n.z]=-1;
    
  }
  }
*/


  
  /*
  ///taglio i archi troppo deboli
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
	int idx=gnode[x][y][z];
	if (idx>=0)
	for (int ct=alist[x][y][z].size()-1;ct>=0;ct--)
	  if (alist[x][y][z][ct].w<avgw/50){
	    alist[x][y][z].erase(alist[x][y][z].begin()+ct);
	  }
      }
  
  rangemaxg=0;
  rangemaxgn=0;
  for (int x=0;x<NB;x++)
  for (int y=0;y<NB;y++)
    for (int z=0;z<NB;z++){
      int idx=gnode[x][y][z];
      float wt=0;
      if (idx>=0){
      for (int ct=0;ct<alist[x][y][z].size();ct++){
	//if ( idx==gnode[alist[x][y][z][ct].x][alist[x][y][z][ct].y][alist[x][y][z][ct].z]) // selfloop
	  wt+=alist[x][y][z][ct].w;
	if ( idx!=gnode[alist[x][y][z][ct].x][alist[x][y][z][ct].y][alist[x][y][z][ct].z] // no selfloop
	     &&
	    rangemaxg<alist[x][y][z][ct].w)
	  rangemaxg=alist[x][y][z][ct].w;
      }
      if (rangemaxgn<wt)
	rangemaxgn=wt;
      pos[x][y][z][0]+=(float)x/sx;
      pos[x][y][z][1]+=(float)y/sy;
      w[idx]=wt;
      }
    }
  printf("avg %f max %f\n",avgw,rangemaxg);
  */

  vector<int> best;
  float bestv=0;


  /*  int st=-1;
  float stv=0;
  for (int st1=0;st1<nnode;st1++){
    float v=0;
    char3 no=gchar3[st1];
    for (int ct=0;ct<alist[no.x][no.y][no.z].size();ct++){
      v+=alist[no.x][no.y][no.z][ct].w;
    }
    //printf("%d (%d %d %d): %f\n",st1,no.x,no.y,no.z,v);
    if (stv<v){
      stv=v;
      st=st1;
    }
  }
  if (st>=0)
*/

  
  //if (0)

  vector<int> nflag;
  for (int i=0;i<nnode;i++)
    nflag.push_back(0);

  long maxv=0;


  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
	C[x][y][z]=-1;
	B[x][y][z]=A[x][y][z];
	if (maxv<A[x][y][z]) maxv=A[x][y][z];
      }

  code=0;
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++)
	if (B[x][y][z]>0) {	  
	  rec_label(x,y,z);
	}

  printf("ok label\n");
  for (int i=0;i<code_map.size();i++){
    printf("codemap %d %d\n",i,code_map[i]);
  }
  
  FILE* grx=fopen("graphx.dot","w+");
  fprintf(grx,"digraph{\nnode [shape = circle, fixedsize=true, style = filled, fillcolor = palegreen];\n");
  for (int i=0;i<nnode;i++){
    char3 no=gchar3[i];
    float w1=w[i]/rangemaxgn;
    if (w1<0) w1=0;
    if (w1>1) w1=1;

    printf("%d %f %f\n",i,w[i],w1);
    //if (w1>0.001)
    fprintf(grx,"%d [width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    i,
	    pow(w1,0.5),
	    no.x*S,no.y*S,no.z*S);
  }

  for (int st1=0;st1<nnode;st1++){
    float v=0;
    char3 no=gchar3[st1];

    float maxw=0;
    for (int ct=0;ct<alist[no.x][no.y][no.z].size();ct++){
      float w=alist[no.x][no.y][no.z][ct].w;
      if (maxw<w) maxw=w;
    }
    
    for (int ct=0;ct<alist[no.x][no.y][no.z].size();ct++){
      float w=1-alist[no.x][no.y][no.z][ct].w/rangemaxg;

      /*
w=1-pow(pow(no.x-alist[no.x][no.y][no.z][ct].x,2.0f)+
	      pow(no.y-alist[no.x][no.y][no.z][ct].y,2.0f)+
	      pow(no.z-alist[no.x][no.y][no.z][ct].z,2.0f),0.5f)/255*S;
      */
      

      if (w<0) w=0;
      if (w>1) w=1;
	int nex=gnode[alist[no.x][no.y][no.z][ct].x][alist[no.x][no.y][no.z][ct].y][alist[no.x][no.y][no.z][ct].z];
	char3 no1=gchar3[nex];
	if (st1<nex)
	{
	w=pow(1-w,0.5);
	fprintf(grx,"%d -> %d [dir=none penwidth=%1.2f color = \"",st1,nex,15*(w)+0.1);
	fprintf(grx,"#%02x%02x%02x len=%f weight=%f",0,0,0,w,10*(w));
	fprintf(grx,"\"];\n");
      }
	if (0 && alist[no.x][no.y][no.z][ct].w<=0.01*maxw && st1<nex)
	{
	w=pow(w,0.5);
	w=0;
	fprintf(grx,"%d -> %d [dir=none penwidth=%1.2f color = \"",st1,nex,1*(1-w)+0.1);
	fprintf(grx,"#%02x%02x%02x len=%f weight=%f",(int)(255),0,0,w,10*(1-w));
	fprintf(grx,"\"];\n");
      }

    }
  }

fprintf(grx,"}\n");

  fclose(grx);

  
  ///////////////////////////////////////////
  // analisi aree contigue
  vector < vector< st > > nei;
  nei.resize(code);
  for (int i=0;i<code;i++)
    nei[i].resize(code);
  /*  struct st{
    long m1,M1;
    long m2,M2;
    float best; // best val di interscambio
    float col;
    int idx;
    float acc;
    float w;

    float arco;

    int x1,y1,z1;
    int x2,y2,z2;
  };*/
  
  printf("ok\n");
  for (int i=0;i<code;i++)
    for (int j=0;j<code;j++){
      nei[i][j].best=0;	  
      nei[i][j].idx=-1;	  
      nei[i][j].m1=1e10;
      nei[i][j].M1=0;
      nei[i][j].m2=1e10;
      nei[i][j].M2=0;
      nei[i][j].best=0;
      nei[i][j].acc=0;
      nei[i][j].w=0;
      nei[i][j].node=-1;
      nei[i][j].ed=-1;
    }

  // per ogni nodo scelgo il migliore arco con cui collegarmi ad altre zone
  // ogni coppia di zone ha alpiu' un arco di collegamento
  
  for (int st1=0;st1<nnode;st1++){
    float v=0;
    char3 no=gchar3[st1];
    //printf("%d/%d: %d %d %d (%d)\n",st1,nnode,no.x,no.y,no.z,alist[no.x][no.y][no.z].size());
    int best1=-1;
    int best2=-1;
    int edge=-1;
    float w1=-1;
    for (int ct1=0;ct1<alist[no.x][no.y][no.z].size();ct1++){ // per ogni arco
      int neigh=gnode[alist[no.x][no.y][no.z][ct1].x][alist[no.x][no.y][no.z][ct1].y][alist[no.x][no.y][no.z][ct1].z];
      char3 no1=gchar3[neigh];
      //printf("ok %d code %d %d\n",neigh,C[no.x][no.y][no.z],C[no1.x][no1.y][no1.z]);
      if (w[neigh]>w[st1] && C[no.x][no.y][no.z]!=C[no1.x][no1.y][no1.z]){ // solo monotono e codici diversi
	int c1=C[no.x][no.y][no.z];
	int c2=C[no1.x][no1.y][no1.z];

	if (alist[no.x][no.y][no.z][ct1].w>w1){
	  w1=alist[no.x][no.y][no.z][ct1].w;
	  best1=c1;
	  best2=c2;
	  edge=ct1;
	}
	}
      }

    // vedo se miglioro coppia proposta
    if (w1>=0 && w1>nei[best1][best2].w){
      nei[best1][best2].w=w1;
      nei[best1][best2].node=st1;
      nei[best1][best2].ed=edge;
    }
    
    }
  printf("ok---------------\n");

  // nuovo grafo
  vector<edge> edges;
  vector<int> nodes;    // 2= root, 1= normale
  nodes.resize(nnode);
  for (int i=0;i<nnode;i++){
    nodes[i]=0;
    char3 no1=gchar3[i];
    if (code_map[C[no1.x][no1.y][no1.z]]==i)  // memorizzo i nodi massimali
      nodes[i]=2;
  }

  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
	int n=gnode[x][y][z];
	for (int ct1=0;ct1<alist[x][y][z].size();ct1++){ // per ogni arco
	  int neigh=gnode[alist[x][y][z][ct1].x][alist[x][y][z][ct1].y][alist[x][y][z][ct1].z];	  
	  edge t;
	  t.x=n;
	  t.y=neigh;
	  t.w=alist[x][y][z][ct1].w;
	  edges.push_back(t);
	  if (nodes[n]!=2)
	  nodes[n]=1;
	  if (nodes[neigh]!=2)
	  nodes[neigh]=1;
	}
      }
  

  int maxim=0;
  for (int i=0;i<nnode;i++){
    if (nodes[i]==2)
      maxim++;
  }
  printf("maxim %d\n",maxim);

  rangemaxgn=0;
  for (int i=0;i<nnode;i++)
    if (nodes[i] && w[i]>rangemaxgn)
      rangemaxgn=w[i];
  rangemaxg=0;

  vector<float> wlist;
  float rangemaxgavg=0;
  int cte=0;
  for (int i=0;i<edges.size();i++)
    if (edges[i].x!=edges[i].y) {
      cte++;
      wlist.push_back(edges[i].w);
      //printf("%d %d %f\n",edges[i].x,edges[i].y,edges[i].w);
      if (edges[i].w>rangemaxg)
	rangemaxg=edges[i].w;
      rangemaxgavg+=edges[i].w;
    }
  rangemaxgavg/=cte;


  sort (wlist.begin(), wlist.begin()+wlist.size());
  float w20=wlist[(int)(wlist.size()*0.9)]; // mediana al 90% delle intensita' di arco

  /*
  // calcolo peso sulla media per selezionare
  // 20% degli edge a massimo
  int ct=0;
  int ct1=0;
  for (int i=0;i<edges.size();i++){
    if (edges[i].x!=edges[i].y && edges[i].w>w20)
      ct++;
    if (edges[i].x!=edges[i].y)
      ct1++;
  }
  printf("tot %f %d\n",((float)ct)/((float)ct1),ct1);
*/
    


  
  // edges raddoppiati
  int ned=edges.size();
  vector <float> nneigh;
  vector <float> nneigh1;
  nneigh.resize(nnode);
  nneigh1.resize(nnode);

  for (int i=0;i<nnode;i++)
    if (nodes[i]){
      int ct=0;
      char3 no=gchar3[i];

      for (int x=-1;x<=1;x++)
	for (int y=-1;y<=1;y++)
	  for (int z=-1;z<=1;z++)
	    if (abs(x)+abs(y)+abs(z)<=1)
	    if (no.x+x>=0 && no.x+x<NB &&
		no.y+y>=0 && no.y+y<NB &&
		no.z+z>=0 && no.z+z<NB){
	      //printf("%d %d %d: %d\n",no.x+x,no.y+y,no.z+z,A[no.x+x][no.y+y][no.z+z]);
	      //if (A[no.x+x][no.y+y][no.z+z]==0)
		ct+=A[no.x+x][no.y+y][no.z+z];
	    }

      if (0)
      for (int i1=0;i1<nnode;i1++){
	  char3 no1=gchar3[i1];
	  
	  if (abs(no.x-no1.x)<=1 &&
	      abs(no.y-no1.y)<=1 &&
	      abs(no.z-no1.z)<=1 && nodes[i1]==0)
	    ct=0;	 
	}
      
      nneigh[i]=1+ct;///A[no.x][no.y][no.z];
    }

  int loo=0;
  //while(loo++<2)
  if (0)
    {
      for (int i=0;i<nnode;i++)
	if (nodes[i]){
	  float avg=0;
	  int ct1=0;
	  int ct2=0;
	  char3 no=gchar3[i];
	  for (int i1=0;i1<nnode;i1++)
	    if (nodes[i1])
	    {
	      char3 no1=gchar3[i1];
	  
	      if (abs(no.x-no1.x)<=1 &&
		  abs(no.y-no1.y)<=1 &&
		  abs(no.z-no1.z)<=1){
		if (nneigh[i]<nneigh[i1]) // vicino piu' alto
		  ct1++;
		if (nneigh[i]>nneigh[i1])
		  ct2++;
		avg+=nneigh[i1];
	      }
	    }
	  if (ct1>ct2) // sono piu' basso --> estremizzo
	    nneigh1[i]=0.95*nneigh[i];
	  if (ct1<ct2) // sono piu' basso --> estremizzo
	    nneigh1[i]=1.05*nneigh[i];
	}

      for (int i=0;i<nnode;i++){
	nneigh[i]=nneigh1[i];
      }

      float maxn=0;  
      for (int i=0;i<nnode;i++)
	if (maxn<nneigh[i]) maxn=nneigh[i];
  
      for (int i=0;i<nnode;i++){
	//	printf("%d %f\n",i,nneigh[i]);
	nneigh[i]/=maxn;
	//printf("%d %f\n",i,nneigh[i]);
      }
    }
    
  if (0)
  for (int i=0;i<ned;i++){
    edge t;
    t.x=edges[i].y;
    t.y=edges[i].x;
    t.w=edges[i].w;
    edges.push_back(t);
  }

  //////////////////////////// analisi frequenza archi consecutivi

  // segui percorso su linea
  /*  int p1x=0;
  int p2x=sx;
  int p1y=sy/2;
  int p2y=sy/2;
  */
    for (int dir=0;dir<2;dir++){
      
  th=divisioni;//sx;
  if (dir==0){
    if (th>sx) th=sx;
  }
  else{
    if (th>sy) th=sy;
  }
  step=1;//+sx/50;
  for (int q=0;q<th;q+=step){
    int p1x,p1y,p2x,p2y;
    if (dir==0){
      p1x=sx*q/th;
      p1y=0;
      p2x=sx*q/th;
      p2y=sy-1;
    }
    else{
      p1x=0;
      p1y=sy*q/th;
      p2x=sx-1;
      p2y=sy*q/th;
    }
    printf("linea %d %d %d %d\n",p1x,p1y,p2x,p2y);
    vector<int> seq;
    
  //  for (int x=0;x<sx;x+=1){
  //  for (int y=0;y<sy;y+=1){
  int steps=abs(p1x-p2x);
  if (steps<abs(p1y-p2y))
    steps=abs(p1y-p2y);
  int prev=-1;
  int curr;
  for (int s=0;s<steps;s++) {
    float s1=(float)s/steps;
    float s2=(float)(s+1)/steps;
    int x=(int)round(s1*p2x+(1-s1)*p1x);
    int y=(int)round(s1*p2y+(1-s1)*p1y);
    int x1=(int)round(s2*p2x+(1-s2)*p1x);
    int y1=(int)round(s2*p2y+(1-s2)*p1y);
    //    printf("%d %d %f, %d %d - %d %d\n",s,steps,(s1*p1y+(1-s1)*p2y),x,y,x1,y1);
    unsigned char r,g,b;
    unsigned char r1,g1,b1;
    r=src(x,y,0,0);
    g=src(x,y,0,1);
    b=src(x,y,0,2);
      
      r1=src(x1,y1,0,0);
      g1=src(x1,y1,0,1);
      b1=src(x1,y1,0,2);

      /// linea di colore
      int stepsc=abs(r-r1);
      if (stepsc<abs(g-g1))
	stepsc=abs(g-g1);
      if (stepsc<abs(b-b1))
	stepsc=abs(b-b1);
  for (int sc=0;sc<stepsc;sc++) {
    float sc1=(float)sc/stepsc;
    float sc2=(float)(sc+1)/stepsc;
    int r2=(int)round(sc1*r1+(1-sc1)*r);
    int g2=(int)round(sc1*g1+(1-sc1)*g);
    int b2=(int)round(sc1*b1+(1-sc1)*b);
    int r3=(int)round(sc2*r1+(1-sc2)*r);
    int g3=(int)round(sc2*g1+(1-sc2)*g);
    int b3=(int)round(sc2*b1+(1-sc2)*b);

    
    if (sc==0)
      curr=gnode[r3/S][g3/S][b3/S];
    curr=gnode[r2/S][g2/S][b2/S];
    
    if (curr!=prev){
      seq.push_back(curr);
      prev=curr;
    }
   
    }
  }

  // statistica per ogni arco entrante, quale e' la prob arco uscente, su cui basare la ricerca di percorsi probabili
  
  for (int i=0;i<(int)seq.size()-2;i++) {
    //    if (seq[i]==157)
    //    printf("%d -> %d -> %d;\n",seq[i],seq[i+1],seq[i+2]);
      // trova edges
      int found1=-1;
      int found2=-1;
      for (int j=0;j<edges.size();j++){
	if (edges[j].x==seq[i] && edges[j].y==seq[i+1])
	  found1=j;
	if (edges[j].x==seq[i+1] && edges[j].y==seq[i+2])
	  found2=j;
      }
      //cerco found2 in found1 (aggiungo statistica)
      int idx=-1;
      if (found1>=0 && found2>=0){
      for (int j=0;j<edges[found1].nedge.size();j++)
	if (edges[found1].nedge[j]==found2)
	  idx=j;
      if (idx==-1){ // non c'e'
	edges[found1].nedge.push_back(found2);
	edges[found1].wedge.push_back(1);	
      }
      else{
	edges[found1].wedge[idx]++;		
      }
      }
      // percorso invertito
      found1=-1;
      found2=-1;
      for (int j=0;j<edges.size();j++){
	if (edges[j].x==seq[i+2] && edges[j].y==seq[i+1])
	  found1=j;
	if (edges[j].x==seq[i+1] && edges[j].y==seq[i])
	  found2=j;
      }
      //cerco found2 in found1 (aggiungo statistica)
      idx=-1;
      if (found1>=0 && found2>=0){
      for (int j=0;j<edges[found1].nedge.size();j++)
	if (edges[found1].nedge[j]==found2)
	  idx=j;
      if (idx==-1){ // non c'e'
	edges[found1].nedge.push_back(found2);
	edges[found1].wedge.push_back(1);	
      }
      else{
	edges[found1].wedge[idx]++;		
      }
      }
	
  }
  }
}
  
  // normalizza stats
  for (int i=0;i<edges.size();i++)
    {
    float sum=0;
    printf("edge %d -> %d: ",edges[i].x,edges[i].y);
    for (int j=0;j<edges[i].nedge.size();j++){
      sum+=edges[i].wedge[j];
    }
    for (int j=0;j<edges[i].nedge.size();j++){
      if (sum>0)
      edges[i].wedge[j]/=sum;
      printf(" %d -> %d (%f) ",edges[edges[i].nedge[j]].x,edges[edges[i].nedge[j]].y,edges[i].wedge[j]);
    }
    printf("\n");
  }

  /*
  char3 nodbg=gchar3[120];
  printf("node %d: %d %d %d\n",120,nodbg.x,nodbg.y,nodbg.z);
  nodbg=gchar3[157];
  printf("node %d: %d %d %d\n",157,nodbg.x,nodbg.y,nodbg.z);
  */
  
  grx=fopen("graphx1.dot","w+");
  fprintf(grx,"digraph{\nnode [shape = circle, fixedsize=true, style = filled, fillcolor = palegreen];\n");
  for (int i=0;i<nnode;i++)
    if (nodes[i])
    {

      /// uso colore massimo
    char3 no1=gchar3[i];
    char3 no=gchar3[code_map[C[no1.x][no1.y][no1.z]]];
    no=gchar3[i];
    
    float w1=w[i]/rangemaxgn;
    if (w1<0) w1=0;
    if (w1>1) w1=1;
    
    if (nodes[i]==3)
      fprintf(grx,"%d [shape = doublecircle ",i);
    else{
      if (nodes[i]==2)
	fprintf(grx,"%d [shape = doublecircle ",i);
      else
	fprintf(grx,"%d [",i);
    }
    int v=(int)(255*nneigh[i]);
    
    fprintf(grx,"width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    pow(w1,0.2),
	    //	    	    v,v,v);
	    no.x*S,no.y*S,no.z*S);
    }


  for (int i=0;i<edges.size();i++)
    if (edges[i].x>edges[i].y && edges[i].nedge.size()>0){
      float w=edges[i].w/w20;
      if (w<0.0) w=0.0;
      if (w>1) w=1;
    
      w=pow(w,0.3);
      fprintf(grx,"%d -> %d [dir=none penwidth=%1.2f color = \"",edges[i].x,edges[i].y,10*(w)+0.1);
      fprintf(grx,"#%02x%02x%02x len=%f weight=%f",0,0,0,w,10*(w));
      fprintf(grx,"\"];\n");
    }
    
  fprintf(grx,"}\n");
  fclose(grx);

  /*
    ///// per ora non lo uso -> lavoro direttamente con grafo semplificato

  // grafo semplificato
  // codemap sono i nodi
  float newedges[code][code];
  float neww[code];
  for (int i=0;i<code;i++)
    neww[i]=0;
  for (int i=0;i<nnode;i++){
    char3 no=gchar3[i];
    int node=C[no.x][no.y][no.z];
    neww[node]+=w[i];
  }
  float maxw=0;
  for (int i=0;i<code;i++){
    printf("node %d: %f\n",code_map[i],neww[i]);
    if (maxw<neww[i])
      maxw=neww[i];
  }
  
  for (int i=0;i<code;i++)
    for (int j=0;j<code;j++){
      newedges[i][j]=0;
    }
  
  for (int i=0;i<edges.size();i++){ // prendo da vecchi edges
    int a=edges[i].x;
    int b=edges[i].y;
    char3 noa=gchar3[a];
    char3 nob=gchar3[b];
    int nodea=C[noa.x][noa.y][noa.z];
    int nodeb=C[nob.x][nob.y][nob.z];
    float v=edges[i].w;
    newedges[nodea][nodeb]+=v;
    newedges[nodeb][nodea]+=v;
  }

  vector<float> newwlist;
  float newrangemaxgavg=0;
  int newcte=0;
  for (int i=0;i<code;i++)
    for (int j=i+1;j<code;j++)
      if (newedges[i][j]>0){
      newcte++;
      newwlist.push_back(newedges[i][j]);
      newrangemaxgavg+=newedges[i][j];
    }
  newrangemaxgavg/=newcte;


  sort (newwlist.begin(), newwlist.begin()+newwlist.size());
  float neww20=newwlist[(int)(newwlist.size()*0.9)]; // mediana al 90% delle intensita' di arco

  
  grx=fopen("graphy1.dot","w+");
  fprintf(grx,"digraph{\nnode [shape = circle, fixedsize=true, style = filled, fillcolor = palegreen];\n");
  for (int i=0;i<code_map.size();i++)    {

      /// uso colore massimo
    char3 no1=gchar3[code_map[i]];
    char3 no=gchar3[code_map[C[no1.x][no1.y][no1.z]]];

    
    float w1=neww[i]/maxw;
    if (w1<0) w1=0;
    if (w1>1) w1=1;
    
    fprintf(grx,"%d [",code_map[i]);
    fprintf(grx,"width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    pow(w1,0.2),
	    no.x*S,no.y*S,no.z*S);
    }

  printf("neww20 %f\n",neww20);
  for (int i=0;i<code;i++)
    for (int j=i+1;j<code;j++)
      if (newedges[i][j]>0){
      printf("New edges %d %d: %f\n",code_map[i],code_map[j],newedges[i][j]);
      float w=newedges[i][j]/neww20;
      if (w<0.0) w=0.0;
      if (w>1) w=1;
    
      w=pow(w,0.3);
      fprintf(grx,"%d -> %d [dir=none penwidth=%1.2f color = \"",code_map[i],code_map[j],10*(w)+0.1);
      fprintf(grx,"#%02x%02x%02x len=%f weight=%f",0,0,0,0,0);//w,10*(w));
      fprintf(grx,"\"];\n");
    }
    
  fprintf(grx,"}\n");
  fclose(grx);


  /// sostituisco grafo semplificato!
  for (int i=0;i<nnode;i++){
    w[i]=0;
    nodes[i]=0;
  }
  for (int i=0;i<code;i++){
    w[code_map[i]]=neww[i];
    nodes[code_map[i]]=2;
  }
  edges.clear();
  for (int i=0;i<code;i++)
    for (int j=0;j<code;j++)
      if (newedges[i][j]>0 && i!=j){
	    edge t;
	    t.x=code_map[i];
	    t.y=code_map[j];
	    t.w=newedges[i][j];
	    edges.push_back(t);
      }
  // ricalcolo limiti
  rangemaxgn=0;
  for (int i=0;i<nnode;i++)
    if (nodes[i] && w[i]>rangemaxgn)
      rangemaxgn=w[i];
  rangemaxg=0;

  wlist.clear();
  rangemaxgavg=0;
  cte=0;
  for (int i=0;i<edges.size();i++)
    if (edges[i].x!=edges[i].y) {
      cte++;
      wlist.push_back(edges[i].w);
      printf("%d %d %f\n",edges[i].x,edges[i].y,edges[i].w);
      if (edges[i].w>rangemaxg)
	rangemaxg=edges[i].w;
      rangemaxgavg+=edges[i].w;
    }
  rangemaxgavg/=cte;


  sort (wlist.begin(), wlist.begin()+wlist.size());
  w20=wlist[(int)(wlist.size()*0.90)]; // mediana al 90% delle intensita' di arco
*/


  // distribuisco su adj
  vector< vector< edge > > adj;
  adj.resize(nnode);

  for (int i=0;i<edges.size();i++){
    adj[edges[i].x].push_back(edges[i]);
  }
  

  

  val.resize(nnode);
  prev.resize(nnode);

  int ngr=-1;
  vector< vector <int> > lists;  
  vector< float > listsw;
  
    for (int st=0;st<nnode;st++)
      if (nodes[st]==2) // solo root
    {
      ngr++;
      
      printf("start %d\n",st);
    fheap_t* H=fh_alloc(nnode);
    
    
    for (long i1=0;i1<nnode;i1++){
      val[i1] = INFTY;
      prev[i1]=-1;
    }
    
    int start=st;
    val[start]= 1;
    fh_insert(H, start, val[start]);
    int dbg=st==445;

  while (H->n>0){
    long u=-1;  
    u=fh_delete_min(H);
    if (dbg) printf("u %d\n",u);

    float localmaxw=0;
    for (int ct=0;ct<adj[u].size();ct++){
      int i1=adj[u][ct].y;
      float v=adj[u][ct].w;
      if (i1!=u && localmaxw<v)
	localmaxw=v;
    }

    if (prev[u]==-1){ // esploro tutti gli archi uscenti (primo nodo
    for (int ct=0;ct<adj[u].size();ct++){
      int i1=adj[u][ct].y;
      float v=adj[u][ct].w;
      if (dbg) printf(" primo --> %d %f\n",i1,v);
      if (i1!=u && v>0){
	char3 no1=gchar3[st];
	char3 no2=gchar3[i1];
	
	char3 x1=gchar3[i1];
	char3 x2=gchar3[u];	      
	float dist=pow((x1.x-x2.x)*(x1.x-x2.x)+
		  (x1.y-x2.y)*(x1.y-x2.y)+
		  (x1.z-x2.z)*(x1.z-x2.z),0.5);
	dist=abs(x1.x-x2.x)+abs(x1.y-x2.y)+abs(x1.z-x2.z);
	
	// w=0 --> veloce, w=1 lento
	float v1=1-v/(w20);
	if (v1<0.0) v1=0.0;
	if (v1>1) v1=1;
	v1=pow(v1,3);
	//	float w=val[u]*(2-pow(v/rangemaxg,0.02))+(1-v/rangemaxg)+0.0001;
	float w;
	//w=val[u]*(1+0.5*v1)+0.001;//+dist/v1;
	w=val[u]*(1+0.1*v1)+0.1*dist*(1+v1);//+dist/v1;
	w=val[u]+0.001*dist*(1+v1);//+dist/v1;
	w=val[u]+0.001;
	
	if (w>=INFTY/2) w=INFTY/2;
	//printf("%f %f %f, %f\n",val[u],v1,w,val[i1]);
	
		// se peggioro, uso nuovo valore
	if (dbg) printf("val %f, newval %f i1 %d u %d\n",val[i1],w,i1,u);
	float newval=w;
	if (val[i1]==INFTY){
	  fh_insert(H, i1, newval);
	  val[i1] = newval;
	  prev[i1]=u;		
	}
	else{
	  if (val[i1]>newval){
	    val[i1] = newval;
	    prev[i1]=u;		
	    fh_decrease_key(H, i1, newval);											
	  }
	}
      }
    }
    }// espando tutti vicini
    else{ /// c'e' prev -> uso statistiche
      
      // cerco arco
      int found=-1;
      for (int ct=0;ct<adj[prev[u]].size();ct++){
	if (adj[prev[u]][ct].y==u)
	  found=ct; 
      }
      
      if (dbg) printf("found %d\n",found);
      if (found>=0)
    for (int ct=0;ct<adj[prev[u]][found].nedge.size();ct++){
      int idxedge=adj[prev[u]][found].nedge[ct];
      int i1=edges[idxedge].y; // successivo
      float v=adj[prev[u]][found].wedge[ct]; // proababilita
      if (dbg) printf(" %d %d --> %d %f\n",prev[u],edges[idxedge].x,edges[idxedge].y,v);
      if (i1!=u && v>0){
	char3 no1=gchar3[st];
	char3 no2=gchar3[i1];
	
	char3 x1=gchar3[i1];
	char3 x2=gchar3[u];	      
	float dist=pow((x1.x-x2.x)*(x1.x-x2.x)+
		  (x1.y-x2.y)*(x1.y-x2.y)+
		  (x1.z-x2.z)*(x1.z-x2.z),0.5);
	dist=abs(x1.x-x2.x)+abs(x1.y-x2.y)+abs(x1.z-x2.z);
	
	// w=0 --> veloce, w=1 lento
	float v1=1-v;
	float w;
	w=val[u]+0.1*dist*(1+4*v1);//*(1-edges[idxedge].w)/rangemaxg;
	w=val[u]*(1+0.5*v1)+0.1*dist;
	
	if (w>=INFTY/2) w=INFTY/2;
	//printf("%f %f %f, %f\n",val[u],v1,w,val[i1]);
	
		// se peggioro, uso nuovo valore
	if (dbg) printf("val %f, newval %f, i1 %d u %d\n",val[i1],w,i1,u);
	float newval=w;
	if (val[i1]==INFTY){
	  fh_insert(H, i1, newval);
	  val[i1] = newval;
	  prev[i1]=u;		
	}
	else{
	  if (val[i1]>newval){
	    val[i1] = newval;
	    prev[i1]=u;		
	    fh_decrease_key(H, i1, newval);											
	  }
	}
      }
    }
      
    }




    
  }

  printf("ok\n");
  
  vector <int> leaf;
  leaf.resize(nnode);
  for (int i1=0;i1<nnode;i1++)
    leaf[i1]=1;
  for (int i1=0;i1<nnode;i1++)
    if (prev[i1]!=-1)
      leaf[prev[i1]]=0;
  
  vector <float> dri;
  dri.resize(nnode);
  
  // stime
    {
   for (int i1=0;i1<nnode;i1++)
     //if (leaf[i1])
       {
       int u=i1;
       int ct=0;
       float dist=0;  // linea d'aria
       float distt=0; // percorsa
       int root=-1;
       float wpath=0;
       while (u!=-1){

	 char3 x1=gchar3[i1];
	 char3 x2=gchar3[u];	      
	 // cerca peso edge
	 for (int ct1=0;ct1<adj[u].size();ct1++)
	   if (adj[u][ct1].y==prev[u]){
	     wpath+=adj[u][ct1].w;
	     char3 x3=gchar3[adj[u][ct1].y];
	     distt+=pow((x3.x-x2.x)*(x3.x-x2.x)+
			(x3.y-x2.y)*(x3.y-x2.y)+
			(x3.z-x2.z)*(x3.z-x2.z),0.5);		  
	   }
	 dist=pow((x1.x-x2.x)*(x1.x-x2.x)+
		  (x1.y-x2.y)*(x1.y-x2.y)+
		  (x1.z-x2.z)*(x1.z-x2.z),0.5);		  
	 ct++;
	 root=u;
	 u=prev[u];
       }

       float dritto=0;
       if (root!=-1){
	 char3 no1=gchar3[root];
	 char3 no2=gchar3[i1];
	// da ottimizz.
	float vec1[3];
	vec1[0]=no2.x-no1.x;
	vec1[1]=no2.y-no1.y;
	vec1[2]=no2.z-no1.z;
	float norm=pow(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2],0.5f);
	vec1[0]/=norm;vec1[1]/=norm;vec1[2]/=norm;
	int tempn=i1;
	float w1=0.01;
	int ctl=0;
	while (tempn!=-1)
	  {
	    char3 no3=gchar3[tempn];
	    float vec2[3];
	    vec2[0]=no3.x-no1.x;
	    vec2[1]=no3.y-no1.y;
	    vec2[2]=no3.z-no1.z;
	    float d=vec2[0]*vec1[0]+vec2[1]*vec1[1]+vec2[2]*vec1[2];
	    float norm2=vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2];
	    float dist=norm2-d*d;
	    if (dist<0) dist=0;
	    //	    norm=1;
	    //	    if (st==38) printf("   %d (%3.3f %f) \n",tempn,proj,dist);
	    //	    if (w<proj/norm) w=proj/norm;

	    dist=pow(dist,0.5);
	    if (w1<dist)
	      w1=dist;
	    
	      //w+=dist;
	      
	    if (0&& st==445)
	      printf("%d %f, %d %d %d, %d %d %d, %d %d %d\n",tempn,dist,no1.x,no1.y,no1.z,no3.x,no3.y,no3.z,no2.x,no2.y,no2.z);
	    //if (w<proj) w=proj;
	    tempn=prev[tempn];
	    ctl++;
	  }
	//	dritto=pow(w/ctl,0.5)/norm;
	dritto=pow(w1,5)/norm;
	
	if (st==38) printf("  w %f n %f d %f\n",pow(w1/ctl,0.5),norm,dritto);
	dritto=dritto;
	
       }
       
       dritto=1+dritto; //evito che molto dritto appiattisca a 0 lo score
       
       if (distt==0)
	 dri[i1]=INFTY;
       else
	 dri[i1]=wpath/rangemaxg/pow(w[i1],0.5)*(dritto)/dist; // divido anche per peso nodo (sotto radice per non dare troppo peso ai picchi forti)
       //	 dri[i1]=val[i1]*(dritto)/dist/pow(w[i1],0.5); // divido anche per peso nodo (sotto radice per non dare troppo peso ai picchi forti)

       if ( st==40)
	 printf("leaf %d: sp %f, le %d, dist %f %f curv %f, dri %f score %f\n",i1,val[i1]/distt,ct,dist,distt,distt/dist,dritto,val[i1]/distt*dritto);
       //       printf("leaf %d: sp %f, w %f, le %d, dist %f %f curv %f, dri %f score %f\n",i1,val[i1]/distt,w/ct,ct,dist,distt,distt/dist,dritto,w/ct/val[i1]*distt);
  }
  }

    
    //
  float drimin=1000;
  float drimax=0;
  for (int i=0;i<nnode;i++)
    if (nodes[i]&& dri[i]<INFTY){
      if (drimin>dri[i]) drimin=dri[i];
      if (drimax<dri[i] && !(dri[i]!=dri[i])) drimax=dri[i];
    }

    
    int mod=1;
    //if (st!=445)
    //    if (0)
    if (1)
    while (mod){
      mod=0;

// toglie foglia se prev e' migliore
      
      for (int i1=0;i1<nnode;i1++)
	if (leaf[i1] && prev[i1]!=-1 && dri[i1]>dri[prev[i1]] && nodes[i1]!=2){
	  prev[i1]=-1;
	  mod=1;
	}
    
      // filter
      for (int i1=0;i1<nnode;i1++)
	leaf[i1]=1;
      for (int i1=0;i1<nnode;i1++)
	if (prev[i1]!=-1)
	  leaf[prev[i1]]=0;

 // toglie foglia se esiste un nodo piu' in alto molto meglio
      if (1)
      for (int i1=0;i1<nnode;i1++)
	if (leaf[i1] && nodes[i1]!=2){
	  int dbg=st==40;
	  float val=dri[i1];
	  int u=i1;
	  while (u!=-1){
	    if (dbg) printf("%d %f %f\n",u,dri[u],dri[u]-drimin);
	    if (val>dri[u])
	      val=dri[u];
	    u=prev[u];
	  }

	  if (dbg) printf("******** check %f %f\n",dri[i1],val);
	  if ((dri[i1]) > 1.2* (val)){
	    prev[i1]=-1;
	    mod=1;
	  }
	}

      // filter
      for (int i1=0;i1<nnode;i1++)
	leaf[i1]=1;
      for (int i1=0;i1<nnode;i1++)
	if (prev[i1]!=-1)
	  leaf[prev[i1]]=0;

    ///// per ogni coppia foglie calcola similarita'
      /// colineari <0.9 angolo
      /// tolgo piu' corta! (
      if (1)
      for (int i1=0;i1<nnode;i1++)
      if (leaf[i1]){
	// cerca root
	vector<int> path1;
	int r1=i1;
	while (r1!=-1 && prev[r1]!=-1){
	  path1.push_back(r1);
	  r1=prev[r1];
	}
	if (r1!=-1 && r1!=i1)
	  for (int i2=0;i2<nnode;i2++)
	    if (i1!=i2 && leaf[i1] &&  leaf[i2]){
	      int dbg=0;
	      if (i1==731 || i2==731) dbg=1;
	      vector<int> path2;
	      int r2=i2;
	      while (r2!=-1 && prev[r2]!=-1){
		path2.push_back(r2);
		r2=prev[r2];
	      }
	      if (r2!=-1 && r2!=i2) {

		// check rette i1-r1  e i2-r2

		char3 na1=gchar3[i1];
		char3 na2=gchar3[r1];
		char3 nb1=gchar3[i2];
		char3 nb2=gchar3[r2];
		float v1[3];
		float v2[3];
		v1[0]=na1.x-na2.x;
		v1[1]=na1.y-na2.y;
		v1[2]=na1.z-na2.z;
		v2[0]=nb1.x-nb2.x;
		v2[1]=nb1.y-nb2.y;
		v2[2]=nb1.z-nb2.z;
		float norma1=norm(v1);
		float norma2=norm(v2);
		normal(v1);
		normal(v2);
		float scal=(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]); // non metto abs, cosi' se sono due percorsi opposti li tengo

		char3 no3=nb1;
		float vec2[3];
		vec2[0]=no3.x-na2.x;
		vec2[1]=no3.y-na2.y;
		vec2[2]=no3.z-na2.z;
		float d=vec2[0]*v1[0]+vec2[1]*v1[1]+vec2[2]*v1[2];
		float norm2=vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2];
		float dist=norm2-d*d;
		if (dist<0) dist=0;
		float proj1=pow(dist,0.5f);

		no3=nb2;
		vec2[0]=no3.x-na2.x;
		vec2[1]=no3.y-na2.y;
		vec2[2]=no3.z-na2.z;
		d=vec2[0]*v1[0]+vec2[1]*v1[1]+vec2[2]*v1[2];
		norm2=vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2];
		dist=norm2-d*d;
		if (dist<0) dist=0;
		float proj2=pow(dist,0.5f);

		int delet=(scal>0.90 &&  (proj1<2 && proj2<2));
		if (dbg)
		  printf("%d - %d --a-- %d - %d: %f ( %d %d %d, %d %d %d - %d %d %d, %d %d %d, %f %f %f - %f %f %f, %f %f, val %f %f)\n",i1,r1,i2,r2,scal,
			 na1.x,na1.y,na1.z,
			 na2.x,na2.y,na2.z,
			 nb1.x,nb1.y,nb1.z,
			 nb2.x,nb2.y,nb2.z,
			 v1[0],v1[1],v1[2],
			 v2[0],v2[1],v2[2],
			 proj1,
			 proj2,
			 dri[i1],
			 dri[i2]);
		if (delet && dri[i1]>dri[i2]){ // norma2>norma1){
		  if (dbg) printf("del %d\n",i1);
		  prev[i1]=-1;
		  mod=1;
		}
		if (delet && dri[i2]>dri[i1]){//norma1>norma2){
		  if (dbg) printf("del %d\n",i2);
		  prev[i2]=-1;
		  mod=1;
		}

      // se condivide troppo path in proporzione, toglie
		/// check path share;
		
		float best=100;
		for (int i=0;i<path1.size();i++){
		  char3 n1=gchar3[path1[i]];
		  char3 n2=gchar3[i2];
		  float val=pow((n1.x-n2.x)*(n1.x-n2.x)+(n1.y-n2.y)*(n1.y-n2.y)+(n1.z-n2.z)*(n1.z-n2.z),0.5f);		  
		  if (val<best) best=val;
		}
		if (dbg) printf("dist %f, ",best);
		if (best<3&& dri[i2]>dri[i1]) {
		  prev[i2]=-1;
		  mod=1;		  
		}

		best=100;
		for (int i=0;i<path2.size();i++){
		  char3 n1=gchar3[path2[i]];
		  char3 n2=gchar3[i1];
		  float val=pow((n1.x-n2.x)*(n1.x-n2.x)+(n1.y-n2.y)*(n1.y-n2.y)+(n1.z-n2.z)*(n1.z-n2.z),0.5f);		  
		  if (val<best) best=val;
		}
		if (dbg) printf("dist %f\n",best);
		if (best<3 && dri[i2]<dri[i1]) {
		  prev[i1]=-1;
		  mod=1;		  
		}

		
		/*
		int shared=0;
		int goon=1;
		for (;goon && shared<path1.size() && shared<path2.size();shared++){
		  if (path1[path1.size()-1-shared]!=path2[path2.size()-1-shared])
		    goon=0;
		}
		shared--;
		if (dbg) printf("sh %d, %f %f\n",shared,(0.0+shared)/path1.size(),(0.0+shared)/path2.size());
		if ((0.0+shared)/path1.size()>0.9){
		  prev[i1]=-1;
		  mod=1;		  
		}
		if ((0.0+shared)/path2.size()>0.9){
		  prev[i2]=-1;
		  mod=1;		  
		}
*/		
	      }
	    }
      }

      // filter
      for (int i1=0;i1<nnode;i1++)
	leaf[i1]=1;
      for (int i1=0;i1<nnode;i1++)
	if (prev[i1]!=-1)
	  leaf[prev[i1]]=0;
    
    }// loop filter


    
  bestv=0;
  //risultati
  for (int i1=0;i1<nnode;i1++)
    if (leaf[i1] && prev[i1]!=-1) {
    int u=i1;//end[i1];
    int ct=1;
    float len=0;
    float ww=0;
    lists.resize(lists.size()+1);
    listsw.resize(listsw.size()+1);
    len=val[u];//w[u];//col[u].w/max1;
    float www=dri[u];//1-(dri[u]-0.4);
    if (www<0) www=0;
    while (u!=-1 && prev[u]!=-1){
      //if (nodes[u]==2)
	{
	lists[lists.size()-1].push_back(u);	
      }
	ww+=w[u]/rangemaxgn;
      u=prev[u];
      ct++;
    }
    if (u!=i1){// && nodes[u]==2)
      printf("leaf %d --> %d %f\n",u,i1,www);
      //len+=w[u];
      ww+=w[u]/rangemaxgn;
      lists[lists.size()-1].push_back(u);
      if (ww>0) ww=1/ww;
      listsw[lists.size()-1]= ww; // serve solo per avere 0 se non c'e' peso
    }

    if (ct>2)
    {
      float v=len;//val[i1];//pow((float)ct,(float)1);//*len;//ct;//pow(2,val[end[i1]])/ct
      //      printf("best %f, %f %f %d\n",bestv,v,len,ct);
      if (bestv< v){
	bestv=v;
	best.clear();
	int u=i1;//end[i1];
	lastnode=i1;
	int ct=0;
	while (u!=-1){
	  best.insert(best.begin(),u);//push_back(u);
	  u=prev[u];
	}
      }
    }
  }


      
  
  //outgraph
  char name[256];
  char num[256];
  strcpy(name,"graphsp");
  sprintf(num,"%d",st);
  strcat(name,num);
  strcat(name,".dot");
    grx=fopen(name,"w+");
  fprintf(grx,"digraph{\nnode [shape = circle, fixedsize=true, style = filled, fillcolor = palegreen];\n");
  for (int i=0;i<nnode;i++)
    if (nodes[i])
    {
    char3 no=gchar3[i];
    float w1=w[i]/rangemaxgn;
    w1=1-val[i]/30;

    w1=1-(dri[i]-drimin)/(drimax-drimin);
    if (st==691) printf(" %f %f %f %f\n",w1,dri[i], drimin, drimax);
    
    if (w1<0) w1=0;
    if (w1>1) w1=1;
    w1=pow(w1,3);
    
      if (nodes[i]==2)
    fprintf(grx,"a%da [shape = doublecircle width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    i,
	    pow(w1,0.9)+0.2,
	    no.x*S,no.y*S,no.z*S);
    else
    fprintf(grx,"a%da [width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    i,
	    pow(w1,0.9)+0.2,
	    no.x*S,no.y*S,no.z*S);
  }


  for (int i=0;i<nnode;i++){
    if (prev[i]!=-1){
      int found=-1;
      for (int j=0;j<edges.size();j++)
	if (edges[j].x==i && edges[j].y==prev[i])
	  found=j;
      if (found>=0){
	float w=1-edges[found].w/rangemaxgavg;
	//      printf("%d: %f --> %f\n",i,edges[found].w,w);
      if (w<0.0) w=0.0;
      w=pow(1-w,0.5);
      fprintf(grx,"a%da -> a%da [penwidth=%1.2f color = \"",edges[found].y,edges[found].x,5*(w)+0.1);
      fprintf(grx,"#%02x%02x%02x len=%f weight=%f",0,0,0,w,10*(w));
      fprintf(grx,"\"];\n");
    }
    }
  }
  fprintf(grx,"}\n");
  fclose(grx);



  


  if (best.size()>1){
    printf("best (%f):\n",bestv);
    for (int i=0;i<best.size();i++)
      printf("%d %f %d (%d %d %d), ",best[i],w[best[i]],prev[best[i]],gchar3[best[i]].x,gchar3[best[i]].y,gchar3[best[i]].z);
    printf("\n");
  }

  if (best.size()>1){
      vector<int> temp;
      vector<float> tempw;

      for (int b=0;b<best.size()-1;b++){
	float w1=0;
	char3 no=gchar3[best[b]];
	for (int ct=0;ct<alist[no.x][no.y][no.z].size();ct++){
	  if (best[b+1]==gnode[alist[no.x][no.y][no.z][ct].x][alist[no.x][no.y][no.z][ct].y][alist[no.x][no.y][no.z][ct].z])
	    w1=alist[no.x][no.y][no.z][ct].w;
	}
	//w=1;
	w1=w[best[b]];
	tempw.push_back(pow(w1,0.9f));
	temp.push_back(best[b]);
	if (b==best.size()-2)
	  temp.push_back(best[b+1]);
      }

      //normalizza w
      float tw=0;
      for (int x=0;x<tempw.size();x++)
	tw+=tempw[x];
      for (int x=0;x<tempw.size();x++)
	tempw[x]/=tw;

      for (int x=0;x<temp.size();x++)
	printf("%d ",temp[x]);
      printf("\n");
      for (int x=0;x<tempw.size();x++)
	printf("%f ",tempw[x]);
      printf("\n");

      if (stro_normal==0)
	stro_normal=sqrt(bestv);
      bestv=0.3;
      int passo=(int)(400*sqrt(bestv)/stro_normal);
      int side=(passo*stro1.height())/stro1.width();

      stro2.resize(passo,(passo*stro1.height())/stro1.width(),-100,-100,5);

      printf("%d x %d\n",passo,side);
	float ct=0;
	//printf("draw %d %d, %f, %d\n",temp.size(),tempw.size(),bestv,base);
	//printf("%d %d\n",side,passo);

	for (int px=base;px<base+side;px++)
	  if (px<1024)
	    for (int k=0;k<passo && k<1024;k++) {
	      //printf("%d %d\n",px,k);
	      float alpha=stro2(k,px-base,0,0)/255.0;
	      float grigio=stro2(k,px-base,0,1);
	      alpha=1;
	      grigio=0;

	      ct=k*1.0/passo;

	      ct+=(alpha-0.5)*(grigio/255.0-0.5)*0.2;

	      //anticipo un po'
	      ct=ct*1.3;

	      if (ct<0) ct=0;
	      if (ct>0.999) ct=0.999;
	      //printf("k %d ct %f del %f\n",k,ct,(grigio/255.0-0.5)*0.5);
	      // cerco indice
	      int r1=0;
	      float r1f=tempw[r1++];
	      float r1fl=0;//last
	      //printf("%d %f %f\n",r1,r1f,r1fl);
	      while (r1f<ct) {r1fl=r1f;r1f+=tempw[r1++];}
	  
	      //printf("%f: %f %f %d (%d)\n",ct,r1fl,r1f,r1,temp.size());
	      //stat1(px,pos,0,0)=(1-se)*temp[k].x*S+se*temp[k+1].x*S;
	      float se=(ct-r1fl)/(r1f-r1fl);
	      float r,g,b;
	      r=((1-se)*gchar3[temp[r1-1]].x*S+se*gchar3[temp[r1]].x*S);
	      g=((1-se)*gchar3[temp[r1-1]].y*S+se*gchar3[temp[r1]].y*S);
	      b=((1-se)*gchar3[temp[r1-1]].z*S+se*gchar3[temp[r1]].z*S);
	      //printf("ok\n");

	      float r2,g1,b1;
	      r2=alpha*r+(1-alpha)*grigio; // grigio
	      g1=alpha*g+(1-alpha)*grigio; // grigio
	      b1=alpha*b+(1-alpha)*grigio; // grigio

	      /*	      stat11(k,px,0,0)=(int)r2;
	      stat11(k,px,0,1)=(int)g1;
	      stat11(k,px,0,2)=(int)b1;
	      */
	      //printf("%d %d %f\n",px,k,r);

	    }
	base+=side;
	printf("oook\n");
	int isolo=0;
	for (int i=2;i<best.size()-2;i++){
	  if (isolo){
	    char3 no=gchar3[best[i]];
	    for (int x=-1;x<=1;x++)
	    for (int y=-1;y<=1;y++)
	    for (int z=-1;z<=1;z++)
	      //	      if (x==0 && y==0 && z==0)
	      if (no.x+x>=0 && no.x+x<NB &&
		  no.y+y>=0 && no.y+y<NB &&
		  no.z+z>=0 && no.z+z<NB)
	    for (int ct=0;ct<alist[no.x+x][no.y+y][no.z+z].size();ct++){
	      alist[no.x+x][no.y+y][no.z+z][ct].w=0.01;
	    }
	  }
	  else{
	    char3 no=gchar3[best[i]];
	    for (int ct=0;ct<alist[no.x][no.y][no.z].size();ct++){
	      if (best[i+1]==gnode[alist[no.x][no.y][no.z][ct].x][alist[no.x][no.y][no.z][ct].y][alist[no.x][no.y][no.z][ct].z])
		alist[no.x][no.y][no.z][ct].w=0;
	    }

	  }
	}
  }

  
  }

    // via corti o peso 0 o nan
    for (int j=lists.size()-1;j>=0;j--){
      if (listsw[j]!=listsw[j] || listsw[j]==0 || lists[j].size()<=2) {
	lists.erase(lists.begin()+j);
	listsw.erase(listsw.begin()+j);
      }
    }
      /*
    // valuto path deviazione linea
// NO: fatto meglio prima
    for (int j=lists.size()-1;j>=0;j--){
      
      printf("%d: ",j);
      float w=0;
      
      int n1,n2;
      n1=lists[j][0];
      n2=lists[j][lists[j].size()-1];
      char3 no1=gchar3[n1];
      char3 no2=gchar3[n2];
      
      float vec1[3];
      vec1[0]=no2.x-no1.x;
      vec1[1]=no2.y-no1.y;
      vec1[2]=no2.z-no1.z;
      float norm=pow(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2],0.5f);
      vec1[0]/=norm;vec1[1]/=norm;vec1[2]/=norm;      
      for (int i=1;i<lists[j].size();i++){
	char3 no3=gchar3[lists[j][i]];
	float vec2[3];
	vec2[0]=no3.x-no1.x;
	vec2[1]=no3.y-no1.y;
	vec2[2]=no3.z-no1.z;
	float d=vec2[0]*vec1[0]+vec2[1]*vec1[1]+vec2[2]*vec1[2];
	float norm2=vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2];
	float proj=pow(norm2-d*d,0.5f);
	printf("%d (%3.3f)- ",lists[j][i],proj/norm);
	if (w<proj/norm) w=proj/norm;
      }
      //printf("| %f\n",w);
      if (w>=0.3 || lists[j].size()<=1) {
	//printf("kill\n");
	lists.erase(lists.begin()+j);
      }
    }
*/
      
 printf("16-"); 
  // rescoring
  for (int j=0;j<lists.size();j++){
    float we=0;
    for (int i=0;i<lists[j].size()-1;i++){
      int a=lists[j][i];
      int b=lists[j][i+1];
    
      //cerco arco corretto
      for (int ct=0;ct<adj[a].size();ct++){
	int i1=adj[a][ct].y;
	float v=adj[a][ct].w;
	if (b==i1)
	  we+=v;
      }
    }
    listsw[j]=we;    
  }

 printf("17-"); 
  /*  for (int j=0;j<lists.size();j++){
    float we=0;
    for (int i=0;i<lists[j].size();i++)
      we+=w[lists[j][i]];
    char3 no1=gchar3[lists[j][0]];
    char3 no2=gchar3[lists[j][lists[j].size()-1]];
    
    float dist=pow(
		   pow(no1.x-no2.x,2.0)+
		   pow(no1.y-no2.y,2.0)+
		   pow(no1.z-no2.z,2.0)
		   ,0.5f);
    we*=dist;
    listsw[j]=we;
    
  }
*/    
 printf("18-"); 
  /// sort n^2
  for (int j1=0;j1<lists.size();j1++) {
    for (int j2=j1+1;j2<lists.size();j2++){
      if (listsw[j1]<listsw[j2]){
	// swap
	float t=listsw[j1];listsw[j1]=listsw[j2];listsw[j2]=t;
	lists[j1].swap(lists[j2]);
      }
    }
  }
   printf("19-"); 

  /// clean liste

    // ragiono solo con nodes=2
    
  int mod=1;
  if (1)
  while(mod) {
    mod=0;
    vector<int> l1;
    vector<int> l2;
 printf("1------------\n");
 // for (int j=0;j<lists.size();j++){   }
    

    //togli doppioni
    for (int j=0;j<lists.size();j++) {

      printf("%d, %d %f: ",j,lists[j].size(),listsw[j]);
      for (int i=0;i<lists[j].size();i++)
  	if (nodes[lists[j][i]]==2)
  	printf("%d - ",lists[j][i]);
      printf("\n");
      
      l1.clear();
      for (int i1=0; i1<lists[j].size();i1++)
  	if (nodes[lists[j][i1]]==2)
  	  l1.push_back(lists[j][i1]);

      int cut=-1; // nodo su cui spezzare
      /// analisi retta l1
      if (l1.size()>2) {
  	float worst=1;
  	for (int i0=0;i0<2 && i0<l1.size();i0++) {
  	  for (int i2=i0+1;i2<l1.size();i2++) {
  	    for (int i1=i2+1;i1<l1.size();i1++)
  	      if (i1>=l1.size()-2)
  	      {
  	      char3 no1=gchar3[l1[i0]];
  	      char3 no3=gchar3[l1[i1]];
  	      char3 no2=gchar3[l1[i2]];
  	      float v1[3];
  	      float v2[3];
  	      v1[0]=no1.x-no2.x;
  	      v1[1]=no1.y-no2.y;
  	      v1[2]=no1.z-no2.z;
  	      v2[0]=no3.x-no2.x;
  	      v2[1]=no3.y-no2.y;
  	      v2[2]=no3.z-no2.z;
  	      normal(v1);
  	      normal(v2);
  	      float ps=-(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
  	      if (worst>ps && ps<0.95){
  		worst=ps;
  		cut=l1[i2];
  	      }
  	      printf("%d: %d %f (%d %d %d)\n",j,i2,ps,l1[i0],l1[i2],l1[i1]);
  	    }
  	  }
  	}
      }

         cut=-1;
      if (cut>=0) {
  	printf("split at node %d\n",cut);
  	lists.resize(lists.size()+1);
  	listsw.resize(listsw.size()+1);
  	int idx=lists.size()-1;
  	/*	printf("pre\n");
  	for (int i=0;i<lists[idx].size();i++){
  	  printf("%d (%d)--> ",lists[idx][i],i);
  	}
  	*/

  	int copy=0;
  	for (int i=0;i<lists[j].size();i++){
  	  if (lists[j][i]==cut || copy){
  	    copy=1;
  	    printf("copy %d (%d) ",i,lists[j][i]);
  	    lists[idx].push_back(lists[j][i]);
  	  }
  	  else
  	    printf("skip %d (%d) ",i,lists[j][i]);
  	}
  	printf("\n");
  	for (int i=0;i<lists[idx].size();i++){
  	  printf("%d (%d)--> ",lists[idx][i],i);
  	}
  	printf("\nkeep: ");
	
  	for (int i=0;i<lists[j].size();i++){
  	  printf("%d (%d)> ", lists[j][i],i);
  	  if (lists[j][i]==cut){
  	    printf("cut\n");
  	    lists[j].erase(lists[j].begin()+i+1,lists[j].end());
  	  }
  	}

 	//rescoring
 	for (int j1=0;j1<2;j1++){
 	  int j2;
 	  if (j1==0) j2=j; else j2=idx;
 	  float we=0;
 	  for (int i=0;i<lists[j].size()-1;i++){
 	    int a=lists[j][i];
 	    int b=lists[j][i+1];
    
 	    //cerco arco corretto
 	    for (int ct=0;ct<adj[a].size();ct++){
 	      int i1=adj[a][ct].y;
 	      float v=adj[a][ct].w;
 	      if (b==i1)
 		we+=v;
 	    }
 	  }
 	  listsw[j]=we;
 	}

  
	
      }

      
      for (int j1=0;j1<lists.size();j1++)
  	if (j!=j1) {

  	l2.clear();
  	for (int i1=0; i1<lists[j1].size();i1++)
  	  if (nodes[lists[j1][i1]]==2)
  	    l2.push_back(lists[j1][i1]);

  	// l1 completamente contenuta in l2
  	if (l1.size()>0){
  	  int match=0;
  	  for (int i1=0; i1<l1.size();i1++){
  	    for (int i2=0; i2<l2.size();i2++){
  	      if (l1[i1]==l2[i2]){
  		match++;
  		i2=l2.size();
  	      }
  	    }
  	  }
  	  if (match>0 && match==l1.size()
  	      && ((match==l2.size() && listsw[j]<listsw[j1])||
  		  (match<l2.size()))){ // l1 contenuta o uguale l2
  	    lists[j].clear();
  	    printf("rimosso %d --> %d\n",j1,j);
  	  }
  	}

  	// l2 completamente contenuta in l1
  	if (l2.size()>0){
  	  int match=0;
  	  for (int i1=0; i1<l2.size();i1++){
  	    for (int i2=0; i2<l1.size();i2++){
  	      if (l2[i1]==l1[i2]){
  		match++;
  		i2=l1.size();
  	      }
  	    }
  	  }
  	  if (match>0 && match==l2.size() &&
  	      ((match==l1.size() && listsw[j1]<listsw[j])||
  	       (match<l1.size()))){ // l2 contiene o uguale in l1
  	    lists[j1].clear();
  	    printf("rimosso %d --> %d\n",j, j1);
  	  }
  	}
	
      }
    }
    
    for (int j=lists.size()-1;j>=0;j--){
      if (lists[j].size()==0){        /// elimino anche se peso troppo alto
  	lists.erase(lists.begin()+j);
  	listsw.erase(listsw.begin()+j);
  	mod=1;
      }
    }
  }















  /*
  // stima nodi interni
  vector<int> nodesint;
  vector<int> nodesbor;
  nodesint.resize(nnode);
  nodesbor.resize(nnode);
  for (int i=0;i<nnode;i++){
    nodesint[i]=0;
    nodesbor[i]=0;
  }
  
    for (int j=0;j<lists.size();j++)
      if (listsw[j]>0)
      for (int b=0;b<lists[j].size();b++){
	if (nodes[lists[j][b]]==2){
	  if (b==0) // non conto estremo destro perche' deriva da sp che parte da quel nodo!
	    nodesbor[lists[j][b]]++;
	  else
	    nodesint[lists[j][b]]++;		
	}	
      }

    for (int i=0;i<nnode;i++)
      if (nodes[i]==2){
	printf("%d: b %d i %d",i,nodesbor[i],nodesint[i]);
	if (nodesbor[i]>0){
	  nodes[i]=3; ///// unici che tengo
	  printf("*");
	}
	printf("\n");
      }
*/	     


  vector< vector<int> > retta;
  retta.resize(lists.size());
  for (int i=0;i<lists.size();i++)
    retta[i].resize(6);
  vector< vector<int> > bezier;
  bezier.resize(lists.size());
  for (int i=0;i<lists.size();i++)
    bezier[i].resize(12);
  /*calcolo regressioni*/
 /*  //lineare (con ellissoide) */
/*   for (int j1=0;j1<lists.size();j1++){ */

/*   double l1[3]; // vettore di regressione */
  
/*   /// restistuisce i vettori base e baricentro */
/* //   float c[3]; // baricentro ellissoide */

/* 	ell_v=(float**)malloc(4*sizeof(float*)); */
/* 	ell_v[0]=(float*)malloc(4*sizeof(float)); */
/* 	ell_v[1]=(float*)malloc(4*sizeof(float)); */
/* 	ell_v[2]=(float*)malloc(4*sizeof(float)); */
/* 	ell_v[3]=(float*)malloc(4*sizeof(float)); */

/* 	int nrot=10; */
	
/* 	float mass=0; */
/* 	ell_c[0]=0; */
/* 	ell_c[1]=0; */
/* 	ell_c[2]=0; */
	
/* 	for (int j2=0;j2<lists[j1].size();j2++){ */
/* 	  char3 no=gchar3[lists[j1][j2]]; */
/* 	  float we=w[lists[j1][j2]]; */
/* 	  ell_c[0]+=we*no.x;		 */
/* 	  ell_c[1]+=we*no.y;		 */
/* 	  ell_c[2]+=we*no.z;		 */
/* 	  mass+=we;		 */
/* 	} */
/* 	ell_c[0]=ell_c[0]/mass; */
/* 	ell_c[1]=ell_c[1]/mass; */
/* 	ell_c[2]=ell_c[2]/mass; */

/* 	float** ell_a; */
/* 	ell_a=(float**)malloc(4*sizeof(float*)); */
/* 	ell_a[0]=(float*)malloc(4*sizeof(float)); */
/* 	ell_a[1]=(float*)malloc(4*sizeof(float)); */
/* 	ell_a[2]=(float*)malloc(4*sizeof(float)); */
/* 	ell_a[3]=(float*)malloc(4*sizeof(float)); */
/* 	for (int idx1=1;idx1<=3;idx1++) */
/* 	for (int idx2=1;idx2<=3;idx2++) */
/* 		ell_a[idx1][idx2]=0; */
	
/* 	for (int j2=0;j2<lists[j1].size();j2++){ */
/* 	  char3 no=gchar3[lists[j1][j2]];	   */
   
/* 	      double t[3]; */
/* 	      t[0]=no.x; */
/* 	      t[1]=no.y; */
/* 	      t[2]=no.z; */
	    
/* 	      for (int idx1=1;idx1<=3;idx1++) */
/* 		for (int idx2=1;idx2<=3;idx2++) */
/* 		  { */
/* 		    double temp=0; */
/* 		    double de=1; */
/* 		    temp=w[lists[j1][j2]]*(t[idx1-1]-ell_c[idx1-1])* */
/* 		      (t[idx2-1]-ell_c[idx2-1]);		 */
/* 		    ell_a[idx1][idx2]+=temp/mass; */
/* 		  } */
/* 	    } */
/* 	//printf("Launch jacobi\n"); */
/*    jacobi(ell_a, 3, ell_d, ell_v, &nrot); */

/*    int min,max,med; */
/* 	min=1; */
/* 	if (ell_d[min]>ell_d[2]) min=2; */
/* 	if (ell_d[min]>ell_d[3]) min=3; */
/* 	max=1; */
/* 	if (ell_d[max]<ell_d[2]) max=2; */
/* 	if (ell_d[max]<ell_d[3]) max=3; */
/* 	if (min!=1 && max!=1) med=1; */
/* 	if (min!=2 && max!=2) med=2; */
/* 	if (min!=3 && max!=3) med=3;	 */

/* 	/\*   normal[0]=ell_v[1][min]; */
/*    normal[1]=ell_v[2][min]; */
/*    normal[2]=ell_v[3][min]; */

/*    l2[0]=ell_v[1][med]; */
/*    l2[1]=ell_v[2][med]; */
/*    l2[2]=ell_v[3][med]; */

/*    l1[0]=ell_v[1][max]; */
/*    l1[1]=ell_v[2][max]; */
/*    l1[2]=ell_v[3][max]; */
/* 	*\/ */

/*    // calcola intersezioni con 6 piani */
/*    float lam[6]; */
/*    lam[0]=-ell_c[0]/l1[0]; */
/*    lam[1]=-ell_c[1]/l1[1]; */
/*    lam[2]=-ell_c[2]/l1[2]; */
/*    lam[3]=(NB-ell_c[0])/l1[0]; */
/*    lam[4]=(NB-ell_c[1])/l1[1]; */
/*    lam[5]=(NB-ell_c[2])/l1[2]; */
/*    float lm=0,lM=0; */
/*    // scelgo spostamenti minimi da baricentro nei due versi */
/*    for (int i=0;i<6;i++) */
/*      if (lam[i]<0){ */
/*        int ok=1; */
/*        for (int j=0;j<6;j++) */
/* 	 if (lam[j]<0 && lam[j]>lam[i]) */
/* 	   ok=0; */
/*        if (ok) lm=lam[i]; */
/*      } */
/*    for (int i=0;i<6;i++) */
/*      if (lam[i]>0){ */
/*        int ok=1; */
/*        for (int j=0;j<6;j++) */
/* 	 if (lam[j]>0 && lam[j]<lam[i]) */
/* 	   ok=0; */
/*        if (ok) lM=lam[i]; */
/*      }  */

/*    retta[j1][0]=(int)round(ell_c[0]+lm*l1[0]); */
/*    retta[j1][1]=(int)round(ell_c[1]+lm*l1[1]); */
/*    retta[j1][2]=(int)round(ell_c[2]+lm*l1[2]); */
/*    retta[j1][3]=(int)round(ell_c[0]+lM*l1[0]); */
/*    retta[j1][4]=(int)round(ell_c[1]+lM*l1[1]); */
/*    retta[j1][5]=(int)round(ell_c[2]+lM*l1[2]); */

/*    if (ell_d[max]/ell_d[med]<10){ /// comincio a sentire l'ellissoide == non e' molto retta */
/*      printf("dovrei cancellare %d %f %f %f\n",j1,ell_d[max],ell_d[med],ell_d[min]); */
/*      //lists[j1].clear(); */
/*        } */
   
/*    if (1) */
/*    printf("%2d: --> %2d %2d %2d - %2d %2d %2d bar %3.3f %3.3f %3.3f, str Max %2.3f Med %2.3f Min %2.3f  line %1.3f %1.3f %1.3f, (%2.3f %2.3f %2.3f, %2.3f %2.3f %2.3f) \n",j1, */
/* 	  (int)round(ell_c[0]+lm*l1[0]), */
/* 	  (int)round(ell_c[1]+lm*l1[1]), */
/* 	  (int)round(ell_c[2]+lm*l1[2]), */
/* 	  (int)round(ell_c[0]+lM*l1[0]), */
/* 	  (int)round(ell_c[1]+lM*l1[1]), */
/* 	  (int)round(ell_c[2]+lM*l1[2]), */
/* 	  ell_c[0],ell_c[1],ell_c[2], */
/* 	  ell_d[max], */
/* 	  ell_d[med], */
/* 	  ell_d[min], */
/* 	  l1[0],l1[1],l1[2], */
/* 	  lam[0],lam[1],lam[2],lam[3],lam[4],lam[5],lam[6] */
/* 	  ); */
   
/*   } */


/*   /// clustering rette con priorita' ai primi (piu' forti) */
/*   if (0) */
/*   for (int j1=0;j1<lists.size();j1++) */
/*     if (lists[j1].size()>0){ // non cancellato */
/*       for (int j2=j1+1;j2<lists.size();j2++) */
/* 	if (lists[j2].size()>0){ // non cancellato */
/* 	  if (abs(retta[j1][0]-retta[j2][0])<=2 && */
/* 	      abs(retta[j1][1]-retta[j2][1])<=2 && */
/* 	      abs(retta[j1][2]-retta[j2][2])<=2 && */
/* 	      abs(retta[j1][3]-retta[j2][3])<=2 && */
/* 	      abs(retta[j1][4]-retta[j2][4])<=2 && */
/* 	      abs(retta[j1][5]-retta[j2][5])<=2){ */
/* 	    lists[j2].clear(); */
/* 	    printf("del %d %d\n",j1,j2); */
/* 	  } */
/* 	} */
/*     } */

  
/*   // cancella */
/*   for (int j=lists.size()-1;j>=0;j--){ */
/*     if (lists[j].size()==0){ */
/*       lists.erase(lists.begin()+j); */
/*       listsw.erase(listsw.begin()+j); */
/*       retta.erase(retta.begin()+j); */
/*       bezier.erase(bezier.begin()+j); */
/*     } */
/*   } */

      
  for (int j1=0;j1<lists.size();j1++){
	/// bezier
	float A1=0,A2=0,A12=0;
	float B0,B1,B2,B3;
	float C1[3];
	float C2[3];
	float P0[3];
	float P1[3];
	float P2[3];
	float P3[3];
	C1[0]=0;
	C1[1]=0;
	C1[2]=0;
	C2[0]=0;
	C2[1]=0;
	C2[2]=0;
	vector<float> temp;
	vector<float> tempw;
	vector<float> tempw1;
	vector<char3> tempcol;
      
      float tw=0;      
      	for (int b=0;b<lists[j1].size();b++){
	float w1=0;
	int nod=lists[j1][b];
	printf("%d ",nod);
	tempw.push_back(pow(w[nod],0.99));
	temp.push_back(lists[j1][b]);
      }
      printf("\n");
      //normalizza w
      for (int x=0;x<tempw.size();x++)
	tw+=tempw[x];
      for (int x=0;x<tempw.size();x++){
	tempw[x]/=tw;
      }
      for (int x=1;x<tempw.size();x++){
	tempw[x]+=tempw[x-1];
	printf("%f ",tempw[x]);
      }
      printf("\n");

      int passo=40;
      for (int i=0;i<passo;i++) {
	float pos=(float)i/passo;
	float pos1=(float)i/passo+1.0/passo;
	// controllo da pos .. pos+1/passo
	float c[3];
	float wc=0;
	c[0]=0;
	c[1]=0;
	c[2]=0;
	for (int j=0;j<tempw.size();j++){
	  if (!(tempw[j]<pos ||              // e' prima
		(j==0 && 0>=pos1) ||	      // primo 
		(j>0 && tempw[j-1]>=pos1) // e' dopo		
		)){
	    float a,b; // estremi colore
	    if (j==0)
	      a=0;
	    else
	      a=tempw[j-1];
	    b=tempw[j];
	    //printf("%d ",j);
	    if ( b>pos && b<pos1 && a<=pos ){ // a cavallo
	      char3 no=gchar3[temp[j]];
	      float w=0;
	      w=b-pos;
	      //printf("add %f %d, %f - %f %f %f\n",w,j,pos,pos1,a,b);
	      c[0]+=w*no.x;
	      c[1]+=w*no.y;
	      c[2]+=w*no.z;
	      wc+=w;
	    }
	    if ( b>=pos1 && a<pos1 && a>pos ){ // a cavallo
	      char3 no=gchar3[temp[j]];
	      float w=0;
	      w=pos1-a;
	      //printf("add1 %f %d, %f - %f %f %f\n",w,j,pos,pos1,a,b);
	      c[0]+=w*no.x;
	      c[1]+=w*no.y;
	      c[2]+=w*no.z;
	      wc+=w;
	    }
	    if ( b<pos1 && a>pos ){ //dentro
	      char3 no=gchar3[temp[j]];
	      float w=0;
	      w=b-a;
	      //printf("add2 %f %d, %f - %f %f %f\n",w,j,pos,pos1,a,b);
	      c[0]+=w*no.x;
	      c[1]+=w*no.y;
	      c[2]+=w*no.z;
	      wc+=w;
	    }	    
	    if ( b>=pos1 && a<=pos ){ //tutto coperto
	      char3 no=gchar3[temp[j]];
	      float w=0;
	      w=pos1-pos;
	      //	      printf("add3 %f %d, %f - %f %f %f\n",w,j,pos,pos1,a,b);
	      c[0]+=w*no.x;
	      c[1]+=w*no.y;
	      c[2]+=w*no.z;
	      wc+=w;
	    }	    
	  }
	}
	c[0]/=wc;
	c[1]/=wc;
	c[2]/=wc;
	char3 temp;
	temp.x=(int)round(c[0]);
	temp.y=(int)round(c[1]);
	temp.z=(int)round(c[2]);
	tempcol.push_back(temp);
      }
	
      int first=0;
      float acc=0;
      //     while (acc+tempw[first]<0.0000001) first++;
      int last=tempw.size();
      acc=0;
      //while (acc+tempw[last]<0.0000001) last--;
      
	char3 no=tempcol[0];
	P0[0]=no.x;
	P0[1]=no.y;
	P0[2]=no.z;
	no=tempcol[tempcol.size()-1];
	P3[0]=no.x;
	P3[1]=no.y;
	P3[2]=no.z;

	for (int i=1;i<tempcol.size()-1;i++){
	  float t=(float)i/(tempcol.size()-1);

	  B0 = pow(1-t,3.0f);//        % Bezeir Basis
	  B1 = pow( 3*t*(1-t),2.0f ) ;
	  B2 = ( 3*t*t*(1-t) ) ;
	  B3 = pow(t,3.0f)                ;
      
	  A1  = A1 +  B1*B1;
	  A2  = A2 +  B2*B2;
	  A12 = A12 + B1*B2;
	  char3 no=tempcol[i];
	  C1[0]+= B1*( no.x - B0*P0[0] - B3*P3[0] );
	  C1[1]+= B1*( no.y - B0*P0[1] - B3*P3[1] );
	  C1[2]+= B1*( no.z - B0*P0[2] - B3*P3[2] );
	  C2[0]+= B2*( no.x - B0*P0[0] - B3*P3[0] );
	  C2[1]+= B2*( no.y - B0*P0[1] - B3*P3[1] );
	  C2[2]+= B2*( no.z - B0*P0[2] - B3*P3[2] );
	}
   
	float den=(A1*A2-A12*A12);//       % common denominator for all points
	if(den==0){
	  P1[0]=P0[0];
	  P1[1]=P0[1];
	  P1[2]=P0[2];
	  P2[0]=P3[0];
	  P2[1]=P3[1];
	  P2[2]=P3[2];
	}
	else{
	  P1[0]=(A2*C1[0]-A12*C2[0])/den;
	  P1[1]=(A2*C1[1]-A12*C2[1])/den;
	  P1[2]=(A2*C1[2]-A12*C2[2])/den;
	  P2[0]=(A1*C2[0]-A12*C1[0])/den;
	  P2[1]=(A1*C2[1]-A12*C1[1])/den;
	  P2[2]=(A1*C2[2]-A12*C1[2])/den;
	}

	bezier[j1][0]=(int)round(P0[0]);
	bezier[j1][1]=(int)round(P0[1]);
	bezier[j1][2]=(int)round(P0[2]);
	bezier[j1][3]=(int)round(P1[0]);
	bezier[j1][4]=(int)round(P1[1]);
	bezier[j1][5]=(int)round(P1[2]);
	bezier[j1][6]=(int)round(P2[0]);
	bezier[j1][7]=(int)round(P2[1]);
	bezier[j1][8]=(int)round(P2[2]);
	bezier[j1][9]=(int)round(P3[0]);
	bezier[j1][10]=(int)round(P3[1]);
	bezier[j1][11]=(int)round(P3[2]);	

	printf("%d: %d %d %d - %d %d %d - %d %d %d - %d %d %d\n",j1,
	       bezier[j1][0],bezier[j1][1],bezier[j1][2],
	       bezier[j1][3],bezier[j1][4],bezier[j1][5],
	       bezier[j1][6],bezier[j1][7],bezier[j1][8],
	       bezier[j1][9],bezier[j1][10],bezier[j1][11]);

	// errori
	for (int i=first;i<last;i++){
	  float t=(float)i/(last-first-1);
	  t=tempw[i];
	  B0 = pow(1-t,3.0f);//        % Bezeir Basis
	  B1 = pow( 3*t*(1-t),2.0f ) ;
	  B2 = ( 3*t*t*(1-t) ) ;
	  B3 = pow(t,3.0f)                ;
      
	  A1  = A1 +  B1*B1;
	  A2  = A2 +  B2*B2;
	  A12 = A12 + B1*B2;
	  char3 no=gchar3[lists[j1][i]];
	  C1[0]+= B1*( no.x - B0*P0[0] - B3*P3[0] );
	  C1[1]+= B1*( no.y - B0*P0[1] - B3*P3[1] );
	  C1[2]+= B1*( no.z - B0*P0[2] - B3*P3[2] );
	  C2[0]+= B2*( no.x - B0*P0[0] - B3*P3[0] );
	  C2[1]+= B2*( no.y - B0*P0[1] - B3*P3[1] );
	  C2[2]+= B2*( no.z - B0*P0[2] - B3*P3[2] );
	}
      }
      
  /// clustering bezier con priorita' ai primi (piu' forti)
  if (1)
  for (int j1=0;j1<lists.size();j1++)
    if (lists[j1].size()>0){ // non cancellato
      for (int j2=j1+1;j2<lists.size();j2++)
	  if (lists[j1].size()>0 && lists[j2].size()>0){ // non cancellato
	    int n=20;
	    int dbg=(j1==5 && j2==61);
	    float hsv1[3*n];
	    float hsv2[3*n];
	    
	    for (int i=0;i<n;i++){
	      float r1,g1,b1,r2,g2,b2;
	      float h1,h2,s1,s2,v1,v2;
	      float t=((float)i)/(n-1);
	      float t1=((float)i)/(n-1);
	      r1=(pow(1-t,3)*bezier[j1][0]+3*t*(1-t)*(1-t)*bezier[j1][3]+3*t*t*(1-t)*bezier[j1][6]+t*t*t*bezier[j1][9]);
	      g1=(pow(1-t,3)*bezier[j1][1]+3*t*(1-t)*(1-t)*bezier[j1][4]+3*t*t*(1-t)*bezier[j1][7]+t*t*t*bezier[j1][10]);
	      b1=(pow(1-t,3)*bezier[j1][2]+3*t*(1-t)*(1-t)*bezier[j1][5]+3*t*t*(1-t)*bezier[j1][8]+t*t*t*bezier[j1][11]);
	      r2=(pow(1-t1,3)*bezier[j2][0]+3*t1*(1-t1)*(1-t1)*bezier[j2][3]+3*t1*t1*(1-t1)*bezier[j2][6]+t1*t1*t1*bezier[j2][9]);
	      g2=(pow(1-t1,3)*bezier[j2][1]+3*t1*(1-t1)*(1-t1)*bezier[j2][4]+3*t1*t1*(1-t1)*bezier[j2][7]+t1*t1*t1*bezier[j2][10]);
	      b2=(pow(1-t1,3)*bezier[j2][2]+3*t1*(1-t1)*(1-t1)*bezier[j2][5]+3*t1*t1*(1-t1)*bezier[j2][8]+t1*t1*t1*bezier[j2][11]);
	      if (r1<0) r1=0;
	      if (g1<0) g1=0;
	      if (b1<0) b1=0;
	      if (r2<0) r2=0;
	      if (g2<0) g2=0;
	      if (b2<0) b2=0;
	      if (r1>=NB) r1=NB-1;
	      if (g1>=NB) g1=NB-1;
	      if (b1>=NB) b1=NB-1;
	      if (r2>=NB) r2=NB-1;
	      if (g2>=NB) g2=NB-1;
	      if (b2>=NB) b2=NB-1;	      
	      RGBtoHSV(r1*S,g1*S,b1*S,&h1,&s1,&v1);
	      RGBtoHSV(r2*S,g2*S,b2*S,&h2,&s2,&v2);
	      hsv1[3*i+0]=h1;
	      hsv1[3*i+1]=s1;
	      hsv1[3*i+2]=v1;
	      hsv2[3*i+0]=h2;
	      hsv2[3*i+1]=s2;
	      hsv2[3*i+2]=v2;
	    }
	    // per ogni punto trova il migliore nell'altro
	    float besta=0;
	    float bestb=0;
	    for (int swap=0;swap<2;swap++)
	      for (int i=0;i<n;i++){
		float bestdist=1000;
		for (int j=0;j<n;j++){
		  float h1,h2,s1,s2,v1,v2;
		  if (swap==0){
		    h1=hsv1[3*i+0];
		    s1=hsv1[3*i+1];
		    v1=hsv1[3*i+2];
		    h2=hsv2[3*j+0];
		    s2=hsv2[3*j+1];
		    v2=hsv2[3*j+2];
		  }else{
		    h1=hsv2[3*i+0];
		    s1=hsv2[3*i+1];
		    v1=hsv2[3*i+2];
		    h2=hsv1[3*j+0];
		    s2=hsv1[3*j+1];
		    v2=hsv1[3*j+2];
		  }
		  //	    float dist=(r1-r2)*(r1-r2)+(g1-g2)*(g1-g2)+(b1-b2)*(b1-b2);
		  float c1x,c1y;
		  float c2x,c2y;
		  c1x=cos(h1*3.1415)*s1;
		  c1y=sin(h1*3.1415)*s1;
		  c2x=cos(h2*3.1415)*s2;
		  c2y=sin(h2*3.1415)*s2;
		  float disth=pow((c1x-c2x)*(c1x-c2x)+(c1y-c2y)*(c1y-c2y),0.5);
		  /*
		    fabs(h1-h2);
		  if (fabs(h1-h2)>0.5)
		    disth=1-disth;
		  */
		  float dist=disth*disth+(v1-v2)*(v1-v2);
		  if (dbg)
		    printf("    %d %d %f: %f %f %f - %f %f %f\n",i,j,dist,h1,s1,v1,h2,s2,v2);

		  if (bestdist>dist) bestdist=dist;
		}
		if (swap==0)
		  besta+=bestdist;
		else
		  bestb+=bestdist;
	      }	    
	    float acc=0;
	    if (besta<bestb)
	      acc=sqrt(besta/n);
	    else
	      acc=sqrt(bestb/n);
	    printf("  dis %d %d %f best %f %f\n",j1,j2,acc,besta,bestb);
	    if (acc<0.1){
	      if (besta<bestb){
	      printf("del %d <-- %d\n",j1,j2);
	      lists[j1].clear();
	      listsw[j2]+=listsw[j1];
	      }
	      if (besta>bestb){
	      printf("del %d --> %d\n",j1,j2);
	      lists[j2].clear();
	      listsw[j1]+=listsw[j2];
	      }
	    }

	}
    }


  
  // cancella
      for (int j=lists.size()-1;j>=0;j--){
      if (lists[j].size()==0){
	lists.erase(lists.begin()+j);
	listsw.erase(listsw.begin()+j);
	retta.erase(retta.begin()+j);
	bezier.erase(bezier.begin()+j);
	
      }
      }


  printf("------------\n");

  // rescoring
  if (0)
  for (int j=0;j<lists.size();j++){
    float we=0;
    for (int i=0;i<lists[j].size()-1;i++){
      int a=lists[j][i];
      int b=lists[j][i+1];
    
      //cerco arco corretto
      for (int ct=0;ct<adj[a].size();ct++){
	int i1=adj[a][ct].y;
	float v=adj[a][ct].w;
	if (b==i1)
	  we+=v;
      }
    }
    listsw[j]=we;    
  }
 

  
  /// sort n^2
  for (int j1=0;j1<lists.size();j1++){
    for (int j2=j1+1;j2<lists.size();j2++){
      if (listsw[j1]<listsw[j2]){
	// swap
	float t=listsw[j1];listsw[j1]=listsw[j2];listsw[j2]=t;
	lists[j1].swap(lists[j2]);
	bezier[j1].swap(bezier[j2]);
      }
    }
  }

  
  base=0;
    for (int j=0;j<lists.size();j++)
      if (listsw[j]>0){
  
      vector<int> temp;
      vector<float> tempw;
      vector<float> tempw1;
      vector<char3> tempcol;

      
      //printf("%d: ",j);
      float tw=0;      
      for (int b=0;b<lists[j].size();b++){
	float w1=0;
	int nod=lists[j][b];
	//printf("%d ",nod);
	tempw.push_back(pow(w[nod],0.99));
	temp.push_back(lists[j][b]);
      }
      //printf("\n");
      //normalizza w
      for (int x=0;x<tempw.size();x++)
	tw+=tempw[x];
      for (int x=0;x<tempw.size();x++){
	tempw[x]/=tw;
      }
      for (int x=1;x<tempw.size();x++){
	tempw[x]+=tempw[x-1];
	//printf("%f ",tempw[x]);
      }
      printf("\n");

      int passo=40;
      for (int i=0;i<passo;i++) {
	float pos=(float)i/passo;
	float pos1=(float)i/passo+1.0/passo;
	// controllo da pos .. pos+1/passo
	//	printf("%f..%f: ",pos,pos1);
	float c[3];
	float wc=0;
	c[0]=0;
	c[1]=0;
	c[2]=0;
	for (int j=0;j<tempw.size();j++){
	  if (!(tempw[j]<pos ||              // e' prima
		(j==0 && 0>=pos1) ||	      // primo 
		(j>0 && tempw[j-1]>=pos1) // e' dopo		
		)){
	    float a,b; // estremi colore
	    if (j==0)
	      a=0;
	    else
	      a=tempw[j-1];
	    b=tempw[j];
	    //printf("%d ",j);
	    if ( b>pos && b<pos1 && a<=pos ){ // a cavallo
	      char3 no=gchar3[temp[j]];
	      float w=0;
	      w=b-pos;
	      //printf("add %f %d, %f - %f %f %f\n",w,j,pos,pos1,a,b);
	      c[0]+=w*no.x;
	      c[1]+=w*no.y;
	      c[2]+=w*no.z;
	      wc+=w;
	    }
	    if ( b>=pos1 && a<pos1 && a>pos ){ // a cavallo
	      char3 no=gchar3[temp[j]];
	      float w=0;
	      w=pos1-a;
	      //printf("add1 %f %d, %f - %f %f %f\n",w,j,pos,pos1,a,b);
	      c[0]+=w*no.x;
	      c[1]+=w*no.y;
	      c[2]+=w*no.z;
	      wc+=w;
	    }
	    if ( b<pos1 && a>pos ){ //dentro
	      char3 no=gchar3[temp[j]];
	      float w=0;
	      w=b-a;
	      //printf("add2 %f %d, %f - %f %f %f\n",w,j,pos,pos1,a,b);
	      c[0]+=w*no.x;
	      c[1]+=w*no.y;
	      c[2]+=w*no.z;
	      wc+=w;
	    }	    
	    if ( b>=pos1 && a<=pos ){ //tutto coperto
	      char3 no=gchar3[temp[j]];
	      float w=0;
	      w=pos1-pos;
	      //	      printf("add3 %f %d, %f - %f %f %f\n",w,j,pos,pos1,a,b);
	      c[0]+=w*no.x;
	      c[1]+=w*no.y;
	      c[2]+=w*no.z;
	      wc+=w;
	    }	    
	  }
	}
	c[0]/=wc;
	c[1]/=wc;
	c[2]/=wc;
	char3 temp;
	temp.x=(int)round(c[0]);
	temp.y=(int)round(c[1]);
	temp.z=(int)round(c[2]);
	tempcol.push_back(temp);
	//printf("%d %d %d\n",temp.x,temp.y,temp.z);
      }

      int side=(passo*stro1.height())/stro1.width();

      stro2.resize(passo,(passo*stro1.height())/stro1.width(),-100,-100,5);

      printf("%d x %d\n",passo,side);
	float ct=0;
	//printf("draw %d %d, %f, %d\n",temp.size(),tempw.size(),bestv,base);
	//printf("%d %d\n",side,passo);

	for (int px=base;px<base+side;px++)
	    for (int k=0;k<passo && k<1024;k++) {
	      //printf("%d %d\n",px,k);
	      float alpha=stro2(k,px-base,0,0)/255.0;
	      float grigio=stro2(k,px-base,0,1);
	      alpha=1;
	      grigio=0;

	      ct=k*1.0/(passo-1);

	      //ct+=(alpha-0.5)*(grigio/255.0-0.5)*0.2;

	      //anticipo un po'
	      //ct=ct*1.3;

	      if (ct<0) ct=0;
	      if (ct>1) ct=1;
	      //printf("k %d ct %f del %f\n",k,ct,(grigio/255.0-0.5)*0.5);
	      float r,g,b;
	      r=tempcol[k].x*S;//((1-se)*gchar3[temp[r1-1]].x*S+se*gchar3[temp[r1]].x*S);
	      g=tempcol[k].y*S;//((1-se)*gchar3[temp[r1-1]].y*S+se*gchar3[temp[r1]].y*S);
	      b=tempcol[k].z*S;//((1-se)*gchar3[temp[r1-1]].z*S+se*gchar3[temp[r1]].z*S);
	      //	      printf("ok\n");

	      float r2,g1,b1;
	      r2=alpha*r+(1-alpha)*grigio; // grigio
	      g1=alpha*g+(1-alpha)*grigio; // grigio
	      b1=alpha*b+(1-alpha)*grigio; // grigio

	      stat11(k+3*passo*(px/1024),px%1024,0,0)=(int)r2;
	      stat11(k+3*passo*(px/1024),px%1024,0,1)=(int)g1;
	      stat11(k+3*passo*(px/1024),px%1024,0,2)=(int)b1;
	      float temp;
	      // bezier
	      float t=(float)k/(passo-1);
	      temp=(pow(1-t,3)*bezier[j][0]+3*t*(1-t)*(1-t)*bezier[j][3]+3*t*t*(1-t)*bezier[j][6]+t*t*t*bezier[j][9])*S;
	      if (temp<0) temp=0;
	      if (temp>255) temp=255;	      
	      stat11(k+passo+3*passo*(px/1024),px%1024,0,0)=(int)temp;//;
	      temp=(pow(1-t,3)*bezier[j][1]+3*t*(1-t)*(1-t)*bezier[j][4]+3*t*t*(1-t)*bezier[j][7]+t*t*t*bezier[j][10])*S;
	      if (temp<0) temp=0;
	      if (temp>255) temp=255;
	      stat11(k+passo+3*passo*(px/1024),px%1024,0,1)=(int)temp;//((((float)k/passo)*(retta[j][1])+(1-(float)k/passo)*(retta[j][4]))*S);
	      temp=(pow(1-t,3)*bezier[j][2]+3*t*(1-t)*(1-t)*bezier[j][5]+3*t*t*(1-t)*bezier[j][8]+t*t*t*bezier[j][11])*S;
	      if (temp<0) temp=0;
	      if (temp>255) temp=255;	      
	      stat11(k+passo+3*passo*(px/1024),px%1024,0,2)=(int)temp;//((((float)k/passo)*(retta[j][2])+(1-(float)k/passo)*(retta[j][5]))*S);

	      // retta
	      /*
	      temp=((((float)k/passo)*(retta[j][0])+(1-(float)k/passo)*(retta[j][3]))*S);
	      if (temp<0) temp=0;
	      if (temp>255) temp=255;	      
	      stat11(k+2*passo+3*passo*(px/1024),px%1024,0,0)=(int)temp;//;
	      temp=((((float)k/passo)*(retta[j][1])+(1-(float)k/passo)*(retta[j][4]))*S);
	      if (temp<0) temp=0;
	      if (temp>255) temp=255;
	      stat11(k+2*passo+3*passo*(px/1024),px%1024,0,1)=(int)temp;//((((float)k/passo)*(retta[j][1])+(1-(float)k/passo)*(retta[j][4]))*S);
	      temp=((((float)k/passo)*(retta[j][2])+(1-(float)k/passo)*(retta[j][5]))*S);
	      if (temp<0) temp=0;
	      if (temp>255) temp=255;	      
	      stat11(k+2*passo+3*passo*(px/1024),px%1024,0,2)=(int)temp;//((((float)k/passo)*(retta[j][2])+(1-(float)k/passo)*(retta[j][5]))*S);
*/
	    }

	base+=side;
    }






    /// ripeso i nodi per cluster su codici

    vector<float> neww;
    neww.resize(nnode);
    for (int i=0;i<nnode;i++)
      neww[i]=0;
      
    for (int i=0;i<nnode;i++)
      if (nodes[i]==2){
	char3 no1=gchar3[i];
	for (int j=0;j<nnode;j++)
	  if (nodes[j]>0){
	    char3 no2=gchar3[j];
	    if (C[no1.x][no1.y][no1.z]==C[no2.x][no2.y][no2.z])
	      neww[i]+=w[j];
	}
      }

    float maxw=0;
    for (int i=0;i<nnode;i++){
      if (nodes[i]==2)
	printf("%d: %f\n",i,neww[i]);
      if (maxw<neww[i])
	maxw=neww[i];
    }
    printf("max %f\n",maxw);

  edges.clear();
  for (int j=0;j<lists.size();j++){
    for (int i=0;i<lists[j].size()-1;i++)


      if (nodes[lists[j][i]]==2){ // c'e' dopo?
	int found=-1;
	for (int j2=i+1;j2<lists[j].size();j2++)
	  if (nodes[lists[j][j2]]==2){
	    found=j2;
	    j2=lists[j].size();
	  }

	// coppia j1 - j2
	if (found>=0){
	  int a,b;
	  a=lists[j][i];
	  b=lists[j][found];
	  if (a>b) {
	    b=lists[j][i];
	    a=lists[j][found];
	  }

	  found=-1; // cerco negli archi
	  for (int k=0;k<edges.size();k++)
	    if (edges[k].x==a && edges[k].y==b)
	      found=k;
	  if (found<0){
	    edge t;
	    t.x=a;
	    t.y=b;
	    t.w=listsw[j];
	    edges.push_back(t);
	  }
	  else{
	    edges[found].w+=listsw[j]; // incremento peso arco
	  }

      
	}
      }
  }
  
  /*
 // tolgo eventuale nodo con 2 archi ( bypassato )
  for (int i=0;i<nnode;i++)
    if (nodes[i]==2) { // per ogni nodo
      printf("node %d\n",i);
      // valuta vicini
      vector<int> vicini;
      for (int j=0;j<edges.size();j++){
	if (edges[j].x==i || edges[j].y==i){
	  vicini.push_back(j);
	}
      }

      // matrice colineari
      for (int j1=0;j1<vicini.size();j1++){
	char3 x1=gchar3[edges[vicini[j1]].x];
	char3 x2=gchar3[edges[vicini[j1]].y];	      
	float v1[3];
	v1[0]=x1.x-x2.x;
	v1[1]=x1.y-x2.y;
	v1[2]=x1.z-x2.z;
	normal(v1);
	for (int j2=j1+1;j2<vicini.size();j2++){
	  printf("%d-%d %d-%d: ",edges[vicini[j1]].x,edges[vicini[j1]].y,
		 edges[vicini[j2]].x,edges[vicini[j2]].y);
	char3 y1=gchar3[edges[vicini[j2]].x];
	char3 y2=gchar3[edges[vicini[j2]].y];
	float v2[3];
	v2[0]=y1.x-y2.x;
	v2[1]=y1.y-y2.y;
	v2[2]=y1.z-y2.z;
	normal(v2);

	if (edges[vicini[j1]].x==edges[vicini[j2]].x ||
	    edges[vicini[j1]].y==edges[vicini[j2]].y){
	  // scambio
	  v2[0]*=-1;
	  v2[1]*=-1;
	  v2[2]*=-1;
	}
	    
	
	printf(" %f %f, %f %f %f - %f %f %f\n",
	       v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2],
	       1-(1-fabs(v1[0]-v2[0]))*(1-fabs(v1[1]-v2[1]))*(1-fabs(v1[2]-v2[2])),
	       v1[0],v1[1],v1[2],
	       v2[0],v2[1],v2[2]);	  
	}
	
      }
      
      
    }
  */
    rangemaxg=0;
  for (int i=0;i<edges.size();i++){
    if (edges[i].w>rangemaxg)
      rangemaxg=edges[i].w;
  }


   //outgraph
    grx=fopen("graphSP.dot","w+");
  fprintf(grx,"digraph{\nnode [shape = circle, fixedsize=true, style = filled, fillcolor = palegreen];\n");
  for (int i=0;i<nnode;i++)
    if (nodes[i]==2)
    {
    char3 no=gchar3[i];
    float w1=neww[i]/maxw;
    if (w1<0) w1=0;
    if (w1>1) w1=1;
    w1=pow(w1,0.5f);
      if (nodes[i]==2)
    fprintf(grx,"%d [shape = doublecircle width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    i,
	    pow(w1,0.5),
	    no.x*S,no.y*S,no.z*S);
    else
    fprintf(grx,"%d [width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    i,
	    pow(w1,0.5),
	    no.x*S,no.y*S,no.z*S);
  }


  
  for (int j=0;j<edges.size();j++){
    float w=1-edges[j].w/rangemaxg;
    //      printf("%d: %f --> %f\n",i,edges[found].w,w);
    if (w<0.0) w=0.0;
    if (w>1) w=1;
    w=pow(w,0.9);
    //w=0;
    fprintf(grx,"%d -> %d [dir=none penwidth=%1.2f color = \"",edges[j].y,edges[j].x,5*(1-w)+0.1);
    fprintf(grx,"#%02x%02x%02x len=%f weight=%f",0,0,0,w,10*(1-w));
    fprintf(grx,"\"];\n");
  }
  fprintf(grx,"}\n");
  fclose(grx);



  FILE *filebez;
  char immagine[NOMEMAX];
  strcpy(immagine,argv[1]);
  strcat(immagine,"-bez.txt");
  filebez=fopen(immagine,"w");
  fprintf(filebez,"%d \n",lists.size());
  for (int j=0;j<lists.size();j++){
    fprintf(filebez,"%f %d %d %d %d %d %d %d %d %d %d %d %d ",listsw[j],
	   bezier[j][0],
	   bezier[j][1],
	   bezier[j][2],
	   bezier[j][3],
	   bezier[j][4],
	   bezier[j][5],
	   bezier[j][6],
	   bezier[j][7],
	   bezier[j][8],
	   bezier[j][9],
	   bezier[j][10],
	   bezier[j][11]
	   );
    /* for (int i=0;i<lists[j].size();i++) { */
    /*   if (nodes[lists[j][i]]==2) */
    /*   fprintf(filebez, "%d -> ",lists[j][i]); */
    /* } */
  fprintf(filebez, "\n");
  }
  fclose(filebez);

  for (int j=0;j<lists.size();j++)
    for (int j1=0;j1<lists[j].size();j1++)
      nodes[lists[j][j1]]=3;

    grx=fopen("graphx2.dot","w+");
  fprintf(grx,"digraph{\nnode [shape = circle, fixedsize=true, style = filled, fillcolor = palegreen];\n");

  /*
  for (int i=0;i<nnode;i++)
    if (nodes[i]==3)
    {
    char3 no=gchar3[i];
    float w1=w[i]/rangemaxgn;
    if (w1<0) w1=0;
    if (w1>1) w1=1;
    
    if (nodes[i]==3)
      fprintf(grx,"%d [shape = doublecircle ",i);
    else{
      if (nodes[i]==2)
	fprintf(grx,"%d [shape = doublecircle ",i);
      else
	fprintf(grx,"%d [",i);
    }
    int v=(int)(255*nneigh[i]);
    
    fprintf(grx,"width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    pow(w1,0.2),
	    //	    	    v,v,v);
	    no.x*S,no.y*S,no.z*S);
    }
  */
  for (int j=0;j<lists.size();j++)
    for (int j1=0;j1<lists[j].size();j1++){
      int i=lists[j][j1];
    char3 no=gchar3[i];
    float w1=w[i]/rangemaxgn;
    if (w1<0) w1=0;
    if (w1>1) w1=1;
    
    fprintf(grx,"\"%d-%d-%d\" [shape = doublecircle ",j,j1,i);
    int v=(int)(255*nneigh[i]);
    
    fprintf(grx,"width=\"%2.4f\" fcolor= \"black\" fillcolor = \"#%02x%02x%02x\"];\n",
	    pow(w1,0.2),
	    //	    	    v,v,v);
	    no.x*S,no.y*S,no.z*S);
    }


  
  for (int j=0;j<lists.size();j++){
    double colw=(listsw[j]+0.0)/listsw[0]*25;
    int r,g,b1;
    palettei(j,r,g,b1);
    
    for (int j1=0;j1<lists[j].size();j1++){
      if (1|| nodes[lists[j][j1]]==2)
	{ // c'e' dopo?
	int found=-1;
	for (int j2=j1+1;j2<lists[j].size();j2++)
	  if (1|| nodes[lists[j][j2]]==2){
	    found=j2;
	    j2=lists[j].size();
	  }

	// coppia j1 - j2
	if (found>=0){
	  int a,b;
	  a=lists[j][j1];
	  b=lists[j][found];
	if (0)
	  if (a>b) {
	    b=lists[j][j1];
	    a=lists[j][found];
	  }
	  fprintf(grx,"\"%d-%d-%d\" -> \"%d-%d-%d\" [dir=none penwidth=%f color = \"",j,j1,a,j,found,b,colw);

	  char3 no1=gchar3[a];
	  
	  fprintf(grx,"#%02x%02x%02x",
		  (int)(no1.x*S),
		  (int)(no1.y*S),
		  (int)(no1.z*S));
	  
	  fprintf(grx,"\" len=\"%f weight=%f",1,10);

	  //	  fprintf(grx,"#%02x%02x%02x len=%f weight=%f",r,g,b1,1,10);
	  
	  fprintf(grx,"\"];\n");
	}
      }
    }
  }

    
  fprintf(grx,"}\n");
  fclose(grx);


  char nomebez[NOMEMAX];
  strcpy(nomebez,argv[1]);
  strcat(nomebez,"-stats11.jpg");

    stat11.save(nomebez);










    
    exit(0);





    //long maxv=0;
  int ctsm=0;
  code=10000;
  //while(ctsm++<10 && code>100)
    {

      // smoothing: 
      if (0)
	for (int ct=0;ct<1;ct++) {
	  for (int x=0;x<NB;x++)
	    for (int y=0;y<NB;y++)
	      for (int z=0;z<NB;z++)
		B[x][y][z]=A[x][y][z];

	  for (int x=0;x<NB;x++)
	    for (int y=0;y<NB;y++)
	      for (int z=0;z<NB;z++){
		int max=1;
		long acc=0;
		int c=0;
		for (int dx=-1;dx<=1;dx++)
		  for (int dy=-1;dy<=1;dy++)
		    for (int dz=-1;dz<=1;dz++)
		      if (x+dx>=0 && y+dy>=0 && z+dz>=0 &&
			  x+dx<NB && y+dy<NB && z+dz<NB){
			int w=1;
			if (dx==0&&dy==0&&dz==0)
			  w=50;
			acc+=w*B[x+dx][y+dy][z+dz];
			c+=w;
		      }
		A[x][y][z]=acc/c;
	      }
	}

  //remap
  maxv=0;
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++){
	C[x][y][z]=0;
	B[x][y][z]=A[x][y][z];
	if (maxv<A[x][y][z]) maxv=A[x][y][z];
      }

  code=0;
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++)
	if (B[x][y][z]>0) {	  
	  rec(x,y,z);
	}

  //calcolo colori
  printf("total codes %d\n",code);

  } ////////loop

  printf("max %d\n",maxv);
  if(0)
    for (int x=0;x<NB;x++)
      for (int y=0;y<NB;y++)
	for (int z=0;z<1 & z<NB;z++)
	  printf("%d %d %d: %d\n",x,y,z,B[x][y][z]);


  int la=1;
  CImg<double> stat(NB*la,NB*NB*la,1,3,255);
    for (int x=0;x<NB;x++)
      for (int y=0;y<NB;y++)
	for (int z=0;z<NB;z++){	  
	  float quanto=(float)B[x][y][z]/maxv;
	  //if (z==0&&quanto>0) printf("%f\n",quanto);
	  quanto*=10;
	  if (quanto>1) quanto=1;

	  for (int px=0;px<la;px++)
	    for (int py=0;py<la;py++){

	      int posx=x*la+px;
	      int posy=y*la+py+z*la*NB;
	      int flag=0;
	      if (z<NB/2)
		flag=1;
	      /*
	      int r=(int)(z*S+y*S);
	      int g=(int)(z*S-0.394*x*S-0.58*y*S);
	      int b=(int)(z*S+2.03*x*S);
	      */
	      int r=x*S;
	      int g=y*S;
	      int b=z*S;
	      stat(posx,posy,0,0)=(int)(r*quanto)+flag*(255*(1-quanto));
	      stat(posx,posy,0,1)=(int)(g*quanto)+flag*(255*(1-quanto));
	      stat(posx,posy,0,2)=(int)(b*quanto)+flag*(255*(1-quanto));

	      if (quanto==0){
	      stat(posx,posy,0,0)=128;
	      stat(posx,posy,0,1)=128;
	      stat(posx,posy,0,2)=128;
	      }
	    }
	}
  stat.save("stats.jpg");
  
  stats();
  exit(0);
  
  int bins=0;
  long thr=0;//tot/10000;
  vector<unsigned int> list;
  vector<long> weight;
  vector<float> px;
  vector<float> py;

  printf("thr %ld\n",thr);
  for (int x=0;x<NB;x++)
    for (int y=0;y<NB;y++)
      for (int z=0;z<NB;z++)
	if (A[x][y][z]>0){
	  // seleziona max locali
	  int max=1;
	  for (int dx=-1;dx<=1;dx++)
	    for (int dy=-1;dy<=1;dy++)
	      for (int dz=-1;dz<=1;dz++)
		if (dx!=0 || dy!=0 || dz!=0)
		  if (x+dx>=0 && y+dy>=0 && z+dz>=0 &&
		      x+dx<NB && y+dy<NB && z+dz<NB){
		    //printf(" %d %d %d d %d %d %d: %d %d\n",x,y,z,dx,dy,dz,A[x][y][z],A[x+dx][y+dy][z+dz]);
		    if (A[x][y][z]<A[x+dx][y+dy][z+dz] ||
			(A[x][y][z]==A[x+dx][y+dy][z+dz] &&
			 (dx<0 || dy<0 || dz<0)))
		      max=0;
		  }	  
	  if (max && A[x][y][z]>thr){
	    //	    printf("max %d %d %d: %ld\n",x,y,z,A[x][y][z]);	  
	    long v=0;
	    v+= ((long)x*S)<<16;
	    v+=((long)y*S)<<8;
	    v+=(long)z*S;
	    list.push_back(v);
	    weight.push_back(A[x][y][z]);	      
	    px.push_back(pos[x][y][z][0]);
	    py.push_back(pos[x][y][z][1]);
	    bins++;
	  }
	}


  int t;


  t=100;  // film 
  //sort: bubble
  int ll=list.size();
  int quanti=ll;
  do{
    printf("quanti %d, t %d\n",quanti,t);
  for (int i=0;i<ll;i++){

    // unisci colori vicini
    for (int j=0;j<ll-1;j++)
      if (i!=j){      
	int r=list[i]>>16;
	int g=(list[i]>>8) & 255;
	int b=(list[i]) & 255;
	int r1=list[j]>>16;
	int g1=(list[j]>>8) & 255;
	int b1=(list[j]) & 255;
	if ((r-r1)*(r-r1)+(g-g1)*(g-g1)+(b-b1)*(b-b1)<t && weight[i]+weight[j]>0){
	  long v=0;
	  v+= (long)((r*weight[i]+r1*weight[j])/(weight[i]+weight[j]))<<16;
	  v+=(long)((g*weight[i]+g1*weight[j])/(weight[i]+weight[j]))<<8;
	  v+=(long)(b*weight[i]+b1*weight[j])/(weight[i]+weight[j]);
	  list[i]=v;
	  weight[i]+=weight[j];
	  weight[j]=0;
	  px[i]+=px[j];
	  py[i]+=py[j];
	  px[j]=0;py[j]=0;
	}
      }
  }

  quanti=0;
  for (int i=0;i<ll;i++)
    if (weight[i]>0) quanti++;
  t=t*1.5;
  
  }while (quanti>20);

  weight.resize(quanti);
  list.resize(quanti);
  px.resize(quanti);
  py.resize(quanti);

  long tot=0;
  for (int i=0;i<ll;i++)
    tot+=weight[i];

  for (int i=0;i<ll;i++){
    for (int j=0;j<ll-1;j++){

      if (weight[j]<weight[j+1]){
	long t=list[j];
	list[j]=list[j+1];
	list[j+1]=t;
	t=weight[j];
	weight[j]=weight[j+1];
	weight[j+1]=t;
	float tt=px[i];
	px[i]=px[j];
	px[j]=tt;
	tt=py[i];
	py[i]=py[j];
	py[j]=tt;
      }
    }
  }

  for (int i=0;i<list.size();i++){
    int r=list[i]>>16;
    int g=(list[i]>>8) & 255;
    int b=(list[i]) & 255;
    printf("%d %d %d: %d %3.2f, %f %f\n",r,g,b,
	   weight[i],(float)weight[i]/tot*100,
	   px[i]/weight[i],py[i]/weight[i]);
  }

  int p=600;
  CImg<double> outi(2*p,p,1,3,255);

  for (int i1=0;i1<list.size();i1++){
    int r=list[i1]>>16;
    int g=(list[i1]>>8) & 255;
    int b=(list[i1]) & 255;
    float ref=pow((float)weight[0]/2,0.5f);
    int bound=(int)(pow((float)weight[i1],0.5f)/ref*p);
    for (int i=0;i<bound;i++)
      for (int j=0;j<bound;j++){
	int posx=2*p*(i1+0.5)/(list.size()-1)-bound/2+i;
	int posy=p/2-bound/2+j;
	if (posx>=0 && posy>=0 &&
	    posx<2*p && posy<p){
	  outi(posx,posy,0,0)=r;
	  outi(posx,posy,0,1)=g;
	  outi(posx,posy,0,2)=b;
	}
      }
  }

  outi.save("out.jpg");
    


  return 0;
}
