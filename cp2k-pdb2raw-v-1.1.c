//  BY FILIPE MATUSALEM, AUG 2020     filipematus@gmail.com 
//  Program to convert CP2K PDB trajectory file to VASP XDATCAR format
//  Compilation: g++ -o cp2k-pdb2xdatcar cp2k-pdb2xdatcar.c   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846


int main(int argc, char *argv[])
{
FILE *pdb,*box,*coord, *type, *type_map, *energy;
double vecx[3],vecy[3],vecz[3],alpha,beta,gamma,x,y,z,a,b,c;
int i,j,k,natoms,nspecies,nsteps,ntype[10],pula;
char str[150],str1[150],ch,species[10][10];


 if( argc < 2 ){
printf("\n\n");
printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR\n");
printf("EE        RRRR       RRRR       OO   OO   RRRR\n");
printf("EE        RR  RR     RR  RR     OO   OO   RR  RR\n");
printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR\n\n");

printf("Enter the name of the cp2k trajectory pdb file \n\n");

 exit(0);}


printf("-----------------------------------------------------------------------------------------------\n\n");
printf("A second argument N (integer) defines if every frame is include or only between each N ones (Default N=1).\n\n");
if( argc > 2 ){pula = atoi(argv[2]);}
else {pula=1;}


strcpy(str1,argv[1]);

pdb = fopen(str1,"r"); /* Arquivo ASCII, para leitura */
if(!pdb)
{
printf( "Error opening argument 1 file\n");
exit(0);
}

coord = fopen("coord.raw","w"); /* Arquivo ASCII, para escrita */
if(!coord)
{
printf( "Error creating coord.raw file\n");
exit(0);
}

box = fopen("box.raw","w"); /* Arquivo ASCII, para escrita */
if(!box)
{
printf( "Error creating box.raw file\n");
exit(0);
}

type = fopen("type.raw","w"); /* Arquivo ASCII, para escrita */
if(!type)
{
printf( "Error creating type.raw file\n");
exit(0);
}

type_map = fopen("type_map.raw","w"); /* Arquivo ASCII, para escrita */
if(!type_map)
{
printf( "Error creating type.raw file\n");
exit(0);
}

energy = fopen("energy_from_pos.raw","w"); /* Arquivo ASCII, para escrita */
if(!energy)
{
printf( "Error creating energy_from_pos.raw file\n");
exit(0);
}

do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

natoms=0;
do
{fscanf(pdb,"%s",str1);  if(strcmp(str1,"ATOM")==0)natoms++;  }                                   
while(strcmp(str1,"END")!=0);

printf("No atoms %d\n",natoms);
rewind(pdb);


do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",species[0]);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

j=k=1;
for(i=0;i<natoms-1;i++){
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",species[j]);
k++;
if(strcmp(species[j-1],species[j])!=0){ntype[j-1]=k-1;k=1;j++;}

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');
}
ntype[j-1]=k;
nspecies=j;

printf("No. Species %d \n\n",nspecies);
printf("Specie  Number\n");
for(i=0;i<nspecies;i++){
printf("   %s      %d \n",species[i],ntype[i]);

fprintf(type_map,"%s\n",species[i]);

for(j=0;j<ntype[i];j++)fprintf(type,"%d ",i);


}





rewind(pdb);

nsteps=0;
while (fscanf(pdb,"%s",str1) != EOF){            /*conta steps*/
if(strcmp(str1,"CRYST1")==0)nsteps++;                      
}

printf("No. steps = %d\n\n",nsteps-1);
printf("Print every %d steps\n",pula);
rewind(pdb);

float M[3][3],V;

do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);




for(i=0;i<nsteps-1;i++){
do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"E")!=0);
fscanf(pdb,"%s",str1);
fscanf(pdb,"%lf",&a);fprintf(energy,"%24.18lE\n",a*2.72113838565563E+01);

do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);

fscanf(pdb,"%lf",&a);
fscanf(pdb,"%lf",&b);
fscanf(pdb,"%lf",&c);
fscanf(pdb,"%lf",&alpha);
fscanf(pdb,"%lf",&beta);
fscanf(pdb,"%lf",&gamma);


alpha=PI*alpha/180;
beta=PI*beta/180;
gamma=PI*gamma/180;

//conversion from fractional coordinates system to cartesian
//http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm

V=a*b*c*sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)+2*cos(alpha)*cos(beta)*cos(gamma));

//transformation matrix M

M[0][0]=a;
M[0][1]=b*cos(gamma);
M[0][2]=c*cos(beta);

M[1][0]=0;
M[1][1]=b*sin(gamma);
M[1][2]=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);

M[2][0]=0;
M[2][1]=0;
M[2][2]=V/(a*b*sin(gamma));

//Vector_x = M * (1 0 0)^t
//Vector_y = M * (0 1 0)^t
//Vector_z = M * (0 0 1)^t

for(j=0;j<3;j++){
vecx[j]=M[j][0]*1+M[j][1]*0+M[j][2]*0;
vecy[j]=M[j][0]*0+M[j][1]*1+M[j][2]*0;
vecz[j]=M[j][0]*0+M[j][1]*0+M[j][2]*1;}

//include or not every step
if(i%pula == 0){



fprintf(box,"%24.18lE %24.18lE %24.18lE ",vecx[0],vecx[1],vecx[2]);
fprintf(box,"%24.18lE %24.18lE %24.18lE ",vecy[0],vecy[1],vecy[2]);
fprintf(box,"%24.18lE %24.18lE %24.18lE\n",vecz[0],vecz[1],vecz[2]);

for(j=0;j<natoms;j++){
do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

fscanf(pdb,"%s",str1);fscanf(pdb,"%s",str1);fscanf(pdb,"%s",str1);
fscanf(pdb,"%lf",&a);fprintf(coord,"%24.18lE ",a);
fscanf(pdb,"%lf",&a);fprintf(coord,"%24.18lE ",a);
fscanf(pdb,"%lf",&a);fprintf(coord,"%24.18lE ",a);
}
}

fprintf(coord,"\n");



}


printf("\n");
printf("coord written!! \n");
printf("-----------------------------------------------------------------------------------------------\n\n");


fclose(type);
fclose(type_map);
fclose(pdb);
fclose(coord);
fclose(box);
}
