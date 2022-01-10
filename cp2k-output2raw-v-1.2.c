/*BY FILIPE MATUSALEM, SEPT 2021     filipematus@gmail.com */
/*Program to extract raw data from CP2K output*/
/*Compilation: g++ -o program.x program.c*/

//units to use in deepMD-kit

//Property	Unit
//Time	ps
//Length	Å
//Energy	eV
//Force	eV/Å
//Virial	eV
//Pressure	Bar





#include <stdio.h>
#include <getopt.h> // *GNU* Para o getopt_long()
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main(int argc, char *argv[])
{
FILE *cp2kout,*stress_tensor,*forces,*energies;
int i,j,k,l,nsteps,fcr,stress,istep;
char str1[150],ch,lixo[150];
double version,gpa2bar,energy,a,b,c,au2ev,au2jperm,j2ev,x[3],y[3],z[3],volume;

fcr=0;
stress=0;
gpa2bar=10000;
au2ev=2.72113838565563E+01;
au2jperm=8.23872205491840E-08;    //1 atomic unit of force = 8.238722514141 x 10-8 joule/meter
j2ev=6.242E+18;

if( argc < 2 ){
printf("\n\n");
printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR\n");
printf("EE        RRRR       RRRR       OO   OO   RRRR\n");
printf("EE        RR  RR     RR  RR     OO   OO   RR  RR\n");
printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR\n\n");

printf("Enter the name of cp2k output file. \n\n");

 exit(0);}

strcpy(str1,argv[1]);

cp2kout = fopen(str1,"r"); /* Arquivo ASCII, para leitura */
if(!cp2kout)
{
printf( "Error opening argument 1 file\n");

exit(0);
}

energies = fopen("energy.raw","w"); /* Arquivo ASCII, para escrita */
if(!energies)
{
printf( "Error creating energy.raw file\n");
exit(0);
}


//find cp2k version
do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"version")!=0);
do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"version")!=0);
fscanf(cp2kout,"%lf",&version);

printf("\nCP2K version: %.2f\n\n",version);
//********************************


//test if forces are written
do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"ENERGY|")!=0);
do
ch = getc(cp2kout);              /*chega ao fim da linha*/
while(ch!='\n');
fscanf(cp2kout,"%s",lixo);
if(strcmp(lixo,"ATOMIC")==0)fcr=1;
//printf("%d\n",fcr);


if(fcr==1){

printf("\nForces found!! force.raw file will be write!!\n");


forces = fopen("force.raw","w"); /* Arquivo ASCII, para escrita */
if(!forces)
{
printf( "Error creating forces.raw file\n");
exit(0);}}
//*********************************** 

//test if stress are written
while (fscanf(cp2kout,"%s",str1) != EOF){            /*conta steps*/
if(strcmp(str1,"STRESS|")==0){stress=1; break;}}
//printf("%d\n",stress);
rewind(cp2kout);

if(stress==1){

printf("\nStress found!! virial.raw file will be write!! \n");


stress_tensor = fopen("virial.raw","w"); /* Arquivo ASCII, para escrita */
if(!stress_tensor)
{
printf( "Error creating virial.raw file\n");
exit(0);
}

do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"volume")!=0);

do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"[ang^3]")!=0);

fscanf(cp2kout,"%lf",&volume);
//printf("%lf\n",volume);

rewind(cp2kout);
}

//test if initial step is set to zero
do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"MD|")!=0);
do
ch = getc(cp2kout);              /*chega ao fim da linha*/
while(ch!='\n');
fscanf(cp2kout,"%s",lixo);
fscanf(cp2kout,"%s",lixo);
fscanf(cp2kout,"%s",lixo);
fscanf(cp2kout,"%d",&istep);

//printf("istep %d\n",istep);
if(istep!=1){printf("\n\nWARNING!! Initial step is not zero!! CP2K do not dump to trajectory the initial step calculated => one step is discarted!!\n\n"); istep=-2;}
else {istep=-1;}
//printf("istep %d\n",istep);
rewind(cp2kout);

nsteps=istep;
while (fscanf(cp2kout,"%s",str1) != EOF){            /*conta steps*/
if(strcmp(str1,"ENERGY|")==0)nsteps++;                      
}

printf("\nTotal number of MD steps = %d\n",nsteps);
rewind(cp2kout);

double tensor[3][3];

for(i=0;i<-1*istep;i++){
do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"ENERGY|")!=0);
}







k=0;
for(j=0;j<nsteps;j++){

do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"ENERGY|")!=0);

fscanf(cp2kout,"%s",lixo);fscanf(cp2kout,"%s",lixo);fscanf(cp2kout,"%s",lixo);fscanf(cp2kout,"%s",lixo);
fscanf(cp2kout,"%s",lixo);fscanf(cp2kout,"%s",lixo);fscanf(cp2kout,"%s",lixo);

//printf("%s\n",lixo);
fscanf(cp2kout,"%lf",&energy);

fprintf(energies,"%24.18lE\n",energy*au2ev); 

//----if forces are present-----
if(fcr==1){


do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"Atom")!=0);

do
ch = getc(cp2kout);              /*chega ao fim da linha*/
while(ch!='\n');


fscanf(cp2kout,"%s",lixo);

do
{fscanf(cp2kout,"%s",lixo);fscanf(cp2kout,"%s",lixo);
fscanf(cp2kout,"%lf",&a); 
fscanf(cp2kout,"%lf",&b); 
fscanf(cp2kout,"%lf",&c);
fprintf(forces,"%24.18lE ",a*au2jperm*j2ev/1E10); 
fprintf(forces,"%24.18lE ",b*au2jperm*j2ev/1E10); 
fprintf(forces,"%24.18lE ",c*au2jperm*j2ev/1E10);
fscanf(cp2kout,"%s",lixo);}
while(strcmp(lixo,"SUM")!=0);



fprintf(forces,"\n");


}
//*************************************

if(stress==1){



//cp2k version >=8.1
if(version>8){
do
fscanf(cp2kout,"%s",str1);                                      /*posiciona o  após a palavra xxx*/
while(strcmp(str1,"STRESS|")!=0);

for(i=0;i<2;i++){
do
ch = getc(cp2kout);              /*chega ao fim da linha*/
while(ch!='\n');
}

for(i=0;i<3;i++){
fscanf(cp2kout,"%s",lixo);fscanf(cp2kout,"%s",lixo);
fscanf(cp2kout,"%lf",&tensor[i][0]);
fscanf(cp2kout,"%lf",&tensor[i][1]);
fscanf(cp2kout,"%lf",&tensor[i][2]);
}
}
else {
for(i=0;i<6;i++){
do
ch = getc(cp2kout);              /*chega ao fim da linha*/
while(ch!='\n');
}

for(i=0;i<3;i++){
fscanf(cp2kout,"%s",lixo);
fscanf(cp2kout,"%lf",&tensor[i][0]);
fscanf(cp2kout,"%lf",&tensor[i][1]);
fscanf(cp2kout,"%lf",&tensor[i][2]);}
}

// stress to virial conversion, default unit for stress in cp2k is GPa
// note the stress is virial (J) = stress (Pa)* volume (m^3) ==> virial (J) = stress (GPa) * 1.0E+9 * volume (Ang^3) * 1.0E-30 ==>  virial (J) = stress (GPa) * volume (Ang^3) * 1.0E-21 ==> virial (eV) = stress (GPa) * volume (Ang^3) * 6.242E-3 

fprintf(stress_tensor,"%24.18lE %24.18lE %24.18lE %24.18lE %24.18lE %24.18lE %24.18lE %24.18lE %24.18lE ",tensor[0][0]*volume*6.242E-3,tensor[0][1]*volume*6.242E-3,tensor[0][2]*volume*6.242E-3,tensor[1][0]*volume*6.242E-3,tensor[1][1]*volume*6.242E-3,tensor[1][2]*volume*6.242E-3,tensor[2][0]*volume*6.242E-3,tensor[2][1]*volume*6.242E-3,tensor[2][2]*volume*6.242E-3);

fprintf(stress_tensor,"\n");

}
}



printf("\n\nEnergy in eV written to energy.raw!\n");
if(stress==1)printf("Virial tensor components in eV (xx xy xz yx yy yz zx zy zz) written to virial.raw!\n\n");
if(fcr==1)printf("Forces in ev/Ang written to force.raw!\n\n");

printf("-----------------------------------------------------------------------------------------------\n\n");


if(fcr==1)fclose(forces);
if(stress==1)fclose(stress_tensor);
fclose(cp2kout);
fclose(energies);

}
