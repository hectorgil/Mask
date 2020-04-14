#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define Pi (4.*atan(1.))
typedef struct {

        double OMEGA_M;
        
}f_params;

double Leg2(double x);

double Leg4(double x);

double Leg6(double x);

double Leg8(double x);

void z_to_r(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

long int countlines(char *filename);

void freeTokens(double** tokens, int N);

int main(int argc, char *argv[])
{
long int i,j,k,l;
double x1,y1,z1;
FILE *f,*g;	
double Omega_matter;	
char path[200],basename1[200],id[200],path_out[200],W_out[200],Norm_out[200],input_fname_ran[200];
char name_file_ini[2000];
char yamamoto[2000];
double* s_x_ran;
double* s_y_ran;
double* s_z_ran;
double* weight_ran;
double* density;
double trash,density1,density2;
double AR,dec,redshift,radial;
double precision;
double weight_fkp,weight_cp, weight_noz, weight_sys,weight_col,weight_tot;

double theta;
double MIN1[1];
double MAX1[1];

long int npar_ran;
long int npar_ran_used;
double npar_ran_eff;
double z_min,z_max;
double random,random_cut;
int veto;
double distance_max;

double s;

double delta_s;
int N_bin;
int index_s;
double mu;

int n_lines_parallel,tid;

double **W0;
double **W2;
double **W4;
double **W6;
double **W8;
double **s_eff;
double **mu_eff;
double **num_eff;
double n_z;
double weight_ij;
double xlos,ylos,zlos;
double norm;
int bin_at_normalization;
char wfkp_cmass_c[200],wfkp_eboss_c[200];
double wfkp_cmass,wfkp_eboss;
char in_eboss_foot[20],iscmass[200];

//inizialization

sprintf(name_file_ini,argv[1]);
f=fopen(name_file_ini,"r");
if(f==NULL){printf("File %s not found...\t exiting now\n",name_file_ini);return 0;}
else{printf("Reading Inizialization file: %s\n\n",name_file_ini);}
fscanf(f,"%*s\n");
fscanf(f,"%*s %*s %*s %*s %s\n",path);
fscanf(f,"%*s %*s %s\n",basename1);
fscanf(f,"%*s %*s %*s %*s %s\n",path_out);
fscanf(f,"%*s %*s %s\n",id);
fscanf(f,"%*s %*s %lf %lf\n",&z_min,&z_max);
if(z_min<0 || z_max<0 || z_min>=z_max){printf("Warning! Unusual values of z_min=%lf, z_max=%lf. Exiting now...\n",z_min,z_max);exit(0);}
fscanf(f,"%*s %*s %*s %lf\n",&Omega_matter);
if(Omega_matter<=0 || Omega_matter>1){printf("Warning! Unusual value of Om=%lf. Exiting now...\n",Omega_matter);exit(0);}
fscanf(f,"%*s %*s %*s %d\n",&bin_at_normalization);
if(bin_at_normalization<1 || bin_at_normalization>100){printf("Warning! Unusual value of bin_norm=%d. Exiting now...\n",bin_at_normalization);exit(0);}
fscanf(f,"%*s %*s %*s %*s %lf\n",&delta_s);
if(delta_s<=0){printf("Warning! Unusual valoe of Delta_S=%lf. Exiting now...\n",delta_s);exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %lf\n",&distance_max);
if(distance_max<0 || distance_max<delta_s*2){printf("Warning! Unusual value of maximum distance %lf. Exiting now.../n",distance_max);exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %lf\n",&random_cut);
if(random_cut<=0 || random_cut>100){printf("Warning. Percentage of randoms (%lf) cannot be negative or larger than 100. Exiting now.../n",random_cut);exit(0);}
fscanf(f,"%*s %*s %*s %s\n",yamamoto);
if(strcmp(yamamoto, "yes") != 0 && strcmp(yamamoto, "no") != 0  ){printf("Yamamoto aproximation option either yes or no: %s. Exiting now...\n",yamamoto);exit(0);}
fclose(f);


N_bin=(int)(distance_max/delta_s);
printf("Nbins=%d\n",N_bin);
printf("Omega_m=%lf\n",Omega_matter);

random_cut=random_cut/100.;
f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_matter;


        sprintf(W_out,"%s/W_%s.txt",path_out,id);

        sprintf(input_fname_ran, "%s/%s", path, basename1);
        printf("Reading %s...\n",input_fname_ran);
        npar_ran=countlines(input_fname_ran);          

	s_x_ran = (double*) calloc(npar_ran, sizeof(double));
	s_y_ran = (double*) calloc(npar_ran, sizeof(double));
	s_z_ran = (double*) calloc(npar_ran, sizeof(double));
	weight_ran = (double*) calloc(npar_ran, sizeof(double));
        density =  (double*) calloc(npar_ran, sizeof(double));
	
	
	f=fopen(input_fname_ran,"r");

//	srand48(time(NULL));
        srand48(12345678910111213);
	npar_ran_used=0;
        npar_ran_eff=0;
	for(i=0;i<npar_ran;i++)
	{

		//You may need to modify this part acccording to your catalogue
		veto=1;weight_col=1;weight_cp=1;weight_noz=1;weight_sys=1;weight_fkp=1;n_z=1;
        
fscanf(f,"%lf %lf %lf %s %lf %lf %lf %lf %*s %s %s %s %*s %lf\n",&AR,&dec,&redshift,wfkp_eboss_c,&weight_sys,&weight_cp,&weight_noz,&n_z,in_eboss_foot,iscmass,wfkp_cmass_c,&weight_fkp);
if(strcmp(iscmass,"T") == 0){weight_col=weight_cp+weight_noz-1;}
else{weight_col=weight_cp*weight_noz;}

         weight_tot=weight_fkp*weight_sys*weight_col;
               
		if(i==0){printf("\n");}
		if(i<5 || i>npar_ran-5){printf("%lf %lf %lf %lf\n",AR,dec,redshift,weight_tot);}
		if(i==5){printf("   ...\n");}
		theta=90.-dec;
                random=drand48();
		if(redshift>z_min && redshift<z_max && veto==1 && weight_tot>0 && random<random_cut)
		{
			
			//passem a radiants
			AR=AR*Pi/180.;
			dec=dec*Pi/180.;
			theta=theta*Pi/180.;	
						
			MAX1[0]=redshift;
			MIN1[0]=0;
			
			precision=1e-6;
			
			adapt_integrate(1, z_to_r , function_parameters, 1, MIN1, MAX1 ,100000, precision, precision, &radial, &trash);
										
			x1=radial*sin(theta)*cos(AR);
			y1=radial*sin(theta)*sin(AR);
			z1=radial*cos(theta);
						
			//Sin Girar
			s_x_ran[npar_ran_used]=x1;
			s_y_ran[npar_ran_used]=y1;
			s_z_ran[npar_ran_used]=z1;
			weight_ran[npar_ran_used]=weight_tot;			
			density[npar_ran_used]=n_z;
			npar_ran_used++;
                        npar_ran_eff=npar_ran_eff+weight_tot;
			
		}
		
	}
	fclose(f);

        
        #pragma omp parallel for private(i,tid) shared(n_lines_parallel,npar_ran_used)
        for(i=0;i<npar_ran_used;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){n_lines_parallel=omp_get_num_threads();}
        }
printf("%d processors used\n",n_lines_parallel);
 


W0=  (double**) calloc(N_bin, sizeof(double*));
W2= (double**) calloc(N_bin, sizeof(double*));
W4= (double**) calloc(N_bin, sizeof(double*));
W6= (double**) calloc(N_bin, sizeof(double*));
W8= (double**) calloc(N_bin, sizeof(double*));
s_eff= (double**) calloc(N_bin, sizeof(double*));
mu_eff= (double**) calloc(N_bin, sizeof(double*));
num_eff= (double**) calloc(N_bin, sizeof(double*));


for(i=0;i<N_bin;i++)
{
W0[i] = (double*) calloc(n_lines_parallel, sizeof(double));
W2[i] = (double*) calloc(n_lines_parallel, sizeof(double));
W4[i] = (double*) calloc(n_lines_parallel, sizeof(double));
W6[i] = (double*) calloc(n_lines_parallel, sizeof(double));
W8[i] = (double*) calloc(n_lines_parallel, sizeof(double));
s_eff[i] = (double*) calloc(n_lines_parallel, sizeof(double));
mu_eff[i] = (double*) calloc(n_lines_parallel, sizeof(double));
num_eff[i] = (double*) calloc(n_lines_parallel, sizeof(double));
}

printf("Starting the parallel loop...\n");
#pragma omp parallel for private(density1,density2,i,j,tid,s,xlos,ylos,zlos,mu,index_s,weight_ij) shared(s_x_ran,s_y_ran,s_z_ran,npar_ran_used,W0,W2,W4,W6,W8,delta_s,N_bin,s_eff,mu_eff,weight_ran,num_eff,yamamoto,density)
for(i=0;i<npar_ran_used;i++)
{
tid=omp_get_thread_num();//thread number
density1=density[i];
for(j=i;j<npar_ran_used;j++)
{
density2=density[j];
weight_ij=weight_ran[i]*weight_ran[j];
s=sqrt( pow(s_x_ran[i]-s_x_ran[j],2)+pow(s_y_ran[i]-s_y_ran[j],2)+pow(s_z_ran[i]-s_z_ran[j],2) );

//observer at 0
if(strcmp(yamamoto, "no") == 0){
xlos=s_x_ran[j]+(s_x_ran[i]-s_x_ran[j])/2.;
ylos=s_y_ran[j]+(s_y_ran[i]-s_y_ran[j])/2.;
zlos=s_z_ran[j]+(s_z_ran[i]-s_z_ran[j])/2.;
}
else{
xlos=s_x_ran[j];//+(s_x_ran[i]-s_x_ran[j])/2.;
ylos=s_y_ran[j];//+(s_y_ran[i]-s_y_ran[j])/2.;
zlos=s_z_ran[j];//+(s_z_ran[i]-s_z_ran[j])/2.;
}

mu=((s_x_ran[i]-s_x_ran[j])*xlos+(s_y_ran[i]-s_y_ran[j])*ylos+(s_z_ran[i]-s_z_ran[j])*zlos)/(s*sqrt(xlos*xlos+ylos*ylos+zlos*zlos));
if(i==j){mu=0;}
index_s=(int)(s/delta_s);
if(index_s<N_bin)
{
num_eff[index_s][tid]=num_eff[index_s][tid]+1.;
W0[index_s][tid]=W0[index_s][tid]+weight_ij/(s*s*delta_s);
W2[index_s][tid]=W2[index_s][tid]+weight_ij*Leg2(mu)*5./(s*s*delta_s);
W4[index_s][tid]=W4[index_s][tid]+weight_ij*Leg4(mu)*9./(s*s*delta_s);
W6[index_s][tid]=W6[index_s][tid]+weight_ij*Leg6(mu)*13./(s*s*delta_s);
W8[index_s][tid]=W8[index_s][tid]+weight_ij*Leg8(mu)*17./(s*s*delta_s);
s_eff[index_s][tid]=s_eff[index_s][tid]+s;
mu_eff[index_s][tid]=mu_eff[index_s][tid]+mu;
}

}

}


free(weight_ran);
free(s_x_ran);
free(s_y_ran);
free(s_z_ran);
printf("End of parallel loop.... Reassigning sectors\n");
        
        for(i=0;i<N_bin;i++)
        {
                for(j=1;j<n_lines_parallel;j++)
                {
                        W0[i][0]=W0[i][0]+W0[i][j];
                        W2[i][0]=W2[i][0]+W2[i][j];
                        W4[i][0]=W4[i][0]+W4[i][j];
                        W6[i][0]=W6[i][0]+W6[i][j];
                        W8[i][0]=W8[i][0]+W8[i][j];
                        s_eff[i][0]=s_eff[i][0]+s_eff[i][j];
                        mu_eff[i][0]=mu_eff[i][0]+mu_eff[i][j];
                        num_eff[i][0]=num_eff[i][0]+num_eff[i][j];

                }
       }

printf("Complete!\n");

norm=W0[bin_at_normalization-1][0];

f=fopen(W_out,"w");
fprintf(f,"#N= %e Sumwtot= %e\n",norm,npar_ran_eff);
for(i=1;i<N_bin;i++)
{
fprintf(f,"%e %e %e %e %e %e %e %e\n",(i+0.5)*delta_s,(s_eff[i][0]/num_eff[i][0]),W0[i][0]/norm,W2[i][0]/norm,W4[i][0]/norm,W6[i][0]/norm,W8[i][0]/norm,(mu_eff[i][0]/num_eff[i][0]));

}
fclose(f);

printf("Exiting...\n");

return 0;
	

}


void z_to_r(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;//cast from void to f_params
        double omega= params_function.OMEGA_M;
        double z=x[0];
        double c=299792.458;

        fval[0]=c/(100.*sqrt(omega*pow(1+z,3)+1.-omega));


}


double Leg2(double x)
{
double f;
f=0.5*(3.*x*x-1.);
return f;
}

double Leg4(double x)
{
double f;
f=1./8.*(35.*pow(x,4)-30.*pow(x,2)+3.);
return f;
}

double Leg6(double x)
{
double f;
f=1./16.*(231.*pow(x,6)-315.*pow(x,4)+105*pow(x,2)-5.);
return f;
}


double Leg8(double x)
{
double f;
f=1./128.*(6435.*pow(x,8)-12012.*pow(x,6)+6930.*pow(x,4)-1260.*pow(x,2)+35.);
return f;
}


long int countlines(char *filename)
{
FILE* myfile = fopen(filename, "r");
if(myfile==NULL){printf("Error reading %s. Exiting now...\n",filename);exit(0);}
long int ch, number_of_lines = -1;

do
{
            ch = fgetc(myfile);
                    if(ch == '\n')
                                        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0)
    number_of_lines++;

        fclose(myfile);
        return number_of_lines;
}

void freeTokens(double **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);

}

