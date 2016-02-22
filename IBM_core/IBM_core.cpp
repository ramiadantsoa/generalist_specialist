#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>
#include <string>
#include <assert.h>

#include "dSFMT/dSFMT.h"

#define PI 3.1415926535897932385

#define chop_0 0.000000000000001

using namespace std;

typedef unsigned long int uInt;

const int rng_seed = (const int)time(NULL);
dsfmt_t rng;

void init_rng() {
    dsfmt_init_gen_rand(&rng, rng_seed);
}

#include "simple_functions.h"

/*CLASS DECLARATION AND DEFINITION*/
#include "class_parameter.h"
#include "class_species.h"
#include "class_grain.h"
#include "initialize_incidence.h"
#include "class_patch_dw.h"
#include "class_grid.h"
#include "class_patch_per_grid.h"
#include "class_colonization_per_grid.h"
#include "class_general_variable.h"

/*COLONIZATION AND TOTAL COLONIZATION FUNCTION AND RANDOM CHOOOSE PATCH*/

#include "colonization_functions.h"
#include "section_birth.h"
#include "section_death.h"
#include "section_colonization.h"

void simulation(Parameter *Initial_Parameter, double _w){
	string valiny; //valiny stores the lambda, tau , gammaH, final result as a fraction of occupancy and the duration

	ostringstream convert;

	//INITIALIZE SPECIES PARAMETER
	double col_a = Initial_Parameter->get_col_a();
	double diff_a = Initial_Parameter->get_diff_a();
	int M = Initial_Parameter->get_M();
	double l_s = Initial_Parameter->get_l_s();
	double muR = Initial_Parameter->get_muR();
	double T = Initial_Parameter->get_T();
	double w = _w;///EPS = 0 and is defined in class species 
	//double disp_d=Initial_Parameter->get_disp_d();
	double tau=Initial_Parameter->get_tau();
	double size =Initial_Parameter->get_size();
	
	//INITIALIZE SPECIES CHARACTERISTICS
	vector<Species*> SpChar(M);

	for(int k=0 ; k < M ; k++){
		Species *sp= new Species(col_a,diff_a,l_s,k,M,w);
		SpChar[k]=sp;
	}

	//double MAX_TAU = 10.0;
	//double INIT_LAMBDA= 0.0625;
	//double MAX_LAMBDA = 5;
	//int N_STEP_T =4;
	//int N_STEP_L = 4;
	int N_STEP_L = 15;
	//double lambda = _lambda; //radius of the grain or patch
	//double tau = _tau;

	// SIMLUATE ALL DIFFERENT LEVEL OF FRAGMENTATION AND DENSITY OF RESOURCE PRODUCTION
	double lambda, gammaH;
	for(int i = 0; i < N_STEP_L ; i++){
		
		lambda = (double) 1/(i*i + 1);//the largest radius of a patch is 3
		gammaH = 1/(lambda*lambda);
		//for(int j = 0; j < N_STEP_T ; j++){			
			//tau =(double) MAX_TAU*(j+1)/N_STEP_T;
			//lambda = (double) MAX_LAMBDA*(i+1)/N_STEP_L;
			convert<<tau<<"," <<lambda<<","<<gammaH; 

			cout<<"\n \n lambda is " <<lambda<< ", tau is "<< tau<< ", gammaH is "<< gammaH<<endl;
            cout.flush();

			//INITIALIZE HABITAT GRAIN: PATCHES LOCATION AND TYPE
			vector<Grain*> *habitat = initialize_habitat(Initial_Parameter,gammaH, lambda);

			//INITIALIZE SPECIES INCIDENCE
			double resource_density = PI*lambda*lambda*tau*gammaH/muR;

			vector<double> ocf = initialize_incidence(Initial_Parameter,resource_density, SpChar);
			
			//INITIALIZE GRID AND PATCH PER GRID 
			Grid *grid = initialize(Initial_Parameter, resource_density, ocf, habitat, SpChar[0]); 
			PatchPerGrid *ppg = initializePatchGrid(grid);

			//INITIALIZE COLONIZATION RATE AND COLONIZATION PER GRID
			double total_colonization = Total_colonization(grid,Initial_Parameter, SpChar);

			ColonizationPerGrid *cpg= initializeColonizationGrid(grid);

			//outputcpg(cpg);
			
			//INITIALIZE VARIABLE CONTAINING BIRTH DEATH AND COLONIZATION RATE
			double total_birth = lambda*lambda*PI*tau*gammaH*size*size;
			uInt npatches = grid->count_patches();			
			double total_death = npatches*muR;

			Variable *var = new Variable(total_birth, total_death, total_colonization, M,muR);
			
			initialize_temp_result(var,grid, M);//this contains the number of patches and species occupancy, it is a private of var.

			cout<<"number of resource units "<<npatches<<" initialization done\n";
			 
			double t = 0.0; 
            int half_T = (int) floor(Initial_Parameter->get_T()/2);
			double choose_event=0.0;
			double total_event=0.0;
			double dt=0;
            //double step = 0.0;
			vector<double> final_result(M+1); 
			for (int m = 0 ; m < M+1; m++) final_result[m]= 0.0;

			if(npatches == 0){
				for( int m = 0; m < M ; m++){
					final_result[m+1]= 0.0;
					cout<<final_result[m+1]<<"\t";
					convert<<","<<final_result[m+1];
				}
				convert<<"\n";	
				goto end;
			}

			while(t<T){
				total_event = var->get_birth()+ var->get_death() + var->get_col();//calculates the rate that an event happens
				dt = get_exp(total_event);
				t+=dt; 
				choose_event= total_event*get_random();

				if(choose_event <= var->get_birth()) {
					//cout<<"birth \n";
					Patch * new_p = new_patch(habitat,Initial_Parameter, lambda);
					update_birth(new_p, grid, ppg, cpg,Initial_Parameter, SpChar, var);
					
					//double totcol = Total_colonization(grid, Initial_Parameter, SpChar);
					//cout<<" total colonization var " << var->get_col() << " and from calculation " << totcol<<endl;
				}else{ 
					if(choose_event <= var->get_birth()+var->get_death()) {
					//cout<<"death \n ";
					vector<uInt> chosen_death = select_patch_death(grid,ppg, var);
					update_death(chosen_death, grid, ppg, cpg, Initial_Parameter, SpChar,var);
					npatches =grid->count_patches();
					if(npatches == 0){
						for( int m = 0; m < M ; m++){
							final_result[m+1]= 0.0;
							cout<<final_result[m+1]<<"\t";
							convert<<","<<final_result[m+1];
						}
						convert<<"\n";	
						goto end;
					}
					
					//double totcol = Total_colonization(grid, Initial_Parameter, SpChar);
					//cout<<" total colonization var " << var->get_col() << " and from calculation " << totcol<<endl;
					} 
					else {
					//cout<<"colonization \n " ;
					vector<uInt> chosen_col = select_patch_col(grid, cpg,var,M);
					update_colonization(chosen_col, grid, cpg, Initial_Parameter, var);
					update_colonization_rest(chosen_col, grid, cpg, Initial_Parameter, SpChar, var);

					//double totcol = Total_colonization(grid, Initial_Parameter, SpChar);
					//cout<<" total colonization var " << var->get_col() << " and from calculation " << totcol<<endl;
					}
				}

				if(t > half_T + final_result[0] ){
					//if (t < half_T+step+1){
						//step+= 1.0;
						for (int m = 0 ; m<M+1 ; m++){
							assert(var->fraction_occupancy(m)<=1);
							final_result[m]+= var->fraction_occupancy(m);//(double) var->get_temp_result(m+1)/var->get_temp_result(0);
						}
					//}
				}

				//for(int m = 0; m < M+1; m++) cout<< final_result[m]<<endl;
	
			}//while closing

			/*start testing if the colonization update works
			
			double tt_col, death_rate;

			cout<< " death rate "<< var->get_death()<< "\n" << "computed death rate "<< var->get_temp_result(0)<<"\n";  

			tt_col= total_colonization_check(grid,Initial_Parameter, SpChar);

			cout<<"total colonization "<< var->get_col()<<"\n"<<"long total colonization "<< tt_col<<"\n";
			

			system("pause");

			end testing*/

            cout << "Cleaning up\n";
            cout.flush();
			
			//cout<< final_result[0]<<"\n";
			output(grid);
			for( int m = 0; m < M ; m++){
				final_result[m+1]= final_result[m+1]/final_result[0];
				cout<<final_result[m+1]<<"\t";
				convert<<","<<final_result[m+1];
			}
			convert<<"\n";

			end:
			
			//grid->export_all_resource_units();

			delete grid;
            /* Delete everything in habitat */
            for(uInt i = 0; i < (*habitat).size(); i++) {
                delete (*habitat)[i];
            }
			delete habitat;
			delete cpg;
			delete ppg;
			delete var;
		//}
	}

    /* Free the species */

    for(int m=0 ; m < M ; m++){
        delete SpChar[m];
    }

	double time_passed;
	time_passed= clock();
	
	cout<<"\n sim time is : " << time_passed/CLOCKS_PER_SEC << "seconds "<<endl;

	convert<<"\n "<<time_passed/CLOCKS_PER_SEC<<endl;

	valiny.append(convert.str());

	ostringstream ltau;
	string addname;

	//ltau<<"_l_"<<lambda<<"_tau_"<<tau<<"_w_"<<w<<".csv";
	ltau<<"_w_"<<w<<".csv";
	addname = ltau.str();
	//string filename = Initial_Parameter->get_name().append(".csv");
	string filename = Initial_Parameter->get_name().append(addname);

	ofstream myfile;
	myfile.open (filename.c_str());
	myfile << valiny;
	myfile.close();
}

int main(){
	double compet_z =0.0, size = 6, time = 1000, muR = 0.1, col_a = 0.5, diff_a = 1, w, tau;
	short int est = 1, aggreg, M;
	double l_s = 0.5, lambda = 1;
	

	init_rng();

	//cout<<"Establishment 1 and Fecundity 0: ";
	//cin>>est;
	//cout<<"Compet z: ";
	//cin>>compet_z;
	cout<<"Aggregation: ";
	cin>>aggreg;
	diff_a = (aggreg == 1)? 1 : 0.8;
	//cout<<"lambda: ";
	//cin>>lambda;
	//cout<<"tau: ";
	//cin>>tau;
	//cout<<"l_s: ";
	//cin>>l_s;
	//cout<<"diff_a:";
	//cin>>diff_a;
	//cout<<"\n";
	//cout<<"Community size: ";
	//cin>> M;
	M = 5;
	w = 2.0;
	//double size0 = 4;
	//double time0 = 2000;
	//double muR0 = 0.1;
	//double gammaH0 = 1;
	
	cout<<"enter tau: ";
	cin>>tau;
	//cout<<"enter muR: ";
	//cin>>muR;
	cout<<"size: ";
	cin>>size;
	cout<<"time: ";
	cin>>time;
	//double time0, size0;
	//time0 = time;
	//size0 = size;

	//for(int ss = 2 ; ss < 3; ss++){
		//size = size0*pow(2.,ss);
		//for(int mm = 0; mm < 1 ; mm++){
		//muR = muR0*pow(2.,mm);
		
			//for(int gg = 0; gg < 1; gg++){
				//gammaH = gammaH0*pow(2.,gg);
				//time=500;
				//for(int tt = 1; tt < 3; tt++){
					//time = time0*pow(2.,tt);
					for(int rep = 0 ; rep < 4 ; rep++){
						Parameter *param = new Parameter(size, M ,l_s, compet_z, time, est, aggreg, muR, rep,col_a,diff_a,tau);
						simulation(param,w);
						delete param;
					}
				//}
			//}
		//}
	//}

	return 0;
}
