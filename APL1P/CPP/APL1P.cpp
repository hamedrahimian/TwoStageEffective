/*
******************************************
* Author:   Hamed Rahimian
*          The Ohio State University
*          March 01, 2016
******************************************
*/
// APL1P


#include <ilcplex/ilocplex.h>
#include <vector>
#include <array>
#include <algorithm> 


ILOSTLBEGIN

typedef IloArray<IloCplex> CPXarray;
typedef IloArray<IloModel> MODELarray;

struct Scenario{
	IloInt No;
	int * Realization;
	IloNum new_ObjCost; //f_new(x^*)
	IloNum Cost;
	IloNum nom_Prob;
	IloNum worst_Prob;
	char * Status;
	IloNum partial_sum;
};


struct SubScenario{
	IloNum WorstCost;
	IloNum nom_Prob;
	IloNum partial_sum;
};

enum availability_demand_matrix{
	low_low=0,
	low_med,
	low_high,
	med_low,
	med_med,
	med_high,
	high_low,
	high_med,
	high_high };

struct Empirical_Analysis{
	double * emp_prob_effective;
	double * emp_prob_ineffective;
};

int main(int argc, char **argv) {

 IloEnv envData;

 try {

	 //read differenet values of rho
	string line;
	ifstream riskpar ("config.txt");
	vector<double> vect;
	int kk=0;
	if (riskpar.is_open()){
		while ( getline (riskpar,line) ){
			kk++;
			stringstream ss(line);

			double i;

			while (ss >> i){
				vect.push_back(i);
				if (ss.peek() == ' ')
					ss.ignore();
			}
			
		}
	}

	

	riskpar.close();
	const int No_instance=kk;

	double config[21];
	
	for(int i=0; i<vect.size(); i++){
		config[i]=vect[i];
	}
	

	 IloNum toler = 0.000001;
	 IloInt i, j; 
	 IloInt w, w1, w2, w3, w4, w5, ww; 

		 
	 IloInt noGen, noDem; 
	 
	 IloNumArray cost(envData), us(envData), ccmin(envData), ccmax(envData);
	 IloNumArray2 nom_prob(envData), demand(envData), f(envData), alpha(envData);

///////////////// DATA FILE READING ////////////////////////////////

	  const char* filename  = "inputs.dat";	      

	  ifstream file(filename);

	  if (!file) {
			cerr << "ERROR: can not open file '" << filename << "' for reading" << endl;
			cerr << "usage:   " << argv[0] << " <file>" << endl;
			throw(-1);
		}

      file >> ccmin >> ccmax >> cost >> f >> us>> nom_prob>> alpha >> demand;

	  noGen = cost.getSize();
	  noDem = us.getSize();
	  
	  IloIntArray noScen(envData, noGen+noDem);
	  
	  for(i=0; i<noGen+noDem; i++) {
		  noScen[i] = nom_prob[i].getSize(); 
	  }

	  IloBool consistentData = (ccmin.getSize() == noGen &&  alpha[0].getSize() == noScen[0] && alpha[1].getSize() == noScen[1] && demand[0].getSize() == noScen[2] && demand[1].getSize() == noScen[3] && demand[2].getSize() == noScen[4]);

      if (!consistentData) { 
 	      cerr << "ERROR: Inconsistent data!" << endl;
          throw(-1);
      }
	  
	  file.close();
	  
	  const IloInt numScen=noScen[0]*noScen[1]*noScen[2]*noScen[3]*noScen[4];
	  IloNumArray prob(envData, numScen);
	  

	  vector<array<int, 5>> realization;
	  w=0;
	  for(w1=0; w1<noScen[0]; w1++) {
		  for(w2=0; w2<noScen[1]; w2++){
			  for(w3=0; w3<noScen[2]; w3++) {	
				  for(w4=0; w4<noScen[3]; w4++) {	
					  for(w5=0; w5<noScen[4]; w5++) {	
						//prob[w]=nom_prob[0][w1]*nom_prob[1][w2]*nom_prob[2][w3]*nom_prob[3][w4]*nom_prob[4][w5];
						array<int, 5> temp={w1, w2, w3, w4, w5};
						realization.push_back(temp);
						
						w++;
					  }
				  }
			  }
		  }
	  }

	  
	  const char* filename2 = "ScenFac.txt";
	  ifstream file2(filename2);
	  file2 >> prob;
	  file2.close();
	  


	  //These are used in checking nestedness
	  

	  //vector to stores each gamma ineffective scenarios
	  vector<int> * Ineffscen= new vector<int> [21];

	  //vector to stores each gamma effective scenarios
	  vector<int> * Effscen= new vector<int> [21];

	  //vector to stores each gamma undetermined scenarios
	  vector<int> * Unknownscen= new vector<int> [21];

	  //vector to stores each gamma undetermined scenarios
	  Empirical_Analysis * Interpretation= new Empirical_Analysis [21];

	  //vector to stores each gamma new cost after removal of an effective scenarios
	  vector<double> * NewCost_Effscen= new vector<double> [21];


	  IloNumArray Count_Effective(envData, numScen);
	  IloNumArray Count_Known(envData, numScen);

	  for (w=0; w<numScen; w++){
		  Count_Effective[w]=0;
		  Count_Known[w]=0;
	  }


	 // const char* Statusfilename  = "Status.txt";
	 // ofstream Summary(Statusfilename);
	 // Summary<<"Status of scenarios.."<<endl;
	 // Summary<< '\t';
	 // 
	 // for (w=0; w<numScen-1; w++){
		//Summary << w<< " ";
	 // }
	 // Summary<<numScen-1<<endl;

	  const char* Detfilename = "EX_Details.txt";
	  ofstream Details(Detfilename);
	  Details << "n= " << numScen << endl;

	  Details << '\t' << "X*" << '\t' << "VaR" << '\t' << "Sup"<< endl;


	  const char* Numberfilename  = "EX_Numbers.txt";
	  ofstream Numbers(Numberfilename);
	  Numbers<<"Number of scenarios.."<<endl;
	  Numbers<< '\t';
	  Numbers<< "I"<< '\t'<<"E"<<'\t'<<"U"<<endl;

	  const char* Number_post_filename  = "EX_Numbers_post.txt";
	  ofstream Numbers_post(Number_post_filename);
	  Numbers_post<<"Number of scenarios.."<<endl;
	  Numbers_post<< '\t';
	  Numbers_post<< "I"<< '\t'<<"E"<<'\t'<<"U"<<endl;
	 
	  const char* Scenariofilename  = "EX_Average.txt";
	  ofstream ScenarioRes(Scenariofilename);
	  
	  const char* Eff_file_name  = "EX_Eff_Change.txt";
	  ofstream EffScenarioRes(Eff_file_name);
	  for (i=0; i<No_instance; i++){
		  EffScenarioRes<< config[i]<< '\t';
	  }
	  EffScenarioRes<<endl;


	  for (int ii=0; ii<No_instance; ii++){
		  IloEnv env;
		  IloNum rho=config[ii];
		  
	  
/////////////////// DECISION VARIABLES WITH NAMES  /////////////////////////////
	        
	  ////////DECISION VARIABLES AND PARAMETERS FOR MASTER PROBLEM///////////

	   IloInt nCut = 0;
	   char varName[100];

	   IloNumArray2 pi_capacity(env, numScen);
	   for(w=0; w<numScen; w++) {
			pi_capacity[w] = IloNumArray(env, noGen);
			for(i=0;i<noGen;i++){
				pi_capacity[w][i]=0;
			}	
	   }

	   
	   
	   IloNumArray2 pi_demand(env, numScen);
	   for(w=0; w<numScen; w++) {
		   pi_demand[w]=IloNumArray(env, noDem); //dual price for demand constraints
		   for(j=0; j<noDem; j++) {
			   pi_demand[w][j]=0;
		   }
	   }

  	   IloNumVarArray  X(env, noGen, 0, IloInfinity);
	   for(i=0; i<noGen; i++){ 
		   sprintf_s(varName, "X_%d", (int) i);
		   X[i].setName(varName);
	   }

	   
  	   IloNumVarArray	Theta(env, numScen, 0, IloInfinity);
	   for(w=0; w<numScen; w++){ 
		   sprintf_s(varName, "Theta_%d",(int) w);
		   Theta[w].setName(varName);
	   }	
	  
	   IloNumVar Alpha(env, 0, IloInfinity);
	   sprintf_s(varName, "Alpha");
	   Alpha.setName(varName);
	   
	   IloNumArray Best_X(env, noGen); 

	   


	 ////////DECISION VARIABLES AND PARAMETERS FOR SUBPROBLEM///////////////

 	   IloNumArray x_hat(env, noGen);     

		IloNumVarArray Pi_Capacity(env, noGen, 0, IloInfinity);//dual variable for capacity constraints
		for(i=0;i<noGen;i++){
			sprintf_s(varName, "Pi_Capacity_%d",(int) i);
			Pi_Capacity[i].setName(varName);
		 }
	  


	   IloNumVarArray Pi_Demand(env, noDem, 0, IloInfinity);  //dual variable for demand constraints
	   for(j=0; j<noDem; j++){ 
		   sprintf_s(varName, "Pi_Demand_%d",(int) j);
		   Pi_Demand[j].setName(varName);
	   }	

	     

	   ////////DECISION VARIABLES AND PARAMETERS FOR Row-generation Subproblem///////////////

	   IloNumVarArray P(env, numScen, 0, IloInfinity);
	   for(w=0;w<numScen;w++){
		   sprintf_s(varName, "P_%d", (int) w);
		   P[w].setName(varName);
	   }

	   IloNumVarArray Z(env, numScen, 0, IloInfinity);
	   for(w=0;w<numScen;w++){
		   sprintf_s(varName, "Z_%d", (int) w);
		   Z[w].setName(varName);
	   }

  
	
	   IloNumArray p(env, numScen);
	   for(w=0; w<numScen; w++){
		   p[w]=prob[w];
	   }

	   ////////PARAMETER to do Post-processing///////////////

	
	   Scenario *Output = new Scenario[numScen];
	   

	   //Scenario Output[numScen];

	   for(w=0; w<numScen; w++){
			Output[w].No=w;
		   Output[w].Status=NULL;   
		   Output[w].new_ObjCost=-1;
		   Output[w].partial_sum=0;
		   Output[w].nom_Prob=prob[w];
		   Output[w].Realization=new int [5];
		   for (ww=0; ww<5; ww++){
			   Output[w].Realization[ww]=realization[w][ww];
		   } 
		   
		   
	   }
	   

////////////BUILD MASTER AND SUBPROBLEM MODELS  //////////////////////////

	    ////// MASTER PROBLEM///////////////////////////
		cout<<"building master and subproblems"<<endl;

	    IloModel model_master(env);
		//Objective function
  	    IloExpr Objective_master(env);
	    Objective_master += Alpha + IloScalProd(cost, X);	
        model_master.add(IloMinimize(env, Objective_master));
        
		//Constraints
		for(i=0; i<noGen; i++){
			model_master.add(X[i] <= ccmax[i]);  //capacity reequired
			model_master.add(X[i] >= ccmin[i]);  
		}

		//Construct model
	    IloCplex cplex_master(model_master);
		

        //////SUBPROBLEM (DUAL FORMULATION)//////////////////////

	    IloModel model_sub(env);
		//Objective function
	    IloObjective Objective_sub(env);
		Objective_sub.setSense(IloObjective::Maximize); 
 	    model_sub.add(Objective_sub);


        //Constraints
 	    for(j=0; j<noDem; j++){
			for(i=0; i<noGen; i++){ 
			    model_sub.add(-Pi_Capacity[i] + Pi_Demand[j] <= f[i][j]);
		    }
		    model_sub.add(Pi_Demand[j] <= us[j]);
	     } 

		 //Construct model
	     IloCplex cplex_sub(model_sub);    

		  //////ROW_GENERATION SUBPROBLEM//////////////////////

		 IloModel  model_row(env);
		 IloObjective Objective_row(env);
		 Objective_row.setSense(IloObjective::Maximize); 
		 model_row.add(Objective_row);
		 
		 IloRangeArray conDistancePos(env);
		 IloRangeArray conDistanceNeg(env);
		 
		 for(w=0;w<numScen; w++){ 
			 conDistancePos.add(P[w]-Z[w]<=prob[w]);
			 conDistanceNeg.add(-P[w]-Z[w]<=-prob[w]);
		 }

		 model_row.add(conDistancePos);
		 model_row.add(conDistanceNeg);
	
		 IloExpr conDistance(env);
		 IloRange RngDistance;
		 IloExpr conProb(env);
		 IloRange RngProb;
		 
		 for(w=0;w<numScen; w++){ 
			 conDistance+=Z[w];
			 conProb+=P[w];
		 }
		 
		 RngDistance=conDistance<=2*rho;
		 model_row.add(RngDistance);
		 RngProb=conProb==1;
		 model_row.add(RngProb);

		 
		 //Constrcut model
		 IloCplex cplex_row(model_row);


//////////BEGIN ITERATIONS/////////////////////////////////

        cout<<"begin iteration..."<<endl;
		cout.precision(30);
		//result file name
		char resName[100];
		//Initialize LB, UB
		IloNum rel_Gap = IloInfinity;
		IloNum LB = 0;
		IloNum UB = IloInfinity;
		IloNum z_hat;
		 //Scenario cost
		IloNumArray subObj_hat(env, numScen);
		
		sprintf_s(resName, "EX_%.2f.txt", rho);
		const char* Resfilename  = resName;
		ofstream Result(Resfilename);

		std::cout.setf( std::ios::fixed, std:: ios::floatfield );
		Result.precision(10);

		Result<<"Number of scenarios:"<< numScen<<endl;
		Result<<endl;
	
	Result<< "Solve the original problem.."<<endl;
	Result <<"Iter"<<'\t'<<"Z_hat"<<'\t'<<"LB"<<'\t'<<"UB"<<'\t'<<"Relative_Gap"<<'\t'<<"Incumbent"<<endl;
		

		

		while (rel_Gap > toler) {

			
			
			

			if( !cplex_master.solve () ){
				env.error() << "Master problem infeasible" << endl;
				throw(-1);
			}  
  
			cplex_master.getValues(x_hat, X); 
			LB = cplex_master.getObjValue();


			IloExpr MaxCut(env);
			
			//root node objective function
			z_hat = 0;

			//Solve subproblem dual for each scenario
			for (w=0; w<numScen; w++){
				IloExpr sub_obj(env);
				for (i=0; i<noGen; i++){
					sub_obj+=-x_hat[i]*Pi_Capacity[i]*alpha[i][realization[w][i]];
				}
				for (j=0; j<noDem; j++){
					sub_obj+=Pi_Demand[j]*demand[j][realization[w][j+noGen]];
				}
				Objective_sub.setExpr(sub_obj);
				sub_obj.end(); 
				// Solve the scenario subproblem
				cplex_sub.solve(); 
						
				subObj_hat[w] = cplex_sub.getObjValue(); 	
						
				cplex_sub.getValues(pi_capacity[w], Pi_Capacity);
				cplex_sub.getValues(pi_demand[w], Pi_Demand);
						
				IloExpr Gx(env);
				for (i=0; i<noGen; i++){
					Gx+=-X[i]*pi_capacity[w][i]*alpha[i][realization[w][i]];
				}
				for (j=0; j<noDem; j++){
					Gx+=pi_demand[w][j]*demand[j][realization[w][j+noGen]];
				}
				model_master.add(Theta[w] >= Gx); 
							
				Gx.end();
				
			}	
				

			//Solve row-generion subproblem and calculationg p
			
			for(w=0;w<numScen;w++){	
				Objective_row.setLinearCoef(P[w],subObj_hat[w]);
			}
			
			for(w=0;w<numScen;w++){		
				conDistancePos[w].setUb(prob[w]);
				conDistanceNeg[w].setUb(-prob[w]);
			}
			//cplex_row.exportModel("prob.lp");
			cplex_row.solve();

			
			cplex_row.getValues(p, P);
					
			

			MaxCut=IloScalProd(Theta, p);
			model_master.add( Alpha >= MaxCut);

			//*****Calculate gap****////
			z_hat+=IloScalProd(subObj_hat, p);		
			z_hat += IloScalProd(cost, x_hat); 
			if (z_hat<= UB || abs(z_hat-UB)/abs(UB)<= toler){
				//store results
				UB = z_hat;
				for(i=0; i<noGen; i++){
					Best_X[i] = x_hat[i]; 
				}

				for(w=0; w<numScen; w++){
					Output[w].Cost=subObj_hat[w];
					Output[w].worst_Prob=p[w];
				}
				
				
				
			}
			
			
			if( (UB - LB)/abs(LB) < rel_Gap){
				rel_Gap = (UB - LB)/abs(LB);
			}
			Result<< nCut+1 << '\t' << z_hat << '\t' << LB <<'\t' << UB <<'\t' << rel_Gap<< '\t' << Best_X <<endl;
			nCut++;
			
		} //end while   

		Details << rho << '\t'<<Best_X << '\t';
		

		Result<<endl;
		IloNum Obj=LB;
		model_master.end();
		cplex_master.end();
		cout<<"DONE"<< endl;

		
		//sort 

		struct CompareScenarios
		{
			bool operator()(const Scenario& lhs, const Scenario& rhs)  const
			{
				if (lhs.Cost == rhs.Cost)
					return lhs.nom_Prob > rhs.nom_Prob;
				else
					return lhs.Cost < rhs.Cost;
			}
		};

		//STEP 0
		//count number of times each realization of demand and availablity repeated for Effective scenarios
		//0: Gen 1; 1: Gen 2; 2: Dem 1; 3: Dem 2; 4: Dem 3
		IloIntArray category(env, 4);
		IloIntArray category_effect(env, 3); //0: Ineffective; 1: Effective; 2: Unknown
		for (i=0; i<3; i++){
			category_effect[i]=0;
		}

		std::sort(Output, Output+numScen, CompareScenarios());
		
		Output[0].partial_sum=Output[0].nom_Prob;
		for(w=1; w<numScen; w++){
			Output[w].partial_sum=Output[w-1].partial_sum+Output[w].nom_Prob;
		}

		vector<int>::iterator it;

		//STEP 1
		vector<int> Sup_index;
		Sup_index.push_back(numScen-1);

		w=numScen-2;
		while (abs(Output[w].Cost-Output[numScen-1].Cost)<toler){
			Sup_index.push_back(w);
			w--;
		}
		const int noSup=Sup_index.size();		
		category[3]=noSup;


		//STEP 2 
		vector<int> VaR_index;
		if(Output[0].partial_sum >= rho )
			VaR_index.push_back(0);

		if (VaR_index.size()==0){
			for (w=1; w<numScen+1; w++){
				if(Output[w].partial_sum>=rho-toler &&  Output[w-1].partial_sum < rho){
					VaR_index.push_back(w);
					break;
				}
			}
		}
		

		
		
		
		//check remval of scenarios with the same cost as index scenarios
		
		vector<int> removal;
		if (VaR_index[0]!=0){
			w=VaR_index[0]-1;
			while (abs(Output[VaR_index[0]].Cost-Output[w].Cost)<toler && w>=0 && w<numScen){
				removal.push_back(w);
				w--;
			}
		}
		if (VaR_index[0]!=0)
			category[0]=w+1;
		else
			category[0]=0;
		if (category[0]!=0){
			for (int ww=0; ww<=w; ww++){
				Ineffscen[ii].push_back(ww);
				Output[ww].Status="Ineffective";
				category_effect[0]++;
			}
		}

		w=VaR_index[0];
		while (abs(Output[VaR_index[0]].Cost-Output[w].Cost)<toler && w>=0 && w<numScen){
			removal.push_back(w);
			w++;
		}
		const int noRemoval=removal.size();
		category[1]=noRemoval;

		//lambda>0
		if (abs(Output[VaR_index[0]].Cost-Output[Sup_index[0]].Cost)> toler){
			for (it=Sup_index.begin(); it!=Sup_index.end(); it++){
				if (Output[*it].nom_Prob <= rho){
					if (Output[*it].nom_Prob>0 && noSup>1){
						category_effect[1]++;// EFFECTIVE
						Effscen[ii].push_back(*it);
						Output[*it].Status="Effective";
					}
					else if (noSup==1){
						category_effect[1]++; // EFFECTIVE
						Effscen[ii].push_back(*it);
						Output[*it].Status="Effective";
					}
					else{
						Unknownscen[ii].push_back(*it);
						category_effect[2]++; //UNDETERMINED
					}
				}
				else{
					category_effect[1]++; // EFFECTIVE
					Output[*it].Status="Effective";
					Effscen[ii].push_back(*it);
				}
			}
		}

		Details << "[";
		for (it = removal.begin(); it != removal.end(); it++) {
			if (++it == removal.end()) {
				it--;
				Details << Output[*it].No << "]" << '\t';
			}
			else {
				it--;
				Details << Output[*it].No << ", ";
			}
		}

		if (abs(Output[VaR_index[0]].Cost - Output[Sup_index[0]].Cost) > toler) {
			Details << "[";
			for (it = Sup_index.begin(); it != Sup_index.end(); it++) {
				if (++it == Sup_index.end()) {
					it--;
					Details << Output[*it].No << "]" << endl;
				}
				else {
					it--;
					Details << Output[*it].No << ", ";
				}
			}
		}
		else
			Details << endl;


		//lambda>0
		if (abs(Output[VaR_index[0]].Cost-Output[Sup_index[0]].Cost)> toler){
			for (it=removal.begin(); it!=removal.end(); it++){
				if (Output[*it].nom_Prob <= rho){
					if (Output[*it].nom_Prob>0 && noRemoval>1){
						category_effect[2]++; //UNDETERMINED
						Unknownscen[ii].push_back(*it);
					}
					else if (noRemoval==1 && Output[*it].worst_Prob==0){
						category_effect[0]++; //INEFFECTIVE
						Output[*it].Status="Ineffective";
						Ineffscen[ii].push_back(*it);
					}
					else if (noRemoval==1 && Output[*it].worst_Prob>0){
						category_effect[1]++; // EFFECTIVE
						Output[*it].Status="Effective";
						Effscen[ii].push_back(*it);
					}
					else if (Output[*it].nom_Prob==0 && noRemoval>1){
						category_effect[0]++; //INEFFECTIVE
						Output[*it].Status="Ineffective";
						Ineffscen[ii].push_back(*it);
					}
				}
				else{
					category_effect[1]++; //EFFECTIVE
					Effscen[ii].push_back(*it);
					Output[*it].Status="Effective";
				}
			}
		}

				//lambda=0
		if (abs(Output[VaR_index[0]].Cost-Output[Sup_index[0]].Cost)< toler){
			for (it=removal.begin(); it!=removal.end(); it++){
				if (Output[*it].nom_Prob <= rho){
					if (noRemoval==1){
						category_effect[1]++; //EFFECTIVE
						Output[*it].Status="Effective";
						Effscen[ii].push_back(*it);
					}
					//else if (Output[*it].nom_Prob>0 && noRemoval>1){
					//	category_effect[1]++; //EFFECTIVE
					//	Output[*it].Status="Effective";
					//	Effscen[ii].push_back(*it);
					//	for (i=0; i<noGen+noDem; i++){
					//		freq_effective[i][Output[*it].Realization[i]]++;
					//	}
					//}
					//else if (Output[*it].nom_Prob==0 && noRemoval>1){
					//	category_effect[2]++; //UNDETERMINED
					//	Unknownscen[ii].push_back(*it);
					//}
					else{
						category_effect[2]++; //UNDETERMINED
						Unknownscen[ii].push_back(*it);
					}
				}
				else{
					category_effect[1]++; //EFFECTIVE
					Output[*it].Status="Effective";
					Effscen[ii].push_back(*it);
				}
			}
		}
	

		//lambda=0
		if (abs(Output[VaR_index[0]].Cost-Output[Sup_index[0]].Cost)< toler){
			category[2]=0;
			category[3]=0;
		}
		else{
			category[2]=numScen-(category[0]+category[1]+category[3]);
		}


		vector<int> third;

		//lambda>0
		if (abs(Output[VaR_index[0]].Cost-Output[Sup_index[0]].Cost)> toler){
			while (abs(Output[w].Cost - Output[Sup_index[0]].Cost)>toler && w>=0 && w<numScen){
				third.push_back(w);
				w++;
			}
			for (it=third.begin(); it!=third.end(); it++){
				if (Output[*it].nom_Prob <= rho){
					if (Output[*it].nom_Prob>0){
						category_effect[1]++; //EFFECTIVE
						Output[*it].Status="Effective";
						Effscen[ii].push_back(*it);
					}
					else {
						category_effect[0]++; //INEFFECTIVE
						Output[*it].Status="Ineffective";
						Ineffscen[ii].push_back(*it);
					}
				}
				else{
					category_effect[1]++; //EFFECTIVE
					Output[*it].Status="Effective";
					Effscen[ii].push_back(*it);
				}
				
			}
		}

		Result<< "Information about each scenario.."<<endl;
		Result << "w"<<'\t'<<"xi"<<'\t'<<"h"<<'\t'<<"q"<<'\t'<<"p"<<'\t'<<"cdf"<<'\t'<<"Status"<<endl;
		for(w=0; w<numScen; w++){
			Result<< Output[w].No<< '\t';
			Result<< "[";
			for (int ww=0; ww<4; ww++){
				Result<<Output[w].Realization[ww]<< ", ";
			}
			Result<<Output[w].Realization[4]<< "]"<<'\t';
			Result<<Output[w].Cost<< '\t'<< Output[w].nom_Prob<< '\t' << Output[w].worst_Prob << '\t' << Output[w].partial_sum <<'\t';
			if (Output[w].Status !=NULL){
				Count_Known[Output[w].No]++;
				Result<< Output[w].Status<<endl;	
			}
			else {
				Result<< endl;
			}

		}

		//IloNum temp_prob=0; 
		//////check the the objective function after an effective scenario is removed
		//for (it=Effscen[ii].begin(); it!=Effscen[ii].end(); it++){
		//	vector<int> temp_VaR_index;
		//	//ww=Output[*it].No;
		//	temp_prob=Output[*it].nom_Prob;
		//	IloNumArray partial_sum(env, numScen);
		//	if (rho==0 && *it!=VaR_index[0]){
		//		temp_VaR_index.push_back(VaR_index[0]);
		//	}
		//	else if (*it==VaR_index[0]){
		//		temp_VaR_index.push_back(2);
		//	}
		//	
		//	if (rho>0){
		//		if (*it!=0){
		//			partial_sum[0]=Output[0].nom_Prob;
		//		}
		//		else{
		//			partial_sum[0]=0;
		//		}
		//		for(w=1; w<numScen; w++){
		//			if (*it!=w){
		//				partial_sum[w]=partial_sum[w-1]+Output[w].nom_Prob;
		//			}
		//			else{
		//				partial_sum[w]=partial_sum[w-1];
		//			}
		//		}

		//		if(partial_sum[0] >= rho/2-temp_prob && *it !=0)
		//			temp_VaR_index.push_back(0);

		//		if (temp_VaR_index.size()==0){
		//			for (w=1; w<numScen+1; w++){
		//				if(partial_sum[w] >= rho/2-temp_prob &&  partial_sum[w-1] < rho/2-temp_prob && *it!=w){
		//					temp_VaR_index.push_back(w);
		//					break;
		//				}
		//			}
		//		}
		//	}

		//	vector<int>::iterator match;
		//	IloNum newcost;
		//	match=std::find(Sup_index.begin(), Sup_index.end(), *it);
		//	if (match ==Sup_index.end()){
		//		newcost=rho/2*Output[numScen-1].Cost;
		//	}
		//	else if (Sup_index.size()>1){
		//		newcost=rho/2*Output[numScen-1].Cost;
		//	}
		//	else{
		//		newcost=rho/2*Output[numScen-2].Cost;
		//	}
		//	newcost+=(1-rho/2)*Output[temp_VaR_index[0]].Cost;
		//	for(w=temp_VaR_index[0]+1; w<numScen; w++){
		//		if (*it !=w){
		//			newcost+=Output[w].nom_Prob*(Output[w].Cost-Output[temp_VaR_index[0]].Cost);
		//		}
		//	}
		//	//store ratio of of the change in optimal value after removal of effective scenario
		//	NewCost_Effscen[ii].push_back(100*(LB-newcost)/LB);
		//
		//	
		//}


		Result<<endl;
		for (i=0; i<3; i++){
			Result<< category_effect[i]<< "& ";
		}
		Result<<endl;
		Result<<endl;

		Numbers<< rho<<'\t';
		for (i=0; i<3; i++){
			Numbers<< category_effect[i]<< '\t';
		}
		Numbers<<endl;


		//
		Result<<"Solve the new problem for unknown scenarios.."<<endl;
		//check the status of undetermined scenarios
		for (it= Unknownscen[ii].begin(); it!= Unknownscen[ii].end(); it++){
			
			ww=Output[*it].No;
			Result<< "w="<< ww <<endl;		
			
			if (Output[*it].nom_Prob<=rho){
				SubScenario * SubOutput= new SubScenario[numScen];

				
				IloModel MODEL(env);
				IloCplex CPX(env);
				IloExpr Objective_master_removal(env);
				Objective_master_removal += Alpha + IloScalProd(cost, X);	
				MODEL.add(IloMinimize(env, Objective_master_removal));
        
				//Constraints
				for(i=0; i<noGen; i++){
					MODEL.add(X[i] <= ccmax[i]);  //capacity reequired
					MODEL.add(X[i] >= ccmin[i]);  
				}
				CPX.extract(MODEL); 
			   for(w=0; w<numScen; w++){
				   SubOutput[w].nom_Prob=prob[w];
			   }
				IloRange RngZero;
				RngZero=P[ww]==0;
				model_row.add(RngZero);

				rel_Gap = IloInfinity;
				LB = -IloInfinity;
				UB = IloInfinity;
				//Scenario cost
				for(w=0; w<numScen; w++){
					subObj_hat[w]=0;
				}
				z_hat=0;
			
				IloBool ctnWhile=IloTrue;
				while (rel_Gap > toler) {
					if( !CPX.solve () ){
						env.error() << "Master problem infeasible" << endl;
						throw(-1);
					}  
					CPX.getValues(x_hat, X); 
					LB = CPX.getObjValue();


					IloExpr MaxCut(env);
			
					//root node objective function
					z_hat = 0;

					//Solve subproblem dual for each scenario
					for (w=0; w<numScen; w++){
						IloExpr sub_obj(env);
						for (i=0; i<noGen; i++){
							sub_obj+=-x_hat[i]*Pi_Capacity[i]*alpha[i][realization[w][i]];
						}
						for (j=0; j<noDem; j++){
							sub_obj+=Pi_Demand[j]*demand[j][realization[w][j+noGen]];
						}
						Objective_sub.setExpr(sub_obj);
						sub_obj.end(); 
						// Solve the scenario subproblem
						cplex_sub.solve(); 
						
						subObj_hat[w] = cplex_sub.getObjValue(); 	
						
						cplex_sub.getValues(pi_capacity[w], Pi_Capacity);
						cplex_sub.getValues(pi_demand[w], Pi_Demand);
						
						IloExpr Gx(env);
						for (i=0; i<noGen; i++){
							Gx+=-X[i]*pi_capacity[w][i]*alpha[i][realization[w][i]];
						}
						for (j=0; j<noDem; j++){
							Gx+=pi_demand[w][j]*demand[j][realization[w][j+noGen]];
						}
						MODEL.add(Theta[w] >= Gx); 
							
						Gx.end();
				
					}
	
				

					//Solve row-generion subproblem and calculationg p
			
					for(w=0;w<numScen;w++){	
						Objective_row.setLinearCoef(P[w],subObj_hat[w]);
					}
			
					for(w=0;w<numScen;w++){		
						conDistancePos[w].setUb(prob[w]);
						conDistanceNeg[w].setUb(-prob[w]);
					}
					//cplex_row.exportModel("prob.lp");
					cplex_row.solve();

					if( cplex_row.solve () ){ 
						cplex_row.getValues(p, P);
					
			

						MaxCut=IloScalProd(Theta, p);
						MODEL.add( Alpha >= MaxCut);

						//*****Calculate gap****////
						z_hat+=IloScalProd(subObj_hat, p);		
						z_hat += IloScalProd(cost, x_hat); 
						if (z_hat<= UB || abs(z_hat-UB)/abs(UB)<= toler){
							//store results
							UB = z_hat;
							for(i=0; i<noGen; i++){
								Best_X[i] = x_hat[i]; 
							}

							for(w=0; w<numScen; w++){
								SubOutput[w].WorstCost=subObj_hat[w];
							}
				
				
				
						}
			
			
						if( (UB - LB)/abs(LB) < rel_Gap){
							rel_Gap = (UB - LB)/abs(LB);
						}
							Result<< z_hat << '\t' << LB <<'\t' << UB <<'\t' << rel_Gap<< '\t' << Best_X  <<endl;
						}
						else{
							rel_Gap=toler;
							ctnWhile=IloFalse;

						}
			
					} //end while   

					IloNum new_Obj=LB;
					Output[*it].new_ObjCost=new_Obj;

					if (ctnWhile) {
						if (abs(Output[*it].new_ObjCost - Obj)>toler) {
							Output[*it].Status = "Effective";
							category_effect[1]++;
							Effscen[ii].push_back(*it);
						}
						else {
							category_effect[0]++;
							Output[*it].Status = "Ineffective";
							Ineffscen[ii].push_back(*it);
						}
					}
					else {
						category_effect[1];
						Output[*it].Status = "Effective";
						Effscen[ii].push_back(*it);
					}

					//cout<<100*(Obj-new_Obj)/Obj<<endl;
					NewCost_Effscen[ii].push_back(100*(Obj-new_Obj)/Obj);
					model_row.remove(RngZero);
					delete[] SubOutput;
					MODEL.end();
					CPX.end();
					}
					else {
						category_effect[1];
						Output[*it].Status = "Effective";
						Effscen[ii].push_back(*it);
					}
					Result << endl;

				}//end for ww

				Result<<"Information about (primarily) unknown scenarios.."<<endl;
				vector<double>::iterator itt;
				Result << "w"<<'\t'<<"f_n(x*)"<<'\t'<<"f(bar{x})"<<'\t'<<"Status"<<endl;

				for (it=Unknownscen[ii].begin(); it!=Unknownscen[ii].end(); it++){
					Result<< Output[*it].No<< '\t';
					Result<< "[";
					for (ww=0; ww<4; ww++){
						Result<<Output[*it].Realization[ww]<< ", ";
					}
					Result<<Output[*it].Realization[4]<< "]"<<'\t';
					Result<< Output[*it].new_ObjCost<< '\t'<< Obj<< '\t'<<Output[*it].Status<< endl;
				}
		
				Result<<endl;

				Result<<"Number of Ineffective & Effective"<<endl;
				for (i=0; i<2; i++){
					Result<< category_effect[i]<< "& ";
				}
			
				Result<<endl;
				Result<<endl;

		

		
		cout<<"DONE"<<endl;

		struct CompareScenarios2
		{
			bool operator()(const Scenario& lhs, const Scenario& rhs)  const
			{
				return lhs.No < rhs.No;
			}
		};

		

		std::sort(Output, Output+numScen, CompareScenarios2());

		Effscen[ii].clear();
		Ineffscen[ii].clear();
		
		//Summary << rho << '\t';
		for (w = 0; w<numScen; w++) {
			if (Output[w].Status == "Effective") {
				Effscen[ii].push_back(w);
				//Summary << "E" << '\t';
			}
			else if (Output[w].Status == "Ineffective") {
				Ineffscen[ii].push_back(w);
				//Summary << "I" << '\t';
			}
		}
		//Summary << endl;
		

		for (w=0; w<numScen; w++){
			if (Output[w].Status=="Effective"){
				Count_Effective[w]++;
			}
		}


		Numbers_post<< rho<<'\t';
		for (i=0; i<2; i++){
			Numbers_post<< category_effect[i]<< '\t';
		}
		Numbers_post<<endl;

		//			
		
		Result.close();
		model_sub.end();
		cplex_sub.end();
		model_row.end();
		cplex_row.end();
		env.end();

}//end instances loop

		for (w=0; w<numScen; w++){
			ScenarioRes<< w<< '\t';
			for (ww=0; ww<2; ww++){
				ScenarioRes<<alpha[ww][realization[w][ww]]<< '\t';
			}
			for (ww=2; ww<4; ww++){
				ScenarioRes<<demand[ww-2][realization[w][ww]]<< '\t';
			}
			ScenarioRes<<demand[2][realization[w][ww]]<<'\t';
			ScenarioRes<<Count_Effective[w]<<'\t'<< Count_Known[w]<<endl;


		}
		/*
		for (w=0; w<numScen; w++){
			for (i=0; i<No_instance; i++){
				if (NewCost_Effscen[i].size()>=w+1){
					EffScenarioRes<<NewCost_Effscen[i][w]<<'\t';
				}
				else{
					EffScenarioRes<<"0"<<'\t';
				}
			}
			EffScenarioRes<<endl;
		}
		*/

		//do some works here to check nestedness
		const char* Nestfilename = "EX_nestedness.txt";
		ofstream Nested(Nestfilename);

		int ii = 0;
		vector<int>::iterator it;
		vector<int>::iterator match; //if found the item in the next ii
		int nestedness = 1;


		Nested << "Check Nestedness of Effective Scenarios" << endl;
		for (ii = 20; ii>0; ii--) {
			if (Effscen[ii].size() > 0) {
				nestedness = 1;
				for (it = Effscen[ii].begin(); it != Effscen[ii].end(); it++) {
					match = std::find(Effscen[ii - 1].begin(), Effscen[ii - 1].end(), *it);
					if (match == Effscen[ii - 1].end()) {
						Nested << "Found " << *it << " with probability " << prob[*it] << " in gamma " << (double)ii / 10 << " not exits in gamma " << (double)(ii - 1) / 10 << endl;
						if (nestedness == 1)
							nestedness = 0;
						//break;
					}
				}
				if (nestedness == 1) {
					Nested << "gamma " << (double)ii / 10 << " is nested in gamma " << (double)(ii - 1) / 10 << endl;
				}

				if (Effscen[ii].size()>Effscen[ii - 1].size()) {
					Nested << "********** gamma " << (double)ii / 10 << " is larger than " << (double)(ii - 1) / 10 << endl;
				}

			}
			/*
			else if (Effscen[ii].size()>Effscen[ii - 1].size()) {
			Nested << "********** gamma " << ii << " is larger than " << ii - 1 << endl;
			}
			*/
		}


		Nested.close();


		cout << "DONE" << endl;

}//end try



 catch(IloException &e) {
	envData.out() << "ERROR: " << e << endl;
 }

 catch(...){
	envData.out() << "Unknown exception" << endl;
 }
 envData.end();
 getchar();
 return 0;

}
