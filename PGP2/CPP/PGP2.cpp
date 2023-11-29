/*
******************************************
* Author:   Hamed Rahimian
*          The Ohio State University
*          January 01, 2015
******************************************
*/

// -------------------------------------------------------------- -*- C++ -*-
//PGP2


//Input data:
//c(j) Annualized capital cost of installing a unit of generator type j ($ per MW)
//f(j) Cost of producing a unit of energy per hour of generator type j ($ per MW per hour)
//Beta(i) Duration of period i (hour)
//d(k,i) Distributional load during period i in point k
//prob(k,i) Probability of having load during period i in point k
//b Budget to install generators ($)
//M Minimum electricity generation capacity to be installed (MW)
//p Unit subcontarcting cost for load during period i ($ per MW)

// Modelling variables:
//X(j) Capacity installed to generator type i (MW)
//Y(l1,l2,l3,i,j) Units of load during period i which is provided by generator type j under realization w (MW)
//S(l1,l2,l3,i)  Units of subcontracted electricity for load during period i under realization w (MW)


#include <ilcplex/ilocplex.h>
#include <vector>
#include <array>
#include <algorithm> 


ILOSTLBEGIN

typedef IloArray<IloNumVarArray> VarArray2;
typedef IloArray<IloNumArray2> NumArray3Matrix;

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
	low=0,
	med,
	high };

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
	 IloInt w, w1, w2, w3, ww; 

	 IloInt budget=0, capacity=0;	
	 
	 IloInt noPeriod, noType; 
	 
	 IloNumArray cost(envData), f(envData), beta(envData), purchase(envData);
	 IloNumArray2 nom_prob(envData), demand(envData);

///////////////// DATA FILE READING ////////////////////////////////

	  const char* filename  = "inputs.dat";	      

	  ifstream file(filename);

	  if (!file) {
			cerr << "ERROR: can not open file '" << filename << "' for reading" << endl;
			cerr << "usage:   " << argv[0] << " <file>" << endl;
			throw(-1);
		}

      file >> budget >> capacity >> cost >> f >> beta >> purchase; 
      file >> nom_prob >> demand; 

	  noPeriod = purchase.getSize();
	  noType = f.getSize();
	  

	  IloIntArray noScen(envData, noPeriod);
	  
	  for(i=0; i<noPeriod; i++) {
		  noScen[i] = nom_prob[i].getSize(); 
	  }

	  IloBool consistentData = (beta.getSize() == noPeriod && cost.getSize() == noType && 
 	      demand[0].getSize() == noScen[0] && demand[1].getSize() == noScen[1] && demand[2].getSize() == noScen[2]);

      if (!consistentData) { 
 	      cerr << "ERROR: Inconsistent data!" << endl;
          throw(-1);
      }
	  
	  file.close();
	  
	  const IloInt numScen=noScen[0]*noScen[1]*noScen[2];
	  IloNumArray prob(envData, numScen);
	  vector<array<int, 3>> realization;
	  w=0;
	  for(w1=0; w1<noScen[0]; w1++) {
		  for(w2=0; w2<noScen[1]; w2++){
			  for(w3=0; w3<noScen[2]; w3++) {	
				  prob[w]=nom_prob[0][w1]*nom_prob[1][w2]*nom_prob[2][w3];
				  array<int, 3> temp={w1, w2, w3};
				  realization.push_back(temp);
				  w++;
			  }
		  }
	  }

	  //vector to stores each gamma ineffective scenarios
	  vector<int> * Ineffscen= new vector<int> [21];

	  //vector to stores each gamma effective scenarios
	  vector<int> * Effscen= new vector<int> [21];

	  //vector to stores each gamma undetermined scenarios
	  vector<int> * Unknownscen= new vector<int> [21];

	  //vector to stores each gamma undetermined scenarios based on Type 1 analysis
	  Empirical_Analysis * Interpretation= new Empirical_Analysis [21];

	  IloNumArray Count_Effective(envData, numScen);
	  IloNumArray Count_Known(envData, numScen);

	  for (w=0; w<numScen; w++){
		  Count_Effective[w]=0;
		  Count_Known[w]=0;
	  }


	   const char* Numberfilename  = "Numbers.txt";
	  ofstream Numbers(Numberfilename);
	  Numbers<<"Number of scenarios.."<<endl;
	  Numbers<< '\t';
	  Numbers<< "I"<< '\t'<<"E"<<'\t'<<"U"<<endl;

	  const char* Number_post_filename  = "Numbers_post.txt";
	  ofstream Numbers_post(Number_post_filename);
	  Numbers_post<<"Number of scenarios.."<<endl;
	  Numbers_post<< '\t';
	  Numbers_post<< "I"<< '\t'<<"E"<<'\t'<<"U"<<endl;
	 
	  const char* Scenariofilename  = "Average.txt";
	  ofstream ScenarioRes(Scenariofilename);



	  for (int ii=0; ii<No_instance; ii++){
		  IloEnv env;
		  IloNum rho=config[ii];
		  
	  
/////////////////// DECISION VARIABLES WITH NAMES  /////////////////////////////
	        
	  ////////DECISION VARIABLES AND PARAMETERS FOR MASTER PROBLEM///////////

	   IloInt nCut = 0;
	   char varName[100];

	   NumArray3Matrix pi_capacity(env, numScen);
	   for(w=0; w<numScen; w++) {
	   	   pi_capacity[w] = IloNumArray2(env, noPeriod);
		   for(i=0; i<noPeriod; i++) {
			   pi_capacity[w][i] = IloNumArray(env, noType);
			   for(j=0;j<noType;j++){
				   pi_capacity[w][i][j]=0;
			   }	
		   }
	   }

	   
	   
	   IloNumArray2 pi_demand(env, numScen);
	   for(w=0; w<numScen; w++) {
		   pi_demand[w]=IloNumArray(env, noPeriod); //dual price for demand constraints
		   for(i=0; i<noPeriod; i++) {
			   pi_demand[w][i]=0;
		   }
	   }

  	   IloNumVarArray  X(env, noType, 0, IloInfinity);
	   for(j=0; j<noType; j++){ 
		   sprintf_s(varName, "X_%d", (int) j);
		   X[j].setName(varName);
	   }

	   
  	   IloNumVarArray	Theta(env, numScen, 0, IloInfinity);
	   for(w=0; w<numScen; w++){ 
		   sprintf_s(varName, "Theta_%d",(int) w);
		   Theta[w].setName(varName);
	   }	
	  
	   IloNumVar Alpha(env, 0, IloInfinity);
	   sprintf_s(varName, "Alpha");
	   Alpha.setName(varName);
	   
	   IloNumArray Best_X(env, noType); 

	   


	 ////////DECISION VARIABLES AND PARAMETERS FOR SUBPROBLEM///////////////

 	   IloNumArray x_hat(env, noType);     

 	   VarArray2 Pi_Capacity(env, noPeriod);  //dual variable for capacity constraints
 	   for(i=0; i<noPeriod; i++) {
		   Pi_Capacity[i] = IloNumVarArray(env, noType, -IloInfinity, 0);
		   for(j=0;j<noType;j++){
			   sprintf_s(varName, "Pi_Capacity_%d_%d",(int) i, (int) j);
			   Pi_Capacity[i][j].setName(varName);
		   }
	   }


	   IloNumVarArray Pi_Demand(env, noPeriod, -IloInfinity, IloInfinity);  //dual variable for demand constraints
	   for(i=0; i<noPeriod; i++){ 
		   sprintf_s(varName, "Pi_Demand_%d",(int) i);
		   Pi_Demand[i].setName(varName);
	   }	

	   //IloNumArray worst_cost(env, numScen);    

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

	   ////////PARAMETER to do Pre-processing///////////////

	
	   Scenario *Output = new Scenario[numScen];
	   

	   //Scenario Output[numScen];

	   for(w=0; w<numScen; w++){
		   Output[w].No=w;
		   Output[w].Status=NULL;   
		   Output[w].new_ObjCost=-1;
		   Output[w].partial_sum=0;
		   Output[w].nom_Prob=prob[w];
		   Output[w].Realization=new int [3];
		   for (ww=0; ww<3; ww++){
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
        model_master.add(IloScalProd(cost, X) <= budget);  //budget constraints
        model_master.add(IloSum(X) >= capacity); //Capacity lower bound

		//Construct model
	    IloCplex cplex_master(model_master);
		

        //////SUBPROBLEM (DUAL FORMULATION)//////////////////////

	    IloModel model_sub(env);
		//Objective function
	    IloObjective Objective_sub(env);
		Objective_sub.setSense(IloObjective::Maximize); 
 	    model_sub.add(Objective_sub);


        //Constraints
 	    for(i=0; i<noPeriod; i++){ 
		    for(j=0; j<noType; j++){
			    model_sub.add(Pi_Capacity[i][j] + Pi_Demand[i] <= beta[i]*f[j]);
		    }
		    model_sub.add(Pi_Demand[i] <= purchase[i]);
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
		
		sprintf_s(resName, "%.2f.txt", rho);
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
						
				for(i=0; i<noPeriod; i++){ 
					for(j=0; j<noType; j++){			
						sub_obj += x_hat[j] * Pi_Capacity[i][j];
					}
				}
				for(i=0; i<noPeriod; i++){ 
					sub_obj += Pi_Demand[i]* demand[i][realization[w][i]];
				}

						
				Objective_sub.setExpr(sub_obj);
				sub_obj.end(); 

				// Solve the scenario subproblem
				cplex_sub.solve(); 
						
				subObj_hat[w] = cplex_sub.getObjValue(); 	
						
						
				for(i=0; i<noPeriod; i++){
					cplex_sub.getValues(pi_capacity[w][i], Pi_Capacity[i]);
				}

				cplex_sub.getValues(pi_demand[w], Pi_Demand);
						
				IloExpr Gx(env);
						
				for(i=0; i<noPeriod; i++) {
					for(j=0; j<noType; j++){	
						Gx += X[j]*pi_capacity[w][i][j];
					}
				}
				
				for(i=0; i<noPeriod; i++){ 
					Gx += pi_demand[w][i]* demand[i][realization[w][i]];
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
				for(j=0; j<noType; j++){
					Best_X[j] = x_hat[j]; 
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
		Result<<endl;
		IloNum Obj=LB;
		
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

		if (rho!=1){
			if (VaR_index.size()==0){
				for (w=1; w<numScen+1; w++){
					if(Output[w].partial_sum+toler >= rho &&  Output[w-1].partial_sum < rho){
						VaR_index.push_back(w);
						break;
					}
				}
			}
		}
		else{
			VaR_index.push_back(numScen-1);
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
			for (int ww=0; ww<2; ww++){
				Result<<Output[w].Realization[ww]<< ", ";
			}
			Result<<Output[w].Realization[2]<< "]"<<'\t';
			Result<<Output[w].Cost<< '\t'<< Output[w].nom_Prob<< '\t' << Output[w].worst_Prob << '\t' << Output[w].partial_sum <<'\t';
			if (Output[w].Status !=NULL){
				Count_Known[Output[w].No]++;
				Result<< Output[w].Status<<endl;	
			}
			else {
				Result<< endl;
			}

		}

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




		//Result<<"Solve the new problem for unknown scenarios.."<<endl;
		////check the status of undetermined scenarios
		//for (it=Unknownscen[ii].begin(); it!=Unknownscen[ii].end(); it++){
		//	ww=Output[*it].No;
		//	Result<< "w="<< ww <<endl;		
		//	
		//	if (Output[*it].nom_Prob<=rho/2){
		//		SubScenario * SubOutput= new SubScenario[numScen];
		//
		//	
		//		IloModel MODEL(env);
		//		IloCplex CPX(env);
		//		IloExpr Objective_master_removal(env);
		//		Objective_master_removal += Alpha + IloScalProd(cost, X);	
		//		MODEL.add(IloMinimize(env, Objective_master_removal));

		//		MODEL.add(IloScalProd(cost, X) <= budget);  //budget constraints
		//		MODEL.add(IloSum(X) >= capacity); //Capacity lower bound

		//		CPX.extract(MODEL); 
		//
		//	

		//   
		//	   for(w=0; w<numScen; w++){
		//		   SubOutput[w].nom_Prob=prob[w];
		//	   }
		//	   IloRange RngZero;
		//	   RngZero=P[ww]==0;
		//	   model_row.add(RngZero);
		//	   
		//	   rel_Gap = IloInfinity;
		//	   LB = -IloInfinity;
		//	   UB = IloInfinity;
		//	//Scenario cost
		//		for(w=0; w<numScen; w++){
		//			subObj_hat[w]=0;
		//		}
		//		z_hat=0;
		//	
		//		IloBool ctnWhile=IloTrue;
		//		while (rel_Gap > toler) {
		//			if( !CPX.solve () ){
		//				env.error() << "Master problem infeasible" << endl;
		//				throw(-1);
		//			}  
		//			CPX.getValues(x_hat, X); 
		//			LB = CPX.getObjValue();


		//			IloExpr MaxCut(env);
		//	
		//			//root node objective function
		//			z_hat = 0;

		//			//Solve subproblem dual for each scenario
		//			for (w=0; w<numScen; w++){

		//				IloExpr sub_obj(env);
		//				
		//				for(i=0; i<noPeriod; i++){ 
		//					for(j=0; j<noType; j++){			
		//						sub_obj += x_hat[j] * Pi_Capacity[i][j];
		//					}
		//				}
		//				for(i=0; i<noPeriod; i++){ 
		//					sub_obj += Pi_Demand[i]* demand[i][realization[w][i]];
		//				}

		//				
		//				Objective_sub.setExpr(sub_obj);
		//				sub_obj.end(); 

		//				// Solve the scenario subproblem
		//				cplex_sub.solve(); 
		//				
		//				subObj_hat[w] = cplex_sub.getObjValue(); 	
		//				
		//				
		//				for(i=0; i<noPeriod; i++){
		//					cplex_sub.getValues(pi_capacity[w][i], Pi_Capacity[i]);
		//				}

		//				cplex_sub.getValues(pi_demand[w], Pi_Demand);
		//				
		//				IloExpr Gx(env);
		//				
		//				for(i=0; i<noPeriod; i++) {
		//					for(j=0; j<noType; j++){	
		//						Gx += X[j]*pi_capacity[w][i][j];
		//					}
		//				}
		//		
		//				for(i=0; i<noPeriod; i++){ 
		//					Gx += pi_demand[w][i]* demand[i][realization[w][i]];
		//				}

		//				MODEL.add(Theta[w] >= Gx); 
		//					
		//				Gx.end();
		//				w++;
		//			}



		//			//Solve row-generion subproblem and calculationg p
		//	
		//			for(w=0;w<numScen;w++){	
		//				Objective_row.setLinearCoef(P[w],subObj_hat[w]);
		//			}
		//	
		//			for(w=0;w<numScen;w++){		
		//				conDistancePos[w].setUb(prob[w]);
		//				conDistanceNeg[w].setUb(-prob[w]);
		//			}
		//			cplex_row.solve();

		//			if( cplex_row.solve () ){ 
		//				cplex_row.getValues(p, P);
		//			
		//	

		//				MaxCut=IloScalProd(Theta, p);
		//				MODEL.add( Alpha >= MaxCut);

		//				//*****Calculate gap****////
		//				z_hat+=IloScalProd(subObj_hat, p);		
		//				z_hat += IloScalProd(cost, x_hat); 
		//				if (z_hat<= UB || abs(z_hat-UB)/abs(UB)<= toler){
		//					//store results
		//					UB = z_hat;
		//					for(j=0; j<noType; j++){
		//						Best_X[j] = x_hat[j]; 
		//					}

		//					for(w=0; w<numScen; w++){
		//						SubOutput[w].WorstCost=subObj_hat[w];
		//					}
		//		
		//		
		//		
		//				}
		//	
		//	
		//				if( (UB - LB)/abs(LB) < rel_Gap){
		//					rel_Gap = (UB - LB)/abs(LB);
		//				}
		//				Result<< z_hat << '\t' << LB <<'\t' << UB <<'\t' << rel_Gap<< '\t' << Best_X  <<endl;
		//			}
		//			else{
		//				rel_Gap=toler;
		//				ctnWhile=IloFalse;

		//			}
		//	
		//		} //end while   

		//		IloNum new_Obj=LB;
		//		Output[*it].new_ObjCost=new_Obj;
		//
		//		if (ctnWhile){
		//				if (abs(Output[*it].new_ObjCost- Obj)>toler){
		//					Output[*it].Status="Effective";
		//					category_effect[1]++;
		//					Effscen[ii].push_back(*it);
		//				}
		//				else{
		//					category_effect[0]++;
		//					Output[*it].Status="Ineffective";
		//					Ineffscen[ii].push_back(*it);
		//				}
		//			}
		//			else{
		//				category_effect[1];
		//				Output[*it].Status="Effective";
		//				Effscen[ii].push_back(*it);
		//			}
		//			model_row.remove(RngZero);
		//			delete[] SubOutput;
		//			MODEL.end();
		//			CPX.end();
		//			}
		//			else{
		//				category_effect[1];
		//				Output[*it].Status="Effective";
		//				Effscen[ii].push_back(*it);
		//			}
		//			Result<<endl;
		//		}//end for ww

		//Result<<"Information about (primarily) unknown scenarios.."<<endl;
		//vector<double>::iterator itt;
		//Result << "w"<<'\t'<<"f_n(x*)"<<'\t'<<"f(bar{x})"<<'\t'<<"Status"<<endl;

		//for (it=Unknownscen[ii].begin(); it!=Unknownscen[ii].end(); it++){
		//	Result<< Output[w].No<< '\t';
		//	Result<< "[";
		//	for (int ww=0; ww<2; ww++){
		//		Result<<Output[*it].Realization[ww]<< ", ";
		//	}
		//	Result<<Output[*it].Realization[2]<< "]"<<'\t';

		//	Result<< Output[*it].No<< '\t'<< Output[*it].new_ObjCost<< '\t'<< Obj<< '\t'<<Output[*it].Status<< endl;
		//}

		//
		//Result<<endl;

		//		Result<<"Number of Ineffective & Effective"<<endl;
		//		for (i=0; i<2; i++){
		//			Result<< category_effect[i]<< "& ";
		//		}
		//	
		//		Result<<endl;
		//		Result<<endl;

		//

		//
		cout<<"DONE"<<endl;
		






		struct CompareScenarios2
		{
			bool operator()(const Scenario& lhs, const Scenario& rhs)  const
			{
				return lhs.No < rhs.No;
			}
		};

		

		std::sort(Output, Output+numScen, CompareScenarios2());
		/*Summary<< rho<<'\t';
		for (w=0; w<numScen; w++){
			if (Output[w].Status=="Effective"){
				Count_Effective[w]++;
				Summary<< "E"<< '\t';	
			}
			else if (Output[w].Status=="Ineffective"){
				Summary<< "I"<< '\t';	
			}
			else{
				Summary<< '\t';
			}
		}
		Summary<<endl;*/

		for (w=0; w<numScen; w++){
			if (Output[w].Status=="Effective"){
				Count_Effective[w]++;
			}
		}

		/*Numbers_post<< rho<<'\t';
		for (i=0; i<2; i++){
			Numbers_post<< category_effect[i]<< '\t';
		}
		Numbers_post<<endl;*/

					
		
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
			ScenarioRes<<demand[ww][realization[w][ww]]<< '\t';
		}
		ScenarioRes<<demand[2][realization[w][ww]]<<'\t';
		ScenarioRes<<Count_Effective[w]<<'\t'<< Count_Known[w]<<endl;

	}


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