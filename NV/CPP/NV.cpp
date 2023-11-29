// obj= W(x-d)+ U(d-x)_{+} -Vd
#include <ilcplex/ilocplex.h>
#include <stdlib.h>   
#include <algorithm> 
#include <time.h>
#include <chrono>
#include <vector>
#include <random>
#include <cmath>
using namespace std;

ILOSTLBEGIN



struct Scenario{
	IloNum partial_sum;
	IloInt No;
	IloNum Cost;
	IloNum nom_Prob;
	IloNum worst_Prob;
	IloNum demand;
	IloNum new_ObjCost; //f_new(x^*)
	char * Status;
};


struct SubScenario{
	IloNum WorstCost;
	IloNum nom_Prob;
	IloNum partial_sum;
};


int main(int argc, char **argv) {
	srand ((unsigned)time(NULL));
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
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
	IloInt i, j, k, w, w2, s, ww;
	const IloInt numStage=2;	
	IloInt nCut=0;
	///////////////// DATA ////////////////////////////////
	//cost=[r, c, s, b]
	//IloNumArray cost (envData, 4, 3.0, 2.0, 0.0, 0.0);
	//IloNumArray cost (envData, 4, 4.0, 3.0, 2.0, 2.0);
	
	//W=cost[0]; U=cost[1]; V=cost[2]; 

	
	//Example 2
	IloInt numScen = 3;
	IloNumArray demand(envData, 3, 2, 5, 1);
	IloNumArray prob(envData, 4, 0.0, 0.3, 0.7);
	IloNumArray cost(envData, 3, 2.0, 1.0, 1.0);

	/*
	//Example 3
	IloInt numScen = 4;
	IloNumArray demand(envData, 4, 1, 2, 3, 4);
	IloNumArray prob(envData, 4, 0.0, 0.5, 0.5, 0.0);
	IloNumArray cost(envData, 3, 9.0, 3.0, -1.0);
	*/

	/*
	//Example 4
	IloInt numScen = 6;
	IloNumArray demand(envData, 6, 1, 2, 3, 4, 5, 6);
	IloNumArray prob(envData, 4, 0.0, 0.2, 0.25, 0.2, 0.35, 0.0);
	IloNumArray cost(envData, 3, 9.0, 1.0, -4.0);
	*/

	//Example 
	
	//read modified scenario tree
	/*const char* filename2  = "ScenFac.txt";	      
	 ifstream file2(filename2);
	 file2>> numScen>> prob>>demand;
	 file2.close();
*/
	
	
	//UNIFORM DISTRIBUTION
	/*
	IloInt numScen=5;
	IloNumArray demand (envData);
	IloNumArray prob(envData);
	
	for (w=0; w<numScen; w++){
		demand.add(w+1);
		prob.add((double)1/numScen);
	}
	*/

	//SKEWED DISTRIBUTION
	/*for (w=0; w<numScen; w++){
		demand.add(w+1);
	}
	for (w=0; w<numScen/10; w++){
		prob.add((double)10/numScen);
	}
	for (w=numScen/10; w<numScen; w++){
		prob.add(0);
	}*/

	/*for (w=0; w<numScen; w++){
		demand.add(w+1);
	}
	for (w=0; w<numScen/10; w++){
		prob.add((double)2/numScen);
	}
	for (w=numScen/10; w<numScen/4; w++){
		prob.add(0);
	}
	for (w=numScen/4; w<numScen/2; w++){
		prob.add((double)1/numScen);
	}
	
	for (w=numScen/2; w<4*numScen/5-1; w++){
		prob.add((double)1/numScen);
	}
	prob.add(0.06);
	for (w=4*numScen/5; w<numScen-1; w++){
		prob.add(0);
	}
	prob.add(0.20);*/
	  //vector to stores each gamma ineffective scenarios
	  vector<int> * Ineffscen= new vector<int> [21];

	  //vector to stores each gamma effective scenarios
	  vector<int> * Effscen= new vector<int> [21];

	  //vector to stores each gamma undetermined scenarios
	  vector<int> * Unknownscen= new vector<int> [21];

	  //vector to stores each gamma All scenarios
	  vector<int> * Allscen = new vector<int>[21];
	  
	 

	  IloNumArray Count_Effective(envData, numScen);
	  IloNumArray Count_Known(envData, numScen);

	  for (w=0; w<numScen; w++){
		  Count_Effective[w]=0;
		  Count_Known[w]=0;
	  }



	  const char* Detfilename  = "EX_Details.txt";
	  ofstream Details(Detfilename);
	  
	  Details <<"n= "<<numScen<<endl;
	  Details <<"[W, U, V]= ["<<cost[0]<<", "<<cost[1]<< ", "<< cost[2]<<"]: ";
	  if ((cost[0] + cost[2] > 0) && (cost[1] - cost[2] > 0))
		  Details << "C1" << endl;
	  else if ((cost[0] + cost[2] > 0) && (cost[1] - cost[2] == 0))
		  Details << "C2a" << endl;
	  else if ((cost[0] + cost[2] > 0) && (cost[1] - cost[2] < 0))
		  Details << "C2b" << endl;
	  else if ((cost[0] + cost[2] == 0) && (cost[1] - cost[2] > 0))
		  Details << "C3a" << endl;
	  else if ((cost[0] + cost[2] < 0) && (cost[1] - cost[2] > 0))
		  Details << "C3b" << endl;

	  Details << "demand=" << demand << endl;
	  Details << "probability=" << prob << endl;
	  Details << endl;

	  Details << '\t' << "X*" << '\t' << "VaR" << '\t' << "Sup" << '\t'<< "Eff"<<endl;


	  const char* Statusfilename  = "EX_Status.txt";
	  ofstream Summary(Statusfilename);
	  
	  Summary<<"Status of scenarios.."<<endl;
	  Summary<< '\t';
	  for (w=0; w<numScen-1; w++){
		Summary << w<< '\t';
	  }
	  Summary<<numScen-1<<endl;


	  const char* Numberfilename  = "EX_Numbers.txt";
	  ofstream Numbers(Numberfilename);
	  Numbers<<"Number of scenarios.."<<endl;
	  Numbers<< '\t';
	  Numbers<< "I"<< '\t'<<"E"<<'\t'<<"U"<<endl;

	  const char* Number_post_filename  = "EX_Numbers_post.txt";
	  ofstream Numbers_post(Number_post_filename);
	  Numbers_post<<"Number of scenarios.."<<endl;
	  Numbers_post<< '\t';
	  Numbers_post<< "I"<< '\t'<<"E"<<endl;
	 
	  const char* Scenariofilename  = "EX_Average.txt";
	  ofstream ScenarioRes(Scenariofilename);

	  const char* Optimafilename  = "EX_Optima.txt";
	  ofstream OptimaRes(Optimafilename);


	  for (int ii=0; ii<No_instance; ii++){
		  IloEnv env;
		  IloNum rho=config[ii];
	
/////////////////// DECISION VARIABLES WITH NAMES  /////////////////////////////
		  
	IloInt MaxCut =10;
	//Upper letters -- Variables
	char varName[100];

	//define X variables
	IloNumVarArray X(env, numStage, 0, 2000);
	for(s=0;s<numStage;s++){
		sprintf_s(varName, "X_%d", (int) s);
		X[s].setName(varName);
	}
	
	//define R variables (stock)
	IloNumVarArray R(env, numStage, 0, 2000);
	for(s=0;s<numStage;s++){
		sprintf_s(varName, "R_%d", (int) s);
		R[s].setName(varName);
	}

	
	
	IloNumVarArray	Theta(env,numScen, -10000, 10000);
	
	for(w=0;w<numScen;w++){
		sprintf_s(varName, "Theta_%d",(int) w);
		Theta[w].setName(varName);
	}	
	

	
	IloNumVar	Alpha(env, -10000, 10000);
	
	sprintf_s(varName, "Alpha");
	Alpha.setName(varName);
	

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

	//Lower case letters -- passed values for the variables 
	IloNumArray x(env, numStage);
	IloNumArray r(env, numStage); 

	IloNum Best_X;

	


	IloNumArray p(env, numScen);
	   for(w=0; w<numScen; w++){
		   p[w]=prob[w];
	   }

    //Dual values 
	IloNumArray pi_inventory(env, numScen);  //dual variable for inventory type constraints
	
	
	

	////////PARAMETER to do Post-processing///////////////

	
	   Scenario *Output = new Scenario[numScen];
	   

	   //Scenario Output[numScen];

	   for(w=0; w<numScen; w++){
		   Output[w].partial_sum=0;
		   Output[w].No=w;
		   Output[w].nom_Prob=prob[w];
		   Output[w].demand=demand[w];
		   Output[w].Status=NULL;
		   Output[w].new_ObjCost=-1;
		   
	   }
	   
	    
////////////BUILD MASTER AND SUBPROBLEM MODELS  //////////////////////////

	cout<<"building master and subproblems"<<endl;
    // arrays for model and lp
    
    IloModel MODEL(env);
	IloCplex CPX(env);
	
		 
	MODEL.add(IloMinimize(env, -cost[2]*X[0]+Alpha)); 
	
	CPX.extract(MODEL); 

	IloModel model_Sub(env);
	IloCplex cplex_Sub(env);

	//IloObjective obj_Sub(env);
	//obj_Sub.setSense(IloObjective::Minimize); 
 	//model_Sub.add(obj_Sub);
	model_Sub.add(IloMinimize(env, X[1] * (cost[0] + cost[2]) + R[1] * (cost[1] - cost[2])));
	IloRange consInventory;
	consInventory=X[1]-R[1]==0;
	model_Sub.add(consInventory);
	cplex_Sub.extract(model_Sub);
		
	
	//model for sub (max) problem
	IloModel  modSub(env);
	IloObjective objSub(env);
	objSub.setSense(IloObjective::Maximize); 
 	modSub.add(objSub);
	IloCplex cplexSub(modSub);
	IloRangeArray conDistancePos(env);
	
	IloRangeArray conDistanceNeg(env);
	for(w=0;w<numScen; w++){ 
		conDistancePos.add(P[w]-Z[w]<=prob[w]);
		conDistanceNeg.add(-P[w]-Z[w]<=-prob[w]);
	}
	modSub.add(conDistancePos);
	modSub.add(conDistanceNeg);
	
	IloExpr conDistance(env);
	IloRange RngDistance;
	IloExpr conProb(env);
	IloRange RngProb;
	for(w=0;w<numScen; w++){ 
		conDistance+=Z[w];
		conProb+=P[w];
	}
	RngDistance=conDistance<=2*rho;
	modSub.add(RngDistance);
	RngProb=conProb==1;
	modSub.add(RngProb);
	
	//////////BEGIN ITERATIONS/////////////////////////////////
		
    cout<<"begin iteration..."<<endl;
	cout.precision(10);
    //Initialize LB, UB
	IloNum rel_Gap = IloInfinity;
	IloNum LB = -IloInfinity;
	IloNum UB = IloInfinity;
    //Scenario cost
	IloNumArray subObj_hat(env, numScen);
	//IloNumArray subObj_est(env, numScenNode);

    IloNum z_hat=0;
     char resName[100];   
    sprintf_s(resName, "EX_%.2f.txt", rho);
	const char* Resfilename  = resName;
	ofstream Result(Resfilename);

	std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	Result.precision(10);
	
	Result<<"Number of scenarios:"<< numScen<<endl;
		Result<<endl;
	
	Result<< "Solve the original problem.."<<endl;
	Result <<"Iter"<<'\t'<<"Z_hat"<<'\t'<<"LB"<<'\t'<<"UB"<<'\t'<<"Relative_Gap"<<'\t'<<"Incumbent"<<endl;
		
			
	nCut=0;
	///////////////////////////////////////////////////////////
	while (rel_Gap > toler) {

	
			if( !CPX.solve () ){
				env.error() << "Master problem infeasible" << endl;
				throw(-1);
			}  
			
  
			x[0]=CPX.getValue(X[0]); 
			LB = CPX.getObjValue();
			
			


			IloExpr MaxCut(env);
			
			//root node objective function
			z_hat = 0;
			
		
		
		//subObj_est[0] = CPX[0].getObjValue();
		
		for(w=0; w<numScen; w++) {
			/*IloExpr OBJ=obj_Sub.getExpr();
			OBJ=X[1]*UU+R[1]*EE-VV*demand[w];
			obj_Sub.setExpr(OBJ);*/
			consInventory.setBounds(x[0]-demand[w], x[0] - demand[w]);
			cplex_Sub.solve();
			x[1]=cplex_Sub.getValue(X[1]);
			r[1]=cplex_Sub.getValue(R[1]);
			subObj_hat[w] = x[1] * (cost[0] + cost[2]) + r[1] * (cost[1] - cost[2]);
			pi_inventory[w] = cplex_Sub.getDual(consInventory);



				IloExpr Gx(env);
						
						
				Gx=pi_inventory[w]*(X[0]- demand[w]);
								
				
				MODEL.add(Theta[w] >= Gx); 
				
		}//end for


		for(w=0;w<numScen;w++){	
			objSub.setLinearCoef(P[w],subObj_hat[w]);
		}
			
		for(w=0;w<numScen;w++){		
			conDistancePos[w].setUb(prob[w]);
			conDistanceNeg[w].setUb(-prob[w]);
		}
		//cplex_row.exportModel("prob.lp");
		cplexSub.solve();

			
		cplexSub.getValues(p, P);
					
			

		MaxCut=IloScalProd(Theta, p);
		MODEL.add( Alpha >= MaxCut);
		

		cout<<"calculating gap"<<endl;
		///*****Calculate UB*****////
		
		
				
		//*****Calculate gap****////
			z_hat=IloScalProd(subObj_hat, p);		
			z_hat+=-cost[2]*x[0];
			if (z_hat<= UB || abs(z_hat-UB)/abs(UB)<= toler){
				//store results
				UB = z_hat;
				Best_X = x[0]; 
				

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
			//grad_X.end();
		} //end while   		    


		Result<<endl;

		Details << rho << '\t';

		IloNum Obj=LB;
		CPX.end();
		MODEL.end();
		cout<<"DONE"<< endl;


		//now obtain multiple optimal solutions
		
		IloNumArray optX(env, 2);
		IloNumVarArray X2 (env, numScen, 0, IloInfinity);
		for (w=0; w<numScen; w++){
			sprintf_s(varName, "Y_%d",(int) w);
			X2[w].setName(varName);
		}
		IloNumVar XX (env, 0, IloInfinity);
		sprintf_s(varName, "X");
		XX.setName(varName);
		IloNumVar lambda (env, 0, IloInfinity);
		sprintf_s(varName, "eta");
		lambda.setName(varName);
		IloNumVarArray V(env, numScen, 0, IloInfinity);
		for (w=0; w<numScen; w++){
			sprintf_s(varName, "V_%d",(int) w);
			V[w].setName(varName);
		}
	
		IloNumVarArray R2(env, numScen, 0, IloInfinity);
		for (w=0; w<numScen; w++){
			sprintf_s(varName, "R_%d",(int) w);
			R2[w].setName(varName);
		}

		IloNumVar	Theta2(env, 0, IloInfinity);
		sprintf_s(varName, "Theta");
		Theta2.setName(varName);
		
		IloModel model_master2(env);
		IloCplex cpx_master2(env);
	
		IloObjective objMaster2(env);
		objMaster2.setSense(IloObjective::Minimize); 
		model_master2.add(objMaster2);
			


		for (w=0; w<numScen; w++){
			model_master2.add(Theta2>=X2[w]*(cost[0]+cost[2])+R2[w]*(cost[1]-cost[2]));
			model_master2.add(X2[w]-R2[w]==XX-demand[w]);
			model_master2.add(V[w]>= X2[w] * (cost[0] + cost[2]) + R2[w] * (cost[1] - cost[2]) -lambda);
		}
	
		cpx_master2.extract(model_master2); 

		IloExpr objExpr;
		objExpr=-cost[2]*XX+rho*Theta2+ (1-rho)*lambda;
		for (w=0; w<numScen; w++){

			objExpr+=prob[w]*V[w];					
		}
		model_master2.add(objExpr==Obj);
		//obtain a in [a,b] range of optimal solutions
		objMaster2.setLinearCoef(XX, 1);
		cpx_master2.solve();
		optX[0]=cpx_master2.getValue(XX);
		//obtain b in [a,b] range of optimal solutions
		objMaster2.setLinearCoef(XX, -1);
		cpx_master2.solve();
		optX[1]=cpx_master2.getValue(XX);	

		model_master2.end();
		cpx_master2.end();

		OptimaRes<< rho<<'\t';
		OptimaRes<< "["<<optX[0]<<", "<<optX[1]<<"]"<<endl;
	
		Details<<Best_X << '\t';

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
		if(Output[0].partial_sum >= rho)
			VaR_index.push_back(0);

		//if (rho!=2){
			if (VaR_index.size()==0){
				for (w=1; w<numScen+1; w++){
					if(Output[w].partial_sum >= rho &&  Output[w-1].partial_sum < rho){
						VaR_index.push_back(w);
						break;
					}
				}
			}
		//}
		//else{
		//	VaR_index.push_back(numScen-1);
		//}
		
		
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
					Details << Output[*it].No << "]" << '\t';
				}
				else {
					it--;
					Details << Output[*it].No << ", ";
				}
			}
		}
		else
			Details << '\t';

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
			Result<<Output[w].demand<<'\t';
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

		//

		/*
		for (it = Effscen[ii].begin(); it != Effscen[ii].end(); it++) {
			Allscen[ii].push_back(*it);
		}
		Effscen[ii].clear();
		for (it = Ineffscen[ii].begin(); it != Ineffscen[ii].end(); it++) {
			Allscen[ii].push_back(*it);
		}
		Ineffscen[ii].clear();
		for (it = Unknownscen[ii].begin(); it != Unknownscen[ii].end(); it++) {
			Allscen[ii].push_back(*it);
		}
		Unknownscen[ii].clear();
		*/
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
				MODEL.add(IloMinimize(env, -cost[2]*X[0]+Alpha)); 
				MODEL.add(X[0]==R[0]);
				
       
				CPX.extract(MODEL); 
			   for(w=0; w<numScen; w++){
				   SubOutput[w].nom_Prob=prob[w];
			   }
				IloRange RngZero;
				RngZero=P[ww]==0;
				modSub.add(RngZero);

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
					x[0]=CPX.getValue(X[0]); 
					LB = CPX.getObjValue();


					IloExpr MaxCut(env);
			
					//root node objective function
					z_hat = 0;

					//Solve subproblem dual for each scenario
					for(w=0; w<numScen; w++) {
						/*IloExpr OBJ=obj_Sub.getExpr();
						OBJ=X[1]*UU+R[1]*EE-VV*demand[w];
						obj_Sub.setExpr(OBJ);*/
						consInventory.setBounds(x[0]-demand[w], x[0]-demand[w]);	
						cplex_Sub.solve();
						x[1]=cplex_Sub.getValue(X[1]);
						r[1]=cplex_Sub.getValue(R[1]);
						subObj_hat[w] =x[1]*(cost[0]+cost[2])+r[1]*(cost[1]-cost[2]);
						pi_inventory[w] = cplex_Sub.getDual(consInventory);


						IloExpr Gx(env);
						
						
						Gx=pi_inventory[w]*(X[0]-demand[w]);
								
				
						MODEL.add(Theta[w] >= Gx); 
				
					}//end for
				

					//Solve row-generion subproblem and calculationg p
			
					for(w=0;w<numScen;w++){	
						objSub.setLinearCoef(P[w],subObj_hat[w]);
					}
			
					for(w=0;w<numScen;w++){		
						conDistancePos[w].setUb(prob[w]);
						conDistanceNeg[w].setUb(-prob[w]);
					}
					//cplex_row.exportModel("prob.lp");
					cplexSub.solve();

					if( cplexSub.solve () ){ 
						cplexSub.getValues(p, P);
					
			

						MaxCut=IloScalProd(Theta, p);
						MODEL.add( Alpha >= MaxCut);

						//*****Calculate gap****////
						z_hat=-cost[2]*x[0];
						z_hat+=IloScalProd(subObj_hat, p);		
						if (z_hat<= UB || abs(z_hat-UB)/abs(UB)<= toler){
							//store results
							UB = z_hat;
							Best_X = x[0]; 
			
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
		
					if (ctnWhile){
						if (abs(Output[*it].new_ObjCost- Obj)>toler){
							Output[*it].Status="Effective";
							category_effect[1]++;
							Effscen[ii].push_back(*it);
						}
						else{
							category_effect[0]++;
							Output[*it].Status="Ineffective";
							Ineffscen[ii].push_back(*it);
						}
					}
					else{
						category_effect[1];
						Output[*it].Status="Effective";
						Effscen[ii].push_back(*it);
					}
					modSub.remove(RngZero);
					delete[] SubOutput;
					MODEL.end();
					CPX.end();
					}
					else{
						category_effect[1];
						Output[*it].Status="Effective";
						Effscen[ii].push_back(*it);
					}
					Result<<endl;
				}//end for ww

				Result<<"Information about (primarily) unknown scenarios.."<<endl;
				vector<double>::iterator itt;
				Result << "w"<<'\t'<<"f_n(x*)"<<'\t'<<"f(bar{x})"<<'\t'<<"Status"<<endl;

				for (it= Unknownscen[ii].begin(); it!= Unknownscen[ii].end(); it++){
					Result<< Output[*it].No<< '\t';
					Result<<Output[*it].demand<<'\t';
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

		Summary<< rho<<'\t';
		for (w=0; w<numScen; w++){
			if (Output[w].Status=="Effective"){
				Effscen[ii].push_back(w);
				Summary<< "E"<< '\t';	
			}
			else if (Output[w].Status=="Ineffective"){
				Ineffscen[ii].push_back(w);
				Summary<< "I"<< '\t';	
			}
		}
		Summary<<endl;

		Details << "[";
		for (it = Effscen[ii].begin(); it != Effscen[ii].end(); it++) {
			if (++it == Effscen[ii].end()) {
				it--;
				Details << *it << "]" << endl;
			}
			else {
				it--;
				Details << *it << ", ";
			}
		}


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




				
		Result.close();
		MODEL.end();
		modSub.end();
		cplexSub.end();
		env.end();
}//end instances loop


		for (w=0; w<numScen; w++){
			ScenarioRes<< w<< '\t'<< demand[w]<<'\t';	
			ScenarioRes<<Count_Effective[w]<<'\t'<< Count_Known[w]<<endl;


		}

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
						Nested << "Found " << *it << " with probability " << prob[*it] << " in gamma " << (double)ii/10 << " not exits in gamma " << (double)(ii - 1)/10<< endl;
						if (nestedness == 1)
							nestedness = 0;
						//break;
					}
				}
				if (nestedness == 1) {
					Nested << "gamma " << (double)ii/10 << " is nested in gamma " << (double)(ii - 1)/10 << endl;
				}
				
				if (Effscen[ii].size()>Effscen[ii - 1].size()) {
					Nested << "********** gamma " << (double)ii/10 << " is larger than " << (double)(ii - 1)/10 << endl;
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