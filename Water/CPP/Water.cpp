/*
******************************************
* Author:   Hamed Rahimian
*          The Ohio State University
*          January 01, 2015
******************************************
*/

#include <ilcplex/ilocplex.h>
#include <stdlib.h>   
#include <vector>
#include <time.h>
#include <algorithm> 
#include <cmath>
using namespace std;

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> VarArray2;
typedef IloArray<VarArray2>      VarArray3;
typedef IloArray<IloRangeArray>  EquationArray;
typedef IloArray<EquationArray>  EquationArray2;
typedef IloArray<IloNumArray2>   IloNumArray3;
typedef IloArray<IloNumArray3>   IloNumArray4;
typedef IloArray<IloExpr>        ObjArray;  
typedef IloArray<IloCplex>       CplexArray;            
typedef IloArray<IloModel>       ModelArray;

struct Scenario{
	IloInt No;
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
	
	IloNum toler = 0.000001, disRate = 0.02;
	IloInt i, j, k, t, w, w2, s,ww;
	const IloInt No_type = 10, No_node = 62;
	IloInt No_stage=4;
	
	// The data below read from input file
	IloNum discRate = 0.04; //discount Rate
	IloIntArray type(envData), No_year(envData);                      //The type of corresponding node
	IloIntArray startNode(envData), endNode(envData);
	IloInt No_link; 	//Number of year in each stage 
	IloNumArray  cost(envData);   //cost on each link
	IloNum   returnRate;
	IloNumArray   loss(envData), capacity(envData), Storage0(envData),  storageUB(envData);
	///////////////// DATA FILE READING ////////////////////////////////
		   
	const char* filename1  = "Data\\input1_800.dat";	      
		
	ifstream file1(filename1);
	if (!file1) {
		cerr << "ERROR: can not open file '" << filename1 << "' for reading" << endl;
		cerr << "usage:   " << argv[0] << " <file>" << endl;
		throw(-1);
	}
		   
    file1 >> No_year;
	file1 >> type >> startNode >> endNode >> cost ;
	file1 >> capacity >>  returnRate; 
	file1 >>  loss >> Storage0>>storageUB;

	No_link = startNode.getSize();
	IloInt No_yearMax = 15 ; 
	 
 	IloBool consistentData = ( No_year.getSize() == No_stage && endNode.getSize() == No_link 
		&& cost.getSize() == No_link && type.getSize() == No_node );
           
	if (!consistentData) { 
 		cerr << "ERROR: Inconsistent data1!" << endl;
		throw(-1);
	}
	   
	file1.close();
	cout<<"Data1 reading down"<<endl;
  ///////////////////Identify type group////////////////////////
	//identify role of each node in the water network, and add node number to their corresponding role vector
	//it acts to know which constraints are needed for each node

	vector<int> userID;
	vector<int> potUserID;
	vector<int> rechargeID;
	vector<int> balanceID;
	vector<int> capacityID;

	for(i = 0; i<No_node; i++){
		if(type[i] == 4){
			rechargeID.push_back(i);
		}
		else if(type[i] == 0 || type[i] == 2 || type[i] == 5 || type[i] == 6 || type[i] == 7){
			balanceID.push_back(i); 
		}
		else if(type[i] == 1){
			capacityID.push_back(i); 
		}
		else if(type[i] == 8 ){
			userID.push_back(i);
			potUserID.push_back(i);
		}
	}
    for(i = 0; i<No_node; i++){
		if(type[i] == 9 ){
			userID.push_back(i);
		}
		else if(type[i] == 4){
			capacityID.push_back(i); 
		}
	}
	const IloInt No_user = int(userID.size());
	const IloInt No_potUser = int(potUserID.size());
	const IloInt No_balance = int(balanceID.size());
    const IloInt No_capacity = int(capacityID.size());
    const IloInt No_recharge = int(rechargeID.size());
	cout<<No_user<<"  "<<No_potUser<<endl;
		 
///////////////////read demeand, capacity, and intial storage//////////////////
          
	IloNumArray3   demand(envData,No_stage);
	IloNumArray2   population(envData,No_stage), TucsonPop(envData,No_stage);
	IloNumArray    demUnit(envData,No_user);
	for(s = 0; s< No_stage; s++){
		demand[s] = IloNumArray2(envData,No_yearMax);
		population[s] = IloNumArray(envData,No_yearMax);
		TucsonPop[s] = IloNumArray(envData,No_yearMax);
		for(t=0; t<No_yearMax; t++)
			demand[s][t] = IloNumArray(envData,No_user);
	}
	const char* filename2  =  "Data\\input2.dat";
	ifstream file2(filename2);
	if (!file2) {
		cerr << "ERROR: can not open file '" << filename2 << "' for reading" << endl;
		cerr << "usage:   " << argv[0] << " <file>" << endl;
		throw(-1);
	}
           
	file2 >> population[0]>>population[1]>>population[2]>>population[3]; 
	file2 >> demUnit; 
	file2 >> TucsonPop[0]>> TucsonPop[1]>> TucsonPop[2]>> TucsonPop[3];
		  
	file2.close();
	cout<<"Data part 2 read done"<<endl;
		 
	const char* filename3  =  "Data\\DemPortion.txt";
	ifstream file3(filename3);
	if (!file3) {
		cerr << "ERROR: can not open file '" << filename3 << "' for reading" << endl;
		cerr << "usage:   " << argv[0] << " <file>" << endl;
		throw(-1);
	}
	IloNumArray2 DembyNode(envData,No_year[0]+No_year[1]+No_year[2]+No_year[3]);
	for(i= 0; i<No_year[0]+No_year[1]+No_year[2]+No_year[3];i++){
		DembyNode[i]=IloNumArray(envData,No_potUser);
		/*file3>>DembyNode[i];*/
	}
	for(i= 0; i<No_year[0]+No_year[1]+No_year[2]+No_year[3];i++){ 
		file3>>DembyNode[i];
	}
       
	file3.close(); 
	cout<<"Data part 3 read done"<<endl; 
          
	for(t=0; t<No_year[0]; t++){
		for(i=0; i<No_potUser;i++){
			demand[0][t][i] = double(135*0.00112*population[0][t]*DembyNode[t][i]*0.8);
		}
	} 
	for(t=0; t<No_year[1]; t++){
		for(i=0; i<No_potUser;i++){
			demand[1][t][i] = double(135*0.00112*population[1][t]*DembyNode[t + No_year[0]][i]*0.8);
		}
	}
	for(t=0; t<No_year[2]; t++){
		for(i=0; i<No_potUser;i++){
			demand[2][t][i] = double(135*0.00112*population[2][t]*DembyNode[t+ No_year[0]+ No_year[1]][i]*0.8);
		}
	}
	for(t=0; t<No_year[3]; t++){
		for(i=0; i<No_potUser;i++){
			demand[3][t][i] = double(135*0.00112*population[3][t]*DembyNode[t+ No_year[0]+ No_year[1]+ No_year[2]][i]*0.8);
		}
	}
	for(s=0;s<No_stage;s++){
		for(t=0; t<No_yearMax; t++){
			for(i=No_potUser; i<No_user;i++){
				demand[s][t][i] = demand[s][t][i-No_potUser]*0.25;
			}
		}
	}
		   
	IloNumArray2  CAPamt(envData,No_stage);
	for(s = 0; s< No_stage; s++){
		CAPamt[s] = IloNumArray(envData,No_yearMax);
		for(t=0; t<No_yearMax; t++){
			CAPamt[s][t] = 144000*population[s][t]/TucsonPop[s][t] ;
		}
	}
    population.end();TucsonPop.end(); demUnit.end();DembyNode.end();
		  
        
	///////////////Define scenario tree nodes 
	//revalue No_stage
	No_stage=2;
	IloInt No_scen =200; //Number of scenarios in each stage

	IloInt No_scenNode= (pow(No_scen,No_stage)-1)/(No_scen-1);//Total number of nodes in the scenario tree 
	IloIntArray No_scen_stage(envData, No_stage);	  
	for(s=0; s<No_stage;s++){ 
		No_scen_stage[s]=pow(No_scen,s);
	}
	IloIntArray No_scen_stage_partial(envData, No_stage);	  
	No_scen_stage_partial[0]=1;
	for(s=1; s<No_stage;s++){ 
		No_scen_stage_partial[s]=No_scen_stage_partial[s-1]+pow(No_scen,s);
	}
	
	IloInt No_scenElm = No_yearMax+2;
	IloNumArray prob(envData, No_scenNode);
	for(w=1;w<No_scenNode; w++){
		prob[w]= (double)1/No_scen;
	}

		  
	IloNumArray2 scenFac(envData);
	const char* filename4  =  "Data\\scenFac_200.txt";
	ifstream file4(filename4);
	file4>>scenFac;
	file4.close();

	//determine stage of each node in the tree
	IloIntArray stageMap(envData, No_scenNode);
	stageMap[0] = 0;
	for(w=1; w<No_scenNode; w++){
		for(s=1; s<No_stage;s++){ 
			if( (w>=No_scen_stage_partial[s-1]) && (w<No_scen_stage_partial[s]) ){
				stageMap[w] = s;
				break;
			}
		}
	}
	
    
	IloIntArray2 ancestor(envData, No_scenNode);
	for(w=0; w<No_scenNode; w++){
		ancestor[w] = IloIntArray(envData, No_scenNode);
		for(w2=0; w2<No_scenNode; w2++){
			ancestor[w][w2] = 0;
		}
	}

	for(w=0; w<No_scen_stage_partial[No_stage-2]; w++){
		for(j= 1+No_scen*w; j<1+No_scen*(w+1);j++){
			//w is ancestor of j
			ancestor[w][j] = 1;
		}
	}
	
	    
	//denotes decendant of each node in the scenario tree by an array, whose elements are the child node number
	IloIntArray2 decendant(envData, No_scenNode);
	for(w=0;w<No_scenNode;w++){
		decendant[w]= IloIntArray(envData,No_scen);
		for(j=0; j<No_scen;j++){
			decendant[w][j] = 0;
		}
	}
	
	for(w=0; w<No_scen_stage_partial[No_stage-2]; w++){
		for(j=0;j<No_scen;j++){
			decendant[w][j] = 1+ No_scen*w+j;
		}
	}
	

	//vector to stores each gamma ineffective scenarios
	  vector<int> * Ineffscen= new vector<int> [21];

	  //vector to stores each gamma effective scenarios
	  vector<int> * Effscen= new vector<int> [21];
	  //vector to stores each gamma new cost after removal of an effective scenarios
	  vector<double> * NewCost_Effscen= new vector<double> [21];

	  //vector to stores each gamma undetermined scenarios
	  vector<int> * Unknownscen= new vector<int> [21];

	  
	  IloNumArray Count_Effective(envData, No_scenNode);
	  for (w=0; w<No_scenNode; w++){
		  Count_Effective[w]=0;
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

	  const char* Eff_file_name  = "Eff_Change.txt";
	  ofstream EffScenarioRes(Eff_file_name);
	  for (i=0; i<21; i++){
		EffScenarioRes<< 0.1*i<< '\t';
	  }
	  EffScenarioRes<<endl;

	

///******
	for (int ii=0; ii<No_instance; ii++){
		  IloEnv env;
		  IloNum rho=config[ii];
	    
 /////////////////// DECISION VARIABLES WITH NAMES  /////////////////////////////
		  
		IloInt nCut = 0;
		IloInt MaxCut =20;

		//Upper letters -- Variables
		
		char varName[100];

		//define x variables
		VarArray3 Q(env, No_link);
		for(i =0; i<No_link; i++){
			Q[i] = VarArray2(env, No_stage);
			for(s = 0; s<No_stage; s++){
				Q[i][s] = IloNumVarArray(env, No_yearMax, 0, IloInfinity, ILOFLOAT); 
				for(t=0;t<No_yearMax;t++){
					sprintf_s(varName, "Q_%d_%d_%d",(int) i, (int) s, (int) t);
					Q[i][s][t].setName(varName);
				}
			}
		}
		//define y variables
		VarArray3 Storage(env, No_recharge);
		for(i =0; i<No_recharge; i++){
			Storage[i] = VarArray2(env, No_stage);
			for(s = 0; s<No_stage; s++){
				Storage[i][s]= IloNumVarArray(env, No_yearMax, 0, IloInfinity, ILOFLOAT);
				for(t=0;t<No_yearMax;t++){
					sprintf_s(varName, "S_%d_%d_%d",(int) i, (int) s, (int)t);
					Storage[i][s][t].setName(varName);
				}
			}
		}

		VarArray2	Theta(env,No_stage);
		for(s=0; s<No_stage; s++){
			Theta[s] = IloNumVarArray(env, No_scen, 0, IloInfinity);
			for(j=0;j<No_scen;j++){
				sprintf_s(varName, "Theta_%d_%d",(int) s, (int) j);
				Theta[s][j].setName(varName);
			}	
		}

	
		IloNumVarArray	Alpha(env, No_stage, 0, IloInfinity);
		for(s=0;s<No_stage;s++){
			sprintf_s(varName, "Alpha_%d", (int) s);
			Alpha[s].setName(varName);
		}

		IloNumArray alpha(env, No_scenNode);
		for(w=0; w<No_scenNode; w++){
			alpha[w]= 0;
		}
		


		//Lower case letters -- passed values for the variables 
		IloNumArray3 storage_hat(env, No_recharge);
		for(i =0; i<No_recharge; i++){
			storage_hat[i] = IloNumArray2(env, No_scenNode); 
			for(w=0;w<No_scenNode;w++){
				storage_hat[i][w] = IloNumArray(env, MaxCut);
				for(nCut=0;nCut<MaxCut;nCut++){
					storage_hat[i][w][nCut] = 0;
				}
			}
		} 


        //Dual values 
		IloNumArray3 pi_balance(env, No_yearMax);  //dual variable for nodal flow balance constraints
		IloNumArray3 pi_capacity(env, No_yearMax);  //dual variable for pump/surface water capacity constraints
		IloNumArray3 pi_demand(env, No_yearMax);  //dual variable for demand satisfication constraints
		IloNumArray3 pi_return(env, No_yearMax);  //dual variable for potable return constraints
		IloNumArray3 pi_safe(env, No_yearMax);  //dual variable for safe yield
		IloNumArray2 pi_storage(env, No_scenNode);  //dual variable for storage balance constraints
		IloNumArray2 pi_flowBound(env, No_scenNode);  //dual variable for RF outflow bound constraints 

		for(t=0;t<No_yearMax;t++){
			pi_balance[t] = IloNumArray2(env, No_scenNode);
			pi_capacity[t] = IloNumArray2(env, No_scenNode);
			pi_demand[t] = IloNumArray2(env, No_scenNode);
			pi_return[t] = IloNumArray2(env, No_scenNode);
			pi_safe[t] = IloNumArray2(env, No_scenNode);
			for(w=0;w<No_scenNode;w++){
				pi_balance[t][w] = IloNumArray(env, No_balance);
				pi_capacity[t][w] = IloNumArray(env, No_capacity);
				pi_demand[t][w] = IloNumArray(env, No_user);
				pi_return[t][w] = IloNumArray(env, No_potUser);
				pi_safe[t][w] = IloNumArray(env, No_recharge);
				pi_storage[w] = IloNumArray(env, No_recharge);
				pi_flowBound[w] = IloNumArray(env, No_recharge); 

				for(i=0; i<No_balance; i++) 
 	   				pi_balance[t][w][i] = 0; 
				
				for(i=0; i<No_capacity; i++) 
 	   				pi_capacity[t][w][i] = 0;
				
				for(i=0; i<No_user; i++) 
 	   				pi_demand[t][w][i] = 0;
				
				for(i=0; i<No_potUser; i++) 
 	   				pi_return[t][w][i] = 0;

				for(i=0; i<No_recharge; i++) 
 	   				pi_safe[t][w][i] = 0;
		   		
				for(i=0; i<No_recharge; i++) 
 	   				pi_storage[w][i] = 0;
		        
				for(i=0; i<No_recharge; i++) 
 	   				pi_flowBound[w][i] = 0;
			}
		}

      //****
		IloNumArray3  pi_optCut(env, No_scenNode);
        for(w=0;w<No_scenNode;w++){
			pi_optCut[w] = IloNumArray2(env,No_scen);
			for(j=0;j<No_scen;j++){
				pi_optCut[w][j] = IloNumArray(env,MaxCut);
				for(nCut=0;nCut<MaxCut;nCut++){
					pi_optCut[w][j][nCut] = 0;
				}
			}
		}
		 
		////////DECISION VARIABLES AND PARAMETERS FOR Row-generation Subproblem///////////////

		IloNumVarArray P(env, No_scen, 0, IloInfinity);
	   for(w=0;w<No_scen;w++){
		   sprintf_s(varName, "P_%d", (int) w);
		   P[w].setName(varName);
	   }

	   IloNumVarArray Z(env, No_scen, 0, IloInfinity);
	   for(w=0;w<No_scen;w++){
		   sprintf_s(varName, "Z_%d", (int) w);
		   Z[w].setName(varName);
	   }

	   IloNumArray2 p(env, No_scenNode);
		for(w=0; w<No_scenNode; w++){
			p[w] = IloNumArray(env, MaxCut); 
			for(nCut=0;nCut<MaxCut;nCut++){
				p[w][nCut] = 0;
			}
		}
	
	   for(nCut=0;nCut<MaxCut;nCut++){
		   p[0][nCut] = 1;
	   }

		

		for(w=1; w<No_scenNode; w++){
			p[w][0]=prob[w];
		}
		
        
		////////PARAMETER to do Post-processing///////////////

	
		Scenario *Output = new Scenario[No_scenNode];
	   

	   //Scenario Output[numScen];

		for(w=0; w<No_scenNode; w++){
		   Output[w].No=w;
		   Output[w].Status=NULL;   
		   Output[w].new_ObjCost=-1;
		   Output[w].partial_sum=0;
		   Output[w].nom_Prob=prob[w];
		   
	   }

////////////BUILD MASTER AND SUBPROBLEM MODELS  //////////////////////////

		cout<<"building master and subproblems"<<endl;
        // arrays for model and lp
        CplexArray CPX(env, No_stage);
        ModelArray MODEL(env, No_stage);

		//Objective function array
		ObjArray      OBJ(env, No_stage);

		//Constraints array
		 EquationArray2 FlowBalance(env, No_stage);
		 EquationArray2 CapBound(env,  No_stage);
         EquationArray2 MeetDemand(env, No_stage);
		 EquationArray2 ReturnFlow(env,  No_stage); 
		 EquationArray2 SafeYield(env, No_stage); 
         EquationArray2 StorageBalance(env, No_stage); 
		 EquationArray2 RFoutflowBound(env,  No_stage);
		 

		 //Optimality Cut
		 //****
		 EquationArray2   optCut(env, No_stage);
		 EquationArray  maxCut(env, No_stage);
		 //****
		 

		 for(s = 0; s<No_stage; s++){
			 MODEL[s] = IloModel(env);
			 CPX[s]   = IloCplex(env);
			
			 OBJ[s]   = IloExpr(env);
			 FlowBalance[s] = EquationArray(env, No_yearMax);
			 CapBound[s]    = EquationArray(env, No_yearMax);
			 MeetDemand[s]  = EquationArray(env, No_yearMax);
			 ReturnFlow[s]  = EquationArray(env, No_yearMax);
			 SafeYield[s]  =  EquationArray(env, No_yearMax);
			 StorageBalance[s] = EquationArray(env, No_yearMax);
			 RFoutflowBound[s] = EquationArray(env, No_yearMax); 
 		     //Obj definition at each stage s 
			 for(i =0; i<No_link; i++){
				for(t=0;t<No_year[s];t++){
					if(s==1){
					OBJ[s] += pow(1+ disRate,-double(No_year[0]+ t))*cost[i]*Q[i][s][t] ;	
					}
					else if(s==2){
					OBJ[s] += pow(1+ disRate,-double(No_year[0]+ No_year[1]+ t))*cost[i]*Q[i][s][t] ;	
					}
					else if(s==3){
					OBJ[s] += pow(1+ disRate,-double(No_year[0]+ No_year[1] + No_year[2]+t))*cost[i]*Q[i][s][t] ;	
					}
					else{
					OBJ[s] += cost[i]*Q[i][s][t] ;	
					}
				}
			 }
			 
			 if( s<No_stage-1){ 
				 OBJ[s]+= Alpha[s];	
			 }

			 MODEL[s].add(IloMinimize(env, OBJ[s])); 
			 CPX[s].extract(MODEL[s]); 
			 OBJ[s].end();

			 char name1[100];
			 for(t = 0;t<No_yearMax;t++){
				 FlowBalance[s][t]     = IloRangeArray(env);
			     CapBound[s][t]        = IloRangeArray(env);
				 MeetDemand[s][t]      = IloRangeArray(env); 
				 ReturnFlow[s][t]      = IloRangeArray(env);
				 SafeYield[s][t]       = IloRangeArray(env);
				 StorageBalance[s][t]  = IloRangeArray(env); 
				 RFoutflowBound[s][t] = IloRangeArray(env); 
                 //define balnace
				 for(i=0;i<No_balance;i++){
					sprintf_s(name1, "flowbalance_%d_%d_%d",(int) t,(int) s, (int)i);
					IloExpr LHS(env), RHS(env);
					for(j=0;j<No_link;j++){
						if(endNode[j] == balanceID[i])
							LHS += Q[j][s][t]*(1-loss[j]);      //inflow
					} 
					for(j=0;j<No_link;j++){
						if(startNode[j] == balanceID[i])       
							RHS += Q[j][s][t] ;                 //outflow
					}
					if(t<No_year[s]){
						FlowBalance[s][t].add(LHS - RHS == 0); 
					    FlowBalance[s][t][i].setName(name1);
					}
					LHS.end(); RHS.end();
				 }
				

                 // CAP allocation + recharge facility capacity
				 //which constraints? Surface Water Supply Allotment and Treatment Plant Capacity Bounds?
                 for(i=0;i<No_capacity; i++){  
					 sprintf_s(name1, "CapBound_%d_%d_%d",(int) t,(int) s, (int)i);
					 IloExpr LHS(env);
					 IloNum RHS ; 
					 for(j=0;j<No_link;j++){
						if(startNode[j] == capacityID[i])
							LHS += Q[j][s][t];      //outflow
				     } 
					 RHS = capacity[i];
					 if(t<No_year[s]){
					    CapBound[s][t].add(LHS <= RHS);
					    CapBound[s][t][i].setName(name1);
					 }

					 LHS.end();
				 }
			

                 for(i=0;i<No_recharge; i++){ // recharge facility storage safe yield
					sprintf_s(name1, "SafeYield_%d_%d_%d",(int) i, (int) s, (int)t);
					if(t<No_year[s]){
						SafeYield[s][t].add(Storage[i][s][t] <= storageUB[i]);
						///>=???
						SafeYield[s][t][i].setName(name1);
					}
				 }
                
				 //potable and non-potable users 
				 for(i=0;i<No_user; i++){  
					sprintf_s(name1, "MeetDemand_%d_%d_%d",(int) s, (int) t, (int)i);
					IloExpr LHS(env);
					IloNum RHS; 
					RHS = demand[s][t][i];// right hand side is the demand; inflow = demand
					for(j=0;j<No_link;j++){
						if(endNode[j] == userID[i])
							LHS += Q[j][s][t]*(1-loss[j]); 
					}
					if(t<No_year[s]){
						MeetDemand[s][t].add(LHS == RHS);
					    MeetDemand[s][t][i].setName(name1);
					}
					LHS.end(); 
				 }
                 //  Return flow from potable users
                 for(i=0;i<No_potUser; i++){  
					sprintf_s(name1, "ReturnFlow_%d_%d_%d",(int) s,(int) t, (int)i);
					IloExpr LHS(env);
					IloNum RHS;	 
					RHS = returnRate*demand[s][t][i];                // 98% of the potable uses return to the facility
					for(j=0;j<No_link;j++){
						//why sum???
						if(startNode[j] == potUserID[i])
							LHS += Q[j][s][t];
					}  
					if(t<No_year[s]){
						ReturnFlow[s][t].add(LHS == RHS);
					    ReturnFlow[s][t][i].setName(name1);
					}
					LHS.end();
				 }
				
				 //Storage balance
                 for(i=0;i<No_recharge; i++){                  
					IloExpr LHS(env), RHS(env);
					sprintf_s(name1, "StorageBalance_%d_%d_%d",(int) s,(int) t, (int)i);
					for(j=0;j<No_link;j++){
						if(endNode[j] == rechargeID[i])
							LHS += Q[j][s][t]*(1-loss[j]);      //inflow
					} 
					for(j=0;j<No_link;j++){
						if(startNode[j] == rechargeID[i])       
							RHS += Q[j][s][t] ;                 //outflow
					} 
					if( t == 0){
						//inflow -Storage[i][t] - outflow  =  -Storage0[i]
						LHS += - Storage[i][s][t] ;
						if(s==0){ 
							StorageBalance[s][t].add(LHS - RHS == - Storage0[i]);
						}
						else{
							//I believe it has been taken care of later since it needs ancestor value (storage)??? later in the code it has been taken care of
							StorageBalance[s][t].add(LHS - RHS == 0);
						}
						StorageBalance[s][t][i].setName(name1);
					}
					else{
						//why it is not taking care of ancestor???!! (4.20)? I beleive the report must change and here is OK! like 4.24 later
						//-Storage[i][t] + Storage[i][t-1] + inflow -outflow=  0
						LHS += - Storage[i][s][t] + Storage[i][s][t-1];
						if(t<No_year[s]){
							StorageBalance[s][t].add(LHS - RHS == 0);
							StorageBalance[s][t][i].setName(name1);}
					} 
					LHS.end(); RHS.end();
				 }
				
                 //Recharge facility outflow bound
                 for(i=0;i<No_recharge; i++){ 
					IloExpr LHS(env);
					for(j=0;j<No_link;j++){
						if(startNode[j] == rechargeID[i])
							LHS += Q[j][s][t];
					}
					if(s>0){
						if(t==0)
							RFoutflowBound[s][t].add(LHS  <= 0);
						else{
							if(t<No_year[s])
								//why is it not having ancestor? later in the code it is clear that there is no need to ancestor, means that
								//equation (4.24) in the report must change
								RFoutflowBound[s][t].add(LHS - Storage[i][s][t-1] <= 0);
						}
					}
					else{
						if(t==0)
							RFoutflowBound[s][t].add(LHS - Storage0[i] <= 0);
						else{
							if(t<No_year[s])
								RFoutflowBound[s][t].add(LHS - Storage[i][s][t-1] <= 0);
						}
					}
					LHS.end(); 
				 }
				 if(t<No_year[s]){
					 MODEL[s].add(FlowBalance[s][t]);
					 MODEL[s].add(CapBound[s][t]);
					 MODEL[s].add(MeetDemand[s][t]);
					 MODEL[s].add(ReturnFlow[s][t]);
					 MODEL[s].add(SafeYield[s][t]);
					 MODEL[s].add(StorageBalance[s][t]);
					 MODEL[s].add(RFoutflowBound[s][t]);  
				 }
				}//end year 
				 
//****
             optCut[s] = EquationArray(env, No_scen);
			 maxCut[s] = IloRangeArray(env);
//****
			 
			 if(s<No_stage-1){
				 for(w = 0; w<No_scen;w++){
					 optCut[s][w] = IloRangeArray(env);
				
					 for(nCut = 0; nCut<MaxCut; nCut++){
						 IloExpr OptCutExpr(env) ;
						 for(i=0;i<No_recharge; i++){ 
							 OptCutExpr += Storage[i][s][No_year[s]-1]; 
						 }
						sprintf_s(name1, "optCut_%d_%d_%d",(int) s, (int)w,(int) nCut);
					 					
						 optCut[s][w].add(Theta[s+1][w] - OptCutExpr >= -IloInfinity );
						 optCut[s][w][nCut].setName(name1);
					 }
				 }
				for(nCut = 0; nCut<MaxCut; nCut++){
					IloExpr maxCutExpr(env) ;
					for(w=0;w<No_scen; w++){ 
						maxCutExpr += prob[w+1]*Theta[s+1][w]; 
					}
					
					 sprintf_s(name1, "maxCut_%d_%d",(int) s,(int) nCut);
					 maxCut[s].add(Alpha[s] -maxCutExpr >= -IloInfinity);
                     maxCut[s][nCut].setName(name1);
				}
			 }
			 CPX[s].extract(MODEL[s]);  

		 }//end stage
		 for(s=0;s<No_stage;s++){
			 CPX[s].extract(MODEL[s]); 
		 }
		
		  //////ROW_GENERATION SUBPROBLEM//////////////////////

		 IloModel  model_row(env);
		 IloObjective Objective_row(env);
		 Objective_row.setSense(IloObjective::Maximize); 
		 model_row.add(Objective_row);
		 
		 IloRangeArray conDistancePos(env);
		 IloRangeArray conDistanceNeg(env);
		 
		 for(w=0;w<No_scen; w++){ 
			 conDistancePos.add(P[w]-Z[w]<=prob[w+1]);
			 conDistanceNeg.add(-P[w]-Z[w]<=-prob[w+1]);
		 }

		 model_row.add(conDistancePos);
		 model_row.add(conDistanceNeg);
	
		 IloExpr conDistance(env);
		 IloRange RngDistance;
		 IloExpr conProb(env);
		 IloRange RngProb;
		 
		 for(w=0;w<No_scen; w++){ 
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
		cout.precision(10);
        //Initialize LB, UB
		IloNum rel_Gap = IloInfinity;
		IloNum LB = 0;
		IloNum UB = IloInfinity;
        //Scenario cost
		IloNumArray subObj_hat(env, No_scenNode);
		
        IloNum z_hat=0;
        char resName[100];
		sprintf_s(resName, "%.2f.txt", rho);
		const char* Resfilename  = resName;
		ofstream Result(Resfilename);

		std::cout.setf( std::ios::fixed, std:: ios::floatfield );
		Result.precision(10);

		Result<<"Number of scenarios:"<< No_scen<<endl;
		Result<<endl;
		
		Result<< "Solve the original problem.."<<endl;
		Result <<"Iter"<<'\t'<<"Z_hat"<<'\t'<<"LB"<<'\t'<<"UB"<<'\t'<<"Relative_Gap"<<endl;

       

		IloNumArray4 G(env,No_scenNode);
		IloNumArray3 g(env,No_scenNode);
		for(w2=0; w2<No_scenNode; w2++){
			G[w2] = IloNumArray3(env,No_scen);
			g[w2] = IloNumArray2(env,No_scen);
			for(j = 0;j<No_scen;j++){
				G[w2][j] = IloNumArray2(env,MaxCut);
				g[w2][j] = IloNumArray(env,MaxCut);
				for(nCut=0;nCut<MaxCut;nCut++){
					G[w2][j][nCut] = IloNumArray(env,No_recharge);
					g[w2][j][nCut] = 0;
					for(i=0;i<No_recharge;i++){
						G[w2][j][nCut][i]=0;
					}
				}
			}
		}

		
        
		///////////////////////////////////////////////////////////
        const clock_t begin_time = clock();
        for(nCut = 0; nCut<MaxCut; nCut++){ //main loop
			cout<<"nCut ="<<nCut<<endl;
			
			for(w=0;w<No_scenNode;w++){
				subObj_hat[w] = 0;
			}

            //***Forward pass***///
            //Step 1: solve master problem -- First stage relaxed problem 
			 

            CPX[0].solve();  //solve first stage relaxed problem
			
		    if( !CPX[0].solve () ){
			   env.error() << "Master problem infeaasible" << endl;
			   throw(-1);
		    } 

			//record storage_hat[i][0]
			//to pass to next stage
			for(i=0;i<No_recharge;i++){
			    storage_hat[i][0][nCut] = CPX[0].getValue(Storage[i][0][No_year[0] - 1]);
		    }
			alpha[0]=CPX[0].getValue(Alpha[0]);
			//record obj_master = CPX[0].getObjValue();
			LB = CPX[0].getObjValue();
						 
			subObj_hat[0] = LB-alpha[0];
			
			
			
		
		
			//Step 2: For each node at each stage, update the sceanrio demands
            for(s=1; s<No_stage;s++){ 
				//update subproblem
	            for(w=1; w<No_scenNode; w++){
					if(stageMap[w]==s){ 
						for(t=0;t<No_year[s];t++){
							//randomness for RHS of supply, for capacityID=0
							//scenFac[w][0] is reserved for \Xi_1
						   CapBound[s][t][0].setBounds(-IloInfinity, CAPamt[s][t]*scenFac[w][0]);
		                   for(i=0;i<No_user;i++)
							   //scenFac[w][t+1] corresponds to \Xi_2 (t)
							   MeetDemand[s][t][i].setBounds(scenFac[w][t+1]*demand[s][t][i], scenFac[w][t+1]*demand[s][t][i]);
						   for(i=0;i<No_potUser; i++)
							   ReturnFlow[s][t][i].setBounds(returnRate*scenFac[w][t+1]*demand[s][t][i],returnRate*scenFac[w][t+1]*demand[s][t][i]);
						}
						
						 
						for(w2=0; w2<No_scenNode; w2++){
							if(ancestor[w2][w]==1){
								//w2 is ancestor of w
							   for(i=0;i<No_recharge; i++){
								   StorageBalance[s][0][i].setBounds(-storage_hat[i][w2][nCut],-storage_hat[i][w2][nCut]);
								   RFoutflowBound[s][0][i].setBounds(-IloInfinity, storage_hat[i][w2][nCut]);
								}
							   //note if ancestor of w is found, there is no need to continue loop, since every node has only one ancestor
							}
						}
						////update cuts
						if(s<No_stage-1 && nCut>0){ 
							for(k=0;k<nCut;k++){ 
								maxCut[s][k].setBounds(0, IloInfinity);
								for(j=0;j<No_scen;j++){
									maxCut[s][k].setLinearCoef(Theta[s+1][j], -p[decendant[w][j]][k]);	 
								}
								
								for(j=0;j<No_scen;j++){
									optCut[s][j][k].setBounds(g[w][j][k], IloInfinity);
									for(i=0;i<No_recharge; i++){ 
										optCut[s][j][k].setLinearCoef(Storage[i][s][No_year[s]-1], -G[w][j][k][i]);
									}
								}
							}//end for nCut
						}//end if
						
						
			   // 
						CPX[s].solve();
		 
						for(i=0;i<No_recharge;i++){
							storage_hat[i][w][nCut] = CPX[s].getValue(Storage[i][s][No_year[s] - 1]);
						} 
						//store dual variables at last stage	
						if(s==No_stage-1){
							//they are the same since there is no Theta
							subObj_hat[w] = CPX[s].getObjValue();
							
							
							for(t=0;t<No_year[s];t++){ 
								for(i=0; i<No_capacity; i++) {
 	   								pi_capacity[t][w][i] = CPX[s].getDual(CapBound[s][t][i]); 
								}
								
								for(i=0; i<No_user; i++) {
 	   								pi_demand[t][w][i] = CPX[s].getDual(MeetDemand[s][t][i]); 
								}
								
								for(i=0; i<No_potUser; i++) {
 	   								pi_return[t][w][i] = CPX[s].getDual(ReturnFlow[s][t][i]); 
								}

								for(i=0; i<No_recharge; i++){ 
 	   								pi_safe[t][w][i] = CPX[s].getDual(SafeYield[s][t][i]);
								}
								
							}
							for(i=0; i<No_recharge; i++){ 
								//we need it to obtain coefficent of the Storage[i][s][No_year[s] - 1] in the cut
 	   							pi_storage[w][i] = CPX[s].getDual(StorageBalance[s][0][i]); 
							}
						        
							for(i=0; i<No_recharge; i++){ 
								//we need it to obtain coefficent of the Storage[i][s][No_year[s] - 1] in the cut
 	   							pi_flowBound[w][i] = CPX[s].getDual(RFoutflowBound[s][0][i]); 
							}
							
						}
						else{   
							//for other stages than last stage
							alpha[w]=CPX[s].getValue(Alpha[s]);
							subObj_hat[w] = CPX[s].getObjValue()-alpha[w];
							
						}//end if
					}
				}//end scenario
			}//end stage 
			//end forward pass

			//Calculating P in a backward procedure
			for(s = No_stage-1; s>0; s--){
				for(w=0; w<No_scenNode; w++){
					if(stageMap[w]==s-1){
						
						for(j=0;j<No_scen;j++){	
							Objective_row.setLinearCoef(P[j], subObj_hat[decendant[w][j]]);
						}
						for(j=0;j<No_scen;j++){		
							conDistancePos[j].setUb(prob[decendant[w][j]]);
							conDistanceNeg[j].setUb(-prob[decendant[w][j]]);
						}
						cplex_row.solve();					
						
					
						for(j=0;j<No_scen;j++){
							p[decendant[w][j]][nCut] = cplex_row.getValue(P[j]);	
						}
						
						
						for(j = 0; j<No_scen; j++){
							subObj_hat[w] +=  p[decendant[w][j]][nCut] *subObj_hat[decendant[w][j]];
						}	
					}
				}
			}

			cout<<"calculating gap"<<endl;
				
			z_hat=subObj_hat[0];
			
			//*****Calculate gap****////
			if(z_hat <= UB || abs(z_hat-UB)/abs(UB)<= toler){
				UB = z_hat;
				for(w=0; w<No_scenNode; w++){
					Output[w].Cost=subObj_hat[w];
					Output[w].worst_Prob=p[w][nCut];
				}
			}

			
            
			
			if( (UB - LB)/LB < rel_Gap){
				rel_Gap = (UB - LB)/LB;
			}
			Result<< nCut+1 << '\t' << z_hat << '\t' << LB <<'\t' << UB <<'\t' << rel_Gap<<endl;
			// Convergence criterion
			if( rel_Gap <= toler ||nCut ==MaxCut-1 || float( clock () - begin_time ) /  CLOCKS_PER_SEC >7200){ 
			   cout<<"DONE"<<endl;
			   IloNum Obj=LB;
			   Result << float( clock () - begin_time ) /  CLOCKS_PER_SEC <<endl;
			   std::cout << "time =" << float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl; 
			   break;
			}


			////////////////////////// 

           for(s = No_stage-1; s>0; s--){
				MODEL[s-1].add(maxCut[s-1][nCut]);		
				for(j=0;j<No_scen;j++){	
					MODEL[s-1].add(optCut[s-1][j][nCut]);	
				}
				
                for(w=0; w<No_scenNode; w++){
					if(stageMap[w]==s-1){ 
						for(k=0;k<nCut;k++){ 
							maxCut[s-1][k].setBounds(0, IloInfinity);
							for(j=0;j<No_scen;j++){
								maxCut[s-1][k].setLinearCoef(Theta[s][j], -p[decendant[w][j]][k]);	 
							}

							for(j=0;j<No_scen;j++){
								optCut[s-1][j][k].setBounds(g[w][j][k], IloInfinity);	
								for(i=0;i<No_recharge; i++){ 
									optCut[s-1][j][k].setLinearCoef(Storage[i][s-1][No_year[s-1]-1], -G[w][j][k][i]);
								}
							}
						}//end for k
						
						
                       //you have done that before
                        if(s-1>0){
							for(t=0;t<No_year[s-1];t++){
								CapBound[s-1][t][0].setBounds(-IloInfinity, CAPamt[s-1][t]*scenFac[w][0]);
							   for(i=0;i<No_user;i++)
								   MeetDemand[s-1][t][i].setBounds(scenFac[w][t+1]*demand[s-1][t][i], scenFac[w][t+1]*demand[s-1][t][i]);
							   for(i=0;i<No_potUser; i++)
								   ReturnFlow[s-1][t][i].setBounds(returnRate*scenFac[w][t+1]*demand[s-1][t][i],returnRate*scenFac[w][t+1]*demand[s-1][t][i]);
							}

						
							for(w2=0; w2<No_scenNode; w2++){
								if(ancestor[w2][w]==1){
									for(i=0;i<No_recharge; i++){
										//update rhs using storage of previous stage
										StorageBalance[s-1][0][i].setBounds(-storage_hat[i][w2][nCut],-storage_hat[i][w2][nCut]);
										RFoutflowBound[s-1][0][i].setBounds(-IloInfinity, storage_hat[i][w2][nCut]);
									}
								}
							}
						}

						///Cut1
						//tage care of g= b + the rest of g
						//take care of b
			
						for(j=0;j<No_scen;j++){
							for(t =0; t<No_year[s];t++){
								for( i=0; i<No_user; i++){
									g[w][j][nCut]+=pi_demand[t][decendant[w][j]][i]*scenFac[decendant[w][j]][t+1]*demand[s][t][i];
								}
								for(i=0; i<No_potUser; i++){
									g[w][j][nCut] += pi_return[t][decendant[w][j]][i]*returnRate*scenFac[decendant[w][j]][t+1]*demand[s][t][i];
								}
								for(i=1;i<No_capacity; i++){ 
									g[w][j][nCut] += pi_capacity[t][decendant[w][j]][i]*capacity[i];
								}
								g[w][j][nCut] += pi_capacity[t][decendant[w][j]][0]*CAPamt[s][t]*scenFac[decendant[w][j]][0];
								for(i=0;i<No_recharge; i++){ 
									g[w][j][nCut] += pi_safe[t][decendant[w][j]][i]*storageUB[i];
								}
							}
						
						//take care of the rest of g
							if (s< No_stage-1){
								for(k=0;k<nCut+1;k++){ 
									for(w2 = 0; w2<No_scen; w2++){
										g[w][j][nCut] += pi_optCut[decendant[w][j]][w2][k]*g[decendant[w][j]][w2][k];	
									}
								}
							}
							optCut[s-1][j][nCut].setBounds(g[w][j][nCut], IloInfinity);

							//take care of G
							for(i=0;i<No_recharge; i++){
								G[w][j][nCut][i]=-pi_storage[decendant[w][j]][i]+ pi_flowBound[decendant[w][j]][i];
								optCut[s-1][j][nCut].setLinearCoef(Storage[i][s-1][No_year[s-1]-1], -G[w][j][nCut][i]); 
							}
						}//end j
						
						maxCut[s-1][nCut].setBounds(0, IloInfinity);
						for(j=0;j<No_scen;j++){				
							maxCut[s-1][nCut].setLinearCoef(Theta[s][j], -p[decendant[w][j]][nCut]);	 
						}
						
						if(s>1){
							CPX[s-1].solve();
							for(t=0;t<No_year[s-1];t++){ 
								for(i=0; i<No_capacity; i++) {
   									pi_capacity[t][w][i] = CPX[s-1].getDual(CapBound[s-1][t][i]); 
								}
								
								for(i=0; i<No_user; i++) {
   									pi_demand[t][w][i] = CPX[s-1].getDual(MeetDemand[s-1][t][i]); 
								}
								
								for(i=0; i<No_potUser; i++) {
   									pi_return[t][w][i] = CPX[s-1].getDual(ReturnFlow[s-1][t][i]); 
								}

								for(i=0; i<No_recharge; i++) 
 	   								pi_safe[t][w][i] = CPX[s-1].getDual(SafeYield[s-1][t][i]);
							}
							for(i=0; i<No_recharge; i++) {
   								pi_storage[w][i] = CPX[s-1].getDual(StorageBalance[s-1][0][i]); 
							}
						        
							for(i=0; i<No_recharge; i++) {
   								pi_flowBound[w][i] = CPX[s-1].getDual(RFoutflowBound[s-1][0][i]); 
							}

							//Dual value for optCut
							for(k=0;k<nCut+1;k++){
								for(j=0;j<No_scen;j++){
									pi_optCut[w][j][k] =  CPX[s-1].getDual(optCut[s-1][j][k]); 
								}
							}  
							 
						} 
					}
				}
			}//end backward  
			
		}//end main loop

		//put other stuff here
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

		std::sort(Output+1, Output+No_scenNode, CompareScenarios());
		
		Output[1].partial_sum=Output[1].nom_Prob;
		for(w=2; w<No_scenNode; w++){
			Output[w].partial_sum=Output[w-1].partial_sum+Output[w].nom_Prob;
		}

		vector<int>::iterator it;

		//STEP 1
		vector<int> Sup_index;
		Sup_index.push_back(No_scenNode-1);

		w=No_scenNode-2;
		while (abs(Output[w].Cost-Output[No_scenNode-1].Cost)<toler){
			Sup_index.push_back(w);
			w--;
		}
		const int noSup=Sup_index.size();		
		category[3]=noSup;


		//STEP 2 
		vector<int> VaR_index;
		if(Output[1].partial_sum >= rho )
			VaR_index.push_back(1);

		if (VaR_index.size()==0){
			for (w=2; w<No_scenNode+1; w++){
				if(Output[w].partial_sum+toler >= rho &&  Output[w-1].partial_sum < rho){
					VaR_index.push_back(w);
					break;
				}
			}
		}
		
		

		
		
		
		//check remval of scenarios with the same cost as index scenarios
		
		vector<int> removal;
		if (VaR_index[0]!=1){
			w=VaR_index[0]-1;
			while (abs(Output[VaR_index[0]].Cost-Output[w].Cost)<toler && w>=1 && w<No_scenNode){
				removal.push_back(w);
				w--;
			}
		}
		if (VaR_index[0]!=1)
			category[0]=w;
		else
			category[0]=0;
		if (category[0]!=0){
			for (int ww=1; ww<=w; ww++){
				Ineffscen[ii].push_back(ww);
				Output[ww].Status="Ineffective";
				category_effect[0]++;
			}
		}

		w=VaR_index[0];
		while (abs(Output[VaR_index[0]].Cost-Output[w].Cost)<toler && w>=1 && w<No_scenNode){
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
			category[2]=No_scenNode-(category[0]+category[1]+category[3])-1;
		}


		vector<int> third;

		//lambda>0
		if (abs(Output[VaR_index[0]].Cost-Output[Sup_index[0]].Cost)> toler){
			while (abs(Output[w].Cost - Output[Sup_index[0]].Cost)>toler && w>=1 && w<No_scenNode){
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
		Result << "w"<<'\t'<<"h"<<'\t'<<"q"<<'\t'<<"p"<<'\t'<<"cdf"<<'\t'<<"Status"<<endl;
		for(w=1; w<No_scenNode; w++){
			Result<< Output[w].No<< '\t';		
			Result<<Output[w].Cost<< '\t'<< Output[w].nom_Prob<< '\t' << Output[w].worst_Prob << '\t' << Output[w].partial_sum <<'\t';
			if (Output[w].Status !=NULL){
				Result<< Output[w].Status<<endl;	
			}
			else {
				Result<< endl;
			}

		}

		Result<<endl;
		

		IloNum temp_prob=0; 
		////check the the objective function after an effective scenario is removed
		for (it=Effscen[ii].begin(); it!=Effscen[ii].end(); it++){
			vector<int> temp_VaR_index;
			//ww=Output[*it].No;
			temp_prob=Output[*it].nom_Prob;
			IloNumArray partial_sum(env, No_scenNode);
			if (rho==0 && *it!=VaR_index[0]){
				temp_VaR_index.push_back(VaR_index[0]);
			}
			else if (*it==VaR_index[0]){
				temp_VaR_index.push_back(2);
			}
			
			if (rho>0){
				if (*it!=1){
					partial_sum[1]=Output[1].nom_Prob;
				}
				else{
					partial_sum[1]=0;
				}
				for(w=2; w<No_scenNode; w++){
					if (*it!=w){
						partial_sum[w]=partial_sum[w-1]+Output[w].nom_Prob;
					}
					else{
						partial_sum[w]=partial_sum[w-1];
					}
				}

				if(partial_sum[1] >= rho-temp_prob && *it !=1)
					temp_VaR_index.push_back(1);

				if (temp_VaR_index.size()==0){
					for (w=2; w<No_scenNode; w++){
						if(partial_sum[w] >= rho-temp_prob &&  partial_sum[w-1] < rho-temp_prob && *it!=w){
							temp_VaR_index.push_back(w);
							break;
						}
					}
				}
			}

			vector<int>::iterator match;
			IloNum newcost;
			match=std::find(Sup_index.begin(), Sup_index.end(), *it);
			if (match ==Sup_index.end()){
				newcost=rho*Output[No_scenNode-1].Cost;
			}
			else if (Sup_index.size()>1){
					newcost=rho*Output[No_scenNode-1].Cost;
			}
			else{
				newcost=rho*Output[No_scenNode-2].Cost;
			}
			newcost+=(1-rho)*Output[temp_VaR_index[0]].Cost;
			for(w=temp_VaR_index[0]+1; w<No_scenNode; w++){
				if (*it !=w){
					newcost+=Output[w].nom_Prob*(Output[w].Cost-Output[temp_VaR_index[0]].Cost);
				}
			}
			//store ratio of of the change in optimal value after removal of effective scenario
			NewCost_Effscen[ii].push_back(100*(LB-newcost)/LB);
		
			
		}

		
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



		cout<<"DONE"<<endl;
		






		struct CompareScenarios2
		{
			bool operator()(const Scenario& lhs, const Scenario& rhs)  const
			{
				return lhs.No < rhs.No;
			}
		};

		

		std::sort(Output+1, Output+No_scenNode, CompareScenarios2());
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

		for (w=1; w<No_scenNode; w++){
			if (Output[w].Status=="Effective"){
				Count_Effective[w]++;
			}
		}

		Numbers_post<< rho<<'\t';
		for (i=0; i<2; i++){
			Numbers_post<< category_effect[i]<< '\t';
		}
		Numbers_post<<endl;

		
		env.end();
		}//end loop on instances

		for (w=1; w<No_scenNode; w++){
			s=1;
			IloNum sumSupply=0;
			IloNum sumDemand=0;
			for(t=0;t<No_year[s];t++){
				sumSupply+=CAPamt[s][t]*scenFac[w][0];
		        for(i=0;i<No_user;i++){
					sumDemand+=scenFac[w][t+1]*demand[s][t][i];
				}
			}
			ScenarioRes<< w<< '\t';
			ScenarioRes<<sumSupply/No_year[s]<< '\t';
			ScenarioRes<<sumDemand/(No_year[s]*No_user)<< '\t';
			ScenarioRes<<Count_Effective[w]<<endl;
		}
		//vector<int>::iterator it;
		for (w=1; w<No_scenNode; w++){
			for (i=0; i<21; i++){
				if (NewCost_Effscen[i].size()>=w){
					EffScenarioRes<<NewCost_Effscen[i][w]<<'\t';
				}
				else{
					EffScenarioRes<<"0"<<'\t';
				}
			}
			EffScenarioRes<<endl;
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