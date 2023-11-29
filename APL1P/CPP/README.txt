"APL1P.cpp" (Table 1)
	a. The inputs are files "inputs.txt" and "config.txt".  
	b. The outputs are 
		1. "%.2f.txt", where f denote the value of gamma. This includes the log of the algorithm iterations, CPU time, and information about the effectiveness of scenarios.
		2. "Numbers.txt". Number of scenairos in each primal category, number of effective, ineffective, and undetermined scenarios by the easy-to-check conditions for each gamma.
		3. "Numbers_post.txt". Number of effective and ineffective scenarios after resolving undetermined scenarios.
		4. "Average.txt". The number of times each sceanrio is effective by easy-to-check conditions out of 21 gammas.
		5. "EX_nestedness.txt". Whether the effective scenarios are nested or not. 
-------------------------------------------------------------------------------

"APL1P_modified.cpp" (Table 2b, Example 4)
	a. The inputs are files "inputs_modified.txt" and "config.txt".  
	b. The outputs are 
		1. "EX_%.2f.txt", where f denote the value of gamma. This includes the log of the algorithms iterations, CPU time, and information about the effectiveness of scenarios.		
		2. "EX_Numbers.txt". Number of scenairos in each primal category, number of effective, ineffective, and undetermined scenarios by the easy-to-check conditions for each gamma.
		3. "EX_Numbers_post.txt". Number of effective and ineffective sceanrios after resolving undetermined scenarios.
		4. "EX_Average.txt". The number of times each sceanrio is effective by easy-to-check conditions out of 21 gammas.
		5. "EX_nestedness.txt". Whether the effective scenarios are nested or not. 		
		6. "EX_Status.txt". Contains the effectivenss of each scenario for each value of gamma.
		7. "EX_Details.txt". Contains optimal solution, the scenarios in the second category, the scenarios in the fourth catgory, and set of effective scenarios for each value of gamma.
		