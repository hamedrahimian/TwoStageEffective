1. "config.txt" contains 21 gamma parameters.
2. "DemPortion.txt" contains demand values over years (rows) for each portable user (columns).
3. "input2.txt" 
	a. First four rows (each row correspond to one stage) contains population of the studied area over years (columns).
	b. Last four rows (each row correspond to one stage) contains Tucson population over years (columns).
4. "input1_%d.txt" in the order of rows
	a. Number of years in each stage (column).
	b. Type of the node.
	c. Start node of the arc.
	d. End node of the arc.
	e. Parameter c_ij.
	f. Parameter H. The first element corresponds to SW, the rest correspond to TP or RF.
	g. Parameter a.
	h. Parameter l_ij.
	i. Parameter S^0.
	j. Parameter U.
5. "scenFac_%d.txt" contain scenario data.