1. "config.txt" contains 21 gamma parameters.
2. "inputs.txt" in the order of rows
	a. min capacity for generator i.
	b. max capacity for generator i.
	c. investment cost per unit capacity  for generator i.
	d. operating cost per unit at operation level j (column) on generator i (row).
	e. cost per unit of unmet demand at operation level j.
	f. probabilities (the first two rows correspond to generators availability, the next three rows correspond to demand).
	g. availability of generator i.
	h. customer demand for operation level j.

3. "inputs_modified.txt" is a modified version of "inputs.txt" where only the first generator and the second demand are stochastic, and the other random paramteres are fixed at their expected values
	a. min capacity for generator i.
	b. max capacity for generator i.
	c. investment cost per unit capacity  for generator i.
	d. operating cost per unit at operation level j (column) on generator i (row).
	e. cost per unit of unmet demand at operation level j.
	f. probabilities (the first two rows correspond to generators availability, the next three rows correspond to demand).
	g. availability of generator 1.
	h. customer demand 2 for operation level j.
	i. availability of generator 2.
	j. customer demand 1 and 3.

