data read;
  do n = 20 to 30;
    l1 = -probit(0.27);
	u2 = 1.92134;
	corr = sqrt(0.5);

	p1 = 0.7;
	p2 = 0.9;
	theta = -log((p1*(1 - p2))/(p2*(1 - p1)));
	pbar = (2*p1 + p2)/3;
	v1 = 2*n*pbar*(1 - pbar)/3;
	v2 = 2*v1;
	delta1 = theta*sqrt(v1);
	delta2 = theta*sqrt(v2);

	pi1 = probnorm(l1 - delta1);
	pproc1 = probbnrm(-l1 + delta1, -u2 + delta2, corr);
	output;
  end;
run;

proc print data = read;
run;
