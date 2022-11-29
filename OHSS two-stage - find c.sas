data read;
  do i = 0 to 10;
	l1 = -probit(0.27);
	u2 = 1.9213 + 0.00001*i;
	corr = sqrt(0.5);
	pi0 = probnorm(l1);
	pproc0 = probbnrm(-l1, -u2, corr);
	output;
  end;
run;

proc print data = read;
run;
