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

data sim;

  seed = 2655780;
  rep = 1000000;
  n = 27;
  n11 = 2*n;
  n21 = n;
  n12 = 4*n;
  n22 = 2*n;
  p1 = 0.7;
  p2 = 0.7;
  l1 = -probit(0.27);
  u2 = 1.92134;

  do i = 1 to rep;
    s11 = ranbin(seed, n11, p1);
    s12 = s11 + ranbin(seed, n11, p1);
    s21 = ranbin(seed, n21, p2);
    s22 = s21 + ranbin(seed, n21, p2);

    z1 = (n21*s11 - n11*s21)/(n11 + n21);
    v1 = n11*n21*(s11 + s21)*(n11 + n21 - s11 - s21)/((n11 + n21)**3);
    z2 = (n22*s12 - n12*s22)/(n12 + n22);
    v2 = n12*n22*(s12 + s22)*(n12 + n22 - s12 - s22)/((n12 + n22)**3);

    stop1 = (z1/sqrt(v1) gt -l1);
    cont1 = 1 - stop1;
    choose2 = cont1*(z2/sqrt(v2) lt -u2);
    nodiff = 1 - choose2;
	nobs = 3*n*(1 + cont1);

    output;
  end;

run;

proc means data = sim mean;
  var rep stop1 cont1 choose2 nodiff nobs;
run;
