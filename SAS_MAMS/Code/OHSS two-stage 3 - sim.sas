data sim;

  seed = 2655780;
  rep = 1000000;

  n = 27;
  n11 = 2*n;   n21 = n;   n31 = n;
  n12 = 4*n;   n22 = 2*n;   n32 = 2*n;
  p1 = 0.7;   p2 = 0.7;   p3 = 0.7;
  l1 = probit(0.27);
  u2 = 1.92134;

  do i = 1 to rep;

	s11 = ranbin(seed, n11, p1);
    s12 = s11 + ranbin(seed, n11, p1);
    s21 = ranbin(seed, n21, p2);
    s22 = s21 + ranbin(seed, n21, p2);
	s31 = ranbin(seed, n31, p3);
    s32 = s31 + ranbin(seed, n31, p3);

    z121 = (n21*s11 - n11*s21)/(n11 + n21);
    v121 = n11*n21*(s11 + s21)*(n11 + n21 - s11 - s21)/((n11 + n21)**3);
    z122 = (n22*s12 - n12*s22)/(n12 + n22);
    v122 = n12*n22*(s12 + s22)*(n12 + n22 - s12 - s22)/((n12 + n22)**3);
	z131 = (n31*s11 - n11*s31)/(n11 + n31);
    v131 = n11*n31*(s11 + s31)*(n11 + n31 - s11 - s31)/((n11 + n31)**3);
    z132 = (n32*s12 - n12*s32)/(n12 + n32);
    v132 = n12*n32*(s12 + s32)*(n12 + n32 - s12 - s32)/((n12 + n32)**3);
	z231 = (n31*s21 - n21*s31)/(n21 + n31);
    v231 = n21*n31*(s21 + s31)*(n21 + n31 - s21 - s31)/((n21 + n31)**3);
    z232 = (n32*s22 - n22*s32)/(n22 + n32);
    v232 = n22*n32*(s22 + s32)*(n22 + n32 - s22 - s32)/((n22 + n32)**3);

    elim21 = (z121/sqrt(v121) gt l1);
	elim31 = (z131/sqrt(v131) gt l1);
	stop1 = elim21*elim31;
	cont121 = elim31*(1 - elim21);
	cont131 = elim21*(1 - elim31);
	contall = (1 - elim21)*(1 - elim31);
	cont1 = cont121 + cont131 + contall;

   	choose2 = (cont121 + contall)*(z122/sqrt(v122) lt -u2);
	choose3 = (cont131 + contall)*(z132/sqrt(v132) lt -u2);
	choose = 1 - (1 - choose2)*(1 - choose3);

	nobs = 4*n*(1 + contall) + 3*n*(cont121 + cont131);

    output;
  end;

run;

proc means data = sim mean;
run;
