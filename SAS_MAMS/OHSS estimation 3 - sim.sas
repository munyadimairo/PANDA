data sim;

  seed = 2655780;
  rep = 1000000;
  n = 27;
  n11 = 2*n;   n21 = n;   n31 = n;
  n12 = 4*n;   n22 = 2*n;   n32 = 2*n;
  p1 = 0.7;   p2 = 0.7;   p3 = 0.7; 
  theta12 = log(p1/(1 - p1)) - log(p2/(1 - p2));
  theta13 = log(p1/(1 - p1)) - log(p3/(1 - p3));
  theta23 = log(p2/(1 - p2)) - log(p3/(1 - p3));
  l1 = -probit(0.27);
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

    elim21 = (z121/sqrt(v121) gt -l1);
	elim31 = (z131/sqrt(v131) gt -l1);
	stop1 = elim21*elim31;
	cont121 = elim31*(1 - elim21);
	cont131 = elim21*(1 - elim31);
	contall = (1 - elim21)*(1 - elim31);
	cont1 = cont121 + cont131 + contall;

   	choose2 = (cont121 + contall)*(z122/sqrt(v122) lt -u2);
	choose3 = (cont131 + contall)*(z132/sqrt(v132) lt -u2);
	choose = 1 - (1 - choose2)*(1 - choose3);

	n1 = n11*stop1 + n12*(cont121 + cont131 + contall);
	n2 = n21*(stop1 + cont131) + n22*(cont121 + contall);
	n3 = n31*(stop1 + cont121) + n32*(cont131 + contall);

	p1hat = (s11/n11)*stop1 + (s12/n12)*(cont121 + cont131 + contall);	
	p2hat = (s21/n21)*(stop1 + cont131) + (s22/n22)*(cont121 + contall);	
	p3hat = (s31/n31)*(stop1 + cont121) + (s32/n32)*(cont131 + contall);	
	theta12hat = (z121/v121)*(stop1 + cont131) + (z122/v122)*(cont121 + contall);
	theta13hat = (z131/v131)*(stop1 + cont121) + (z132/v132)*(cont131 + contall);
	theta23hat = (z231/v231)*(stop1 + cont121 + cont131) + (z232/v232)*contall;
	var12 = (1/v121)*(stop1 + cont131) + (1/v122)*(cont121 + contall);
	var13 = (1/v131)*(stop1 + cont121) + (1/v132)*(cont131 + contall);
	var23 = (1/v231)*(stop1 + cont121 + cont131) + (1/v232)*contall;

	p1lo = p1hat - 1.96*sqrt((p1hat*(1 - p1hat))/n1);
	p2lo = p2hat - 1.96*sqrt((p2hat*(1 - p2hat))/n2);
	p3lo = p3hat - 1.96*sqrt((p3hat*(1 - p3hat))/n3);
	p1hi = p1hat + 1.96*sqrt((p1hat*(1 - p1hat))/n1);
	p2hi = p2hat + 1.96*sqrt((p2hat*(1 - p2hat))/n2);
	p3hi = p3hat + 1.96*sqrt((p3hat*(1 - p3hat))/n3);
	incp1 = (p1 ge p1lo)*(p1 le p1hi);
	incp2 = (p2 ge p2lo)*(p2 le p2hi);
	incp3 = (p3 ge p3lo)*(p3 le p3hi);

    theta12lo = theta12hat - 1.96*sqrt(var12);
	theta13lo = theta13hat - 1.96*sqrt(var13);
    theta23lo = theta23hat - 1.96*sqrt(var23);
	theta12hi = theta12hat + 1.96*sqrt(var12);
	theta13hi = theta13hat + 1.96*sqrt(var13);
    theta23hi = theta23hat + 1.96*sqrt(var23);
	inctheta12 = (theta12 ge theta12lo)*(theta12 le theta12hi);
	inctheta13 = (theta13 ge theta13lo)*(theta13 le theta13hi);
	inctheta23 = (theta23 ge theta23lo)*(theta23 le theta23hi);

    output;
  end;

run;

proc means data = sim mean;
run;
