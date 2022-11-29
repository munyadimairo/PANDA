%macro interim;
data interim;
  set input;

  z12 = (n2*s1 - n1*s2)/(n1 + n2);
  z13 = (n3*s1 - n1*s3)/(n1 + n3);
  z23 = (n3*s2 - n2*s3)/(n2 + n3);
  v12 = n1*n2*(s1 + s2)*(n1 + n2 - s1 - s2)/((n1 + n2)**3);
  v13 = n1*n3*(s1 + s3)*(n1 + n3 - s1 - s3)/((n1 + n3)**3);
  v23 = n2*n3*(s2 + s3)*(n2 + n3 - s2 - s3)/((n2 + n3)**3);

  p1hat = (s1/n1);	
  p2hat = (s2/n2);	
  p3hat = (s3/n3);	
  theta12hat = (z12/v12);
  theta13hat = (z13/v13);
  theta23hat = (z23/v23);

  elim2 = (z12/sqrt(v12) gt u);
  elim3 = (z13/sqrt(v13) gt u);
  stop1 = elim2*elim3;
  cont12 = elim3*(1 - elim2);
  cont13 = elim2*(1 - elim3);
  contall = (1 - elim2)*(1 - elim3);
  cont1 = cont12 + cont13 + contall;

%mend;

%macro analyse;

data analyse;
  set interim;

  varp1hat = p1hat*(1 - p1hat)/n1;
  varp2hat = p2hat*(1 - p2hat)/n2;
  varp3hat = p3hat*(1 - p3hat)/n3;
  vartheta12hat = 1/v12;
  vartheta13hat = 1/v13;
  vartheta23hat = 1/v23;

  p1lo = p1hat - 1.96*sqrt(varp1hat);
  p2lo = p2hat - 1.96*sqrt(varp2hat);
  p3lo = p3hat - 1.96*sqrt(varp3hat);
  p1hi = p1hat + 1.96*sqrt(varp1hat);
  p2hi = p2hat + 1.96*sqrt(varp2hat);
  p3hi = p3hat + 1.96*sqrt(varp3hat);

  theta12lo = theta12hat - 1.96*sqrt(vartheta12hat);
  theta13lo = theta13hat - 1.96*sqrt(vartheta13hat);
  theta23lo = theta23hat - 1.96*sqrt(vartheta23hat);
  theta12hi = theta12hat + 1.96*sqrt(vartheta12hat);
  theta13hi = theta13hat + 1.96*sqrt(vartheta13hat);
  theta23hi = theta23hat + 1.96*sqrt(vartheta23hat);

run;

%mend;

%macro rbanalyse;

data input;
  set results;

  srb_11_lo = max(0, s12 - (n12 - n11));
  srb_11_hi = min(s12, n11);
  srb_21_lo = max(0, s22 - (n22 - n21));
  srb_21_hi = min(s22, n21);
  srb_31_lo = max(0, s32 - (n32 - n31));
  srb_31_hi = min(s32, n31);

  n1 = n11;   n2 = n21;   n3 = n31;
  u = l1;

  do s1 = srb_11_lo to srb_11_hi;
    do s2 = srb_21_lo to srb_21_hi;
	  do s3 = srb_31_lo to srb_31_hi;

	    prob1 = pdf('hyper', s1, n12, s12, n11, 1);
	    prob2 = pdf('hyper', s2, n22, s22, n21, 1);
		prob3 = pdf('hyper', s3, n32, s32, n31, 1);
	    prob = prob1*prob2*prob3;  

        cont121 = (n22 > n21)*(n32 = n31); 
		cont131 = (n22 = n21)*(n32 > n31); 
		contall1 = (n22 > n21)*(n32 > n31); 
 
	    output;
	  end;
    end;
  end;

run;

%mend;

%macro onesim;

data repno;
  set repno;
  rep = rep + 1;
run;

data results;
  set repno;
  seed = 88697453 + rep;
  seed = int(1000000*ranuni(seed));
  l1 = probit(0.27);   u2 = 1.92134;
  p1 = 0.7;   p2 = 0.9;   p3 = 0.9;
  n11 = 54;	   n21 = 27;	n31 = 27;
  n12 = 108;   n22 = 54;    n32 = 54;
  s11 = ranbin(seed, n11, p1);
  s21 = ranbin(seed, n21, p2);
  s31 = ranbin(seed, n31, p3);
  s12 = s11 + ranbin(seed, (n12 - n11), p1);
  s22 = s21 + ranbin(seed, (n22 - n21), p2);
  s32 = s31 + ranbin(seed, (n32 - n31), p3);
run;

data input;
  set results;
  u = l1;
  n1 = n11;   n2 = n21;   n3 = n31;
  s1 = s11;   s2 = s21;   s3 = s31;
run;

%interim;

%analyse;

data results;
  set interim;
  n12 = (stop1 = 0)*n12 + (stop1 = 1)*n11;
  n22 = (elim2 = 0)*n22 + (elim2 = 1)*n21;
  n32 = (elim3 = 0)*n32 + (elim3 = 1)*n31;
  s12 = (stop1 = 0)*s12 + (stop1 = 1)*s11;
  s22 = (elim2 = 0)*s22 + (elim2 = 1)*s21;
  s32 = (elim3 = 0)*s32 + (elim3 = 1)*s31;
run;

%rbanalyse;

%interim;

data valid;
  set interim;
  valid = (cont12 = cont121)*(cont13 = cont131)*(contall = contall1);
  if valid = 0 then delete;
run;

data rb;
  set valid;
  one = 1;

  p1hat = p1hat*prob;
  p2hat = p2hat*prob;
  p3hat = p3hat*prob;
  theta12hat = theta12hat*prob;
  theta13hat = theta13hat*prob;
  theta23hat = theta23hat*prob;

  p1hatsq = p1hat*p1hat/prob;
  p2hatsq = p2hat*p2hat/prob;
  p3hatsq = p3hat*p3hat/prob;
  theta12hatsq = theta12hat*theta12hat/prob;
  theta13hatsq = theta13hat*theta13hat/prob;
  theta23hatsq = theta23hat*theta23hat/prob;

run;

proc means data = rb sum noprint;
  var one prob1 prob2 prob3 prob p1hat p2hat p3hat theta12hat theta13hat theta23hat
      p1hatsq p2hatsq p3hatsq theta12hatsq theta13hatsq theta23hatsq;
  output out = rbmeans 
	sum = sumone sumprob1 sumprob2 sumprob3 sumprob sump1hat sump2hat  sump3hat 
          sumtheta12hat sumtheta13hat sumtheta23hat sump1hatsq sump2hatsq  sump3hatsq 
          sumtheta12hatsq sumtheta13hatsq sumtheta23hatsq;
run;

data rbresults;
  merge rbmeans analyse;
  prb1hat = sump1hat/sumprob;
  prb2hat = sump2hat/sumprob;
  prb3hat = sump3hat/sumprob;
  thetarb12hat = sumtheta12hat/sumprob;
  thetarb13hat = sumtheta13hat/sumprob;
  thetarb23hat = sumtheta23hat/sumprob;

  prb1hatsq = sump1hatsq/sumprob;
  prb2hatsq = sump2hatsq/sumprob;
  prb3hatsq = sump3hatsq/sumprob;
  thetarb12hatsq = sumtheta12hatsq/sumprob;
  thetarb13hatsq = sumtheta13hatsq/sumprob;
  thetarb23hatsq = sumtheta23hatsq/sumprob;

  varrbp1hat = varp1hat - (prb1hatsq - prb1hat*prb1hat);
  varrbp2hat = varp2hat - (prb2hatsq - prb2hat*prb2hat);
  varrbp3hat = varp3hat - (prb3hatsq - prb3hat*prb3hat);
  varthetarb12hat = vartheta12hat - (thetarb12hatsq - thetarb12hat*thetarb12hat);
  varthetarb13hat = vartheta13hat - (thetarb13hatsq - thetarb13hat*thetarb13hat);
  varthetarb23hat = vartheta23hat - (thetarb23hatsq - thetarb23hat*thetarb23hat);

  prb1lo = prb1hat - 1.96*sqrt(varrbp1hat);
  prb2lo = prb2hat - 1.96*sqrt(varrbp2hat);
  prb3lo = prb3hat - 1.96*sqrt(varrbp3hat);
  prb1hi = prb1hat + 1.96*sqrt(varrbp1hat);
  prb2hi = prb2hat + 1.96*sqrt(varrbp2hat);
  prb3hi = prb3hat + 1.96*sqrt(varrbp3hat);

  thetarb12lo = thetarb12hat - 1.96*sqrt(varthetarb12hat);
  thetarb13lo = thetarb13hat - 1.96*sqrt(varthetarb13hat);
  thetarb23lo = thetarb23hat - 1.96*sqrt(varthetarb23hat);
  thetarb12hi = thetarb12hat + 1.96*sqrt(varthetarb12hat);
  thetarb13hi = thetarb13hat + 1.96*sqrt(varthetarb13hat);
  thetarb23hi = thetarb23hat + 1.96*sqrt(varthetarb23hat);

run;

data input;
  set results;
  u = u2;
  n1 = n12;   n2 = n22;   n3 = n32;
  s1 = s12;   s2 = s22;   s3 = s32;
run;

%interim;

%analyse;

data naive;
  set analyse;
  pn1lo = p1lo;   pn1hat = p1hat;   pn1hi = p1hi;
  pn2lo = p2lo;   pn2hat = p2hat;   pn2hi = p2hi;
  pn3lo = p3lo;   pn3hat = p3hat;   pn3hi = p3hi; 
  thetan12lo = theta12lo;   thetan12hat = theta12hat;   thetan12hi = theta12hi; 
  thetan13lo = theta13lo;   thetan13hat = theta13hat;   thetan13hi = theta13hi;
  thetan23lo = theta23lo;   thetan23hat = theta23hat;   thetan23hi = theta23hi;
  drop p1lo p1hat p1hi p2lo p2hat p2hi p3lo p3hat p3hi;
  drop theta12lo theta12hat theta12hi theta13lo theta13hat theta13hi theta23lo theta23hat theta23hi;
run;

data allres;
  merge rbresults naive;

  thetarbop12lo  = (1 - (n12 > n11)*(n22 > n21))*theta12lo  + (n12 > n11)*(n22 > n21)*thetarb12lo;
  thetarbop12hat = (1 - (n12 > n11)*(n22 > n21))*theta12hat + (n12 > n11)*(n22 > n21)*thetarb12hat;
  thetarbop12hi  = (1 - (n12 > n11)*(n22 > n21))*theta12hi  + (n12 > n11)*(n22 > n21)*thetarb12hi;
  thetarbop13lo  = (1 - (n12 > n11)*(n32 > n31))*theta13lo  + (n12 > n11)*(n32 > n31)*thetarb13lo;
  thetarbop13hat = (1 - (n12 > n11)*(n32 > n31))*theta13hat + (n12 > n11)*(n32 > n31)*thetarb13hat;
  thetarbop13hi  = (1 - (n12 > n11)*(n32 > n31))*theta13hi  + (n12 > n11)*(n32 > n31)*thetarb13hi;
  thetarbop23lo  = (1 - (n22 > n21)*(n32 > n31))*theta23lo  + (n22 > n21)*(n32 > n31)*thetarb23lo;
  thetarbop23hat = (1 - (n22 > n21)*(n32 > n31))*theta23hat + (n22 > n21)*(n32 > n31)*thetarb23hat;
  thetarbop23hi  = (1 - (n22 > n21)*(n32 > n31))*theta23hi  + (n22 > n21)*(n32 > n31)*thetarb23hi;

  thetanop12lo  = (1 - (n12 > n11)*(n22 > n21))*theta12lo  + (n12 > n11)*(n22 > n21)*thetan12lo;
  thetanop12hat = (1 - (n12 > n11)*(n22 > n21))*theta12hat + (n12 > n11)*(n22 > n21)*thetan12hat;
  thetanop12hi  = (1 - (n12 > n11)*(n22 > n21))*theta12hi  + (n12 > n11)*(n22 > n21)*thetan12hi;
  thetanop13lo  = (1 - (n12 > n11)*(n32 > n31))*theta13lo  + (n12 > n11)*(n32 > n31)*thetan13lo;
  thetanop13hat = (1 - (n12 > n11)*(n32 > n31))*theta13hat + (n12 > n11)*(n32 > n31)*thetan13hat;
  thetanop13hi  = (1 - (n12 > n11)*(n32 > n31))*theta13hi  + (n12 > n11)*(n32 > n31)*thetan13hi;
  thetanop23lo  = (1 - (n22 > n21)*(n32 > n31))*theta23lo  + (n22 > n21)*(n32 > n31)*thetan23lo;
  thetanop23hat = (1 - (n22 > n21)*(n32 > n31))*theta23hat + (n22 > n21)*(n32 > n31)*thetan23hat;
  thetanop23hi  = (1 - (n22 > n21)*(n32 > n31))*theta23hi  + (n22 > n21)*(n32 > n31)*thetan23hi;

  keep 	rep p1 p2 p3 n11 n21 n31 n12 n22 n32
  		p1lo p1hat p1hi p2lo p2hat p2hi p3lo p3hat p3hi
  		theta12lo theta12hat theta12hi theta13lo theta13hat theta13hi theta23lo theta23hat theta23hi
  		prb1lo prb1hat prb1hi prb2lo prb2hat prb2hi prb3lo prb3hat prb3hi
  		thetarb12lo thetarb12hat thetarb12hi thetarb13lo thetarb13hat thetarb13hi thetarb23lo thetarb23hat thetarb23hi
  		pn1lo pn1hat pn1hi pn2lo pn2hat pn2hi pn3lo pn3hat pn3hi
  		thetan12lo thetan12hat thetan12hi thetan13lo thetan13hat thetan13hi thetan23lo thetan23hat thetan23hi
		thetarbop12lo thetarbop12hat thetarbop12hi thetarbop13lo thetarbop13hat thetarbop13hi thetarbop23lo thetarbop23hat thetarbop23hi
  		thetanop12lo thetanop12hat thetanop12hi thetanop13lo thetanop13hat thetanop13hi thetanop23lo thetanop23hat thetanop23hi;

run;

data resfile;
  set allres resfile;
run;

%mend;

%macro simrun;

%do simulation = 1 %to 10000;
  %onesim;
%end;

%mend;

proc printto log = _null_;
run;

data repno;
  rep = 0;
run;

data resfile;
  blank = 0;
run;

%simrun;

data resfile;
  set resfile;
  if blank = 0 then delete;
run;

data cover;
  set resfile;

  theta12 = log(p1/(1 - p1)) - log(p2/(1 - p2));
  theta13 = log(p1/(1 - p1)) - log(p3/(1 - p3));
  theta23 = log(p2/(1 - p2)) - log(p3/(1 - p3));

  incp1ci = (p1 > p1lo)*(p1 < p1hi);
  incp2ci = (p2 > p2lo)*(p2 < p2hi);
  incp3ci = (p3 > p3lo)*(p3 < p3hi);
  incprb1ci = (p1 > prb1lo)*(p1 < prb1hi);
  incprb2ci = (p2 > prb2lo)*(p2 < prb2hi);
  incprb3ci = (p3 > prb3lo)*(p3 < prb3hi);
  incpn1ci = (p1 > pn1lo)*(p1 < pn1hi);
  incpn2ci = (p2 > pn2lo)*(p2 < pn2hi);
  incpn3ci = (p3 > pn3lo)*(p3 < pn3hi);

  inctheta12ci = (theta12 > theta12lo)*(theta12 < theta12hi);
  inctheta13ci = (theta13 > theta13lo)*(theta13 < theta13hi);
  inctheta23ci = (theta23 > theta23lo)*(theta23 < theta23hi);
  incthetarb12ci = (theta12 > thetarb12lo)*(theta12 < thetarb12hi);
  incthetarb13ci = (theta13 > thetarb13lo)*(theta13 < thetarb13hi);
  incthetarb23ci = (theta23 > thetarb23lo)*(theta23 < thetarb23hi);
  incthetan12ci = (theta12 > thetan12lo)*(theta12 < thetan12hi);
  incthetan13ci = (theta13 > thetan13lo)*(theta13 < thetan13hi);
  incthetan23ci = (theta23 > thetan23lo)*(theta23 < thetan23hi);
  incthetarbop12ci = (theta12 > thetarbop12lo)*(theta12 < thetarbop12hi);
  incthetarbop13ci = (theta13 > thetarbop13lo)*(theta13 < thetarbop13hi);
  incthetarbop23ci = (theta23 > thetarbop23lo)*(theta23 < thetarbop23hi);
  incthetanop12ci = (theta12 > thetanop12lo)*(theta12 < thetanop12hi);
  incthetanop13ci = (theta13 > thetanop13lo)*(theta13 < thetanop13hi);
  incthetanop23ci = (theta23 > thetanop23lo)*(theta23 < thetanop23hi);

  tot = n12 + n22 + n32;

  drop 	blank p1lo p1hi p2lo p2hi p3lo p3hi
  		theta12lo theta12hi theta13lo theta13hi theta23lo theta23hi
  		prb1lo prb1hi prb2lo prb2hi prb3lo prb3hi
  		thetarb12lo thetarb12hi thetarb13lo thetarb13hi thetarb23lo thetarb23hi
  		pn1lo pn1hi pn2lo pn2hi pn3lo pn3hi
  		thetan12lo thetan12hi thetan13lo thetan13hi thetan23lo thetan23hi
		thetarbop12lo thetarbop12hi thetarbop13lo thetarbop13hi thetarbop23lo thetarbop23hi
  		thetanop12lo thetanop12hi thetanop13lo thetanop13hi thetanop23lo thetanop23hi;

run;

proc printto;
run;

proc means data = cover mean;
run;







