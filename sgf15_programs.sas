/***********************************************************
	
   Programs used in SGF 2015 paper #2440
   "Permit Me to Permute: A Basic Introduction to Permutation Tests with SAS/IML" 
   http://support.sas.com/resources/papers/proceedings15/2440-2015.pdf

************************************************************/

/* create dataset */
data permtest;
   input prod1-prod3 @@;
datalines;
0 0 6       8 3 8         3 4 4         10 10 10     2 6 4        7 7 10
10 10 10    0 0 0         10 10 10      3 4 9        2 4 5        9 10 10
2 0 0       10 10 10      2 2 0         3 3 9        4 4 7        2 1 10
3 0 4       2 10 10       6 5 7         7 7 0        5 4 4        2 1 0
0 10 3      9 10 9        1 1 4         6 6 8        10 10 10     4 3 6
10 10 9     8 7 4         3 7 8         4 1 1        4 5 6        4 4 7
2 4 3       0 2 3         9 1 6         10 10 10     7 0 9        10 10 10
0 0 3       1 1 5         10 0 10       6 9 1        4 3 5        10 4 10
10 6 10     6 2 6         9 9 7         10 10 10     10 10 10     1 4 6
7 1 7       9 9 0         7 4 5         4 7 7        1 0 5        5 3 6
7 2 1       6 2 6         3 7 4         0 1 1        7 10 10      3 4 6
1 0 10      3 2 10        10 2 10       6 1 5        1 0 5        0 0 0
0 0 4       7 8 10        10 6 10       3 10 10      5 5 5        5 0 0
0 0 1       4 5 7         0 0 0         0 0 1        9 10 5       9 10 8
7 0 8       10 10 10      0 3 1         3 0 6        9 7 6        0 6 4
0 0 0       10 10 10      1 2 4         7 8 9        9 9 10       8 9 9
8 9 7       2 0 0         10 9 10       6 8 10       0 8 7        9 10 10
8 10 10     9 10 10       10 10 10      1 1 8        9 8 8        1 2 9
10 10 0     1 3 8         0 0 0         4 2 2        10 10 10     8 10 2
3 5 7       0 2 1         5 7 2         1 3 6        2 1 3        0 0 0
10 10 10    1 0 4         3 4 0         0 1 0        4 7 7        2 2 1
0 0 5       0 0 0         4 5 6         2 3 4        4 5 3        0 0 8
3 7 5       8 10 0        6 3 2         3 0 2        8 10 3       10 10 4
0 0 3
;
run;


/* Matched Pairs example */
proc iml;
   varnames = {'prod1' 'prod2'};   
   use permtest;    
   read all var varnames;
   close permtest;            

   prod12 = prod1-prod2;

   ObsDiff = prod12[:];
   call symputx('obsdiff',ObsDiff); 

   call randseed(1234);             /* initialize seed for random numbers */
   y = prod1 || prod2;              /* concatenate values into n x 2 matrix */

/* METHOD #1 - USING LOC FUNCTION */
/* see page 11 and 12 of Wicklin's "Rediscovering SAS/IML Software: Modern Data Analysis for the Practicing Statistician" */
/* http://support.sas.com/resources/papers/proceedings10/329-2010.pdf */
   u = j(nrow(y),1);                /* allocate vector to hold random numbers */
   B = 10000;                       /* number of bootstrap resamples */
   s = j(B,1);                      /* allocate vector to hold bootstrap dist */
   do i = 1 to B;
      call randgen(u, "Uniform");
      z = y;                       /* copy original data */
      rows = loc(u > 0.5);         /* locate rows for which u[i] > 0.5 */
      z[rows,] = z[rows,{2 1}];    /* swap values of these rows */
      x = z[,1] - z[,2];           /* difference of permuted data */
      s[i] = x[:];                 /* mean of difference */
   end;

/* compute empirical two-sided p-value */
   pval = sum(s > abs(ObsDiff)) / B + sum(s < -abs(ObsDiff)) / B;
   print pval[label='PROD1-PROD2 P Value - USING LOC FUNCTION'];

/* METHOD #2 - USING RANPERM FUNCTION AND PERMUTEWITHINROWS MODULE */
/* module to independently permute elements of each row of a matrix */
/* see SAS Communities: https://communities.sas.com/message/212330 */
/* see "The DO Loop": http://blogs.sas.com/content/iml/2014/05/29/permute-elements-within-each-row-of-a-matrix.html */
   start PermuteWithinRows(m);
      colIdx = ranperm(1:ncol(m), nrow(m));
      f = (row(m) - 1) * ncol(m);               
      matIdx = f + colIdx;                    
      return( shape(m[matIdx], nrow(m)) );
   finish;

   B = 10000;                                 
   mpdist = j(B,1);
   do i = 1 to B;  
      x = PermuteWithinRows(y);
      x = x[,1] - x[,2];
      mpdist[i,] = x[:];
   end;

   create bootmp var {mpdist};     /* create dataset from mp vector*/ 
   append;
   close bootmp;
   
/* compute empirical two-sided p-value */
   pval = sum(mpdist>abs(ObsDiff)) / B + sum(mpdist<-abs(ObsDiff)) / B;
   print pval[label='PROD1-PROD2 P Value - USING RANPERM AND MODULE'];
   call symputx('prod12p',pval);

quit;



/* Repeated Measures ANOVA example */
proc iml;
   use permtest;
   read all var{prod1 prod2 prod3} into xobt;
   close permtest;

/* module to calculate F ratio */
   start fmod(x);
      gmean = x[:];                         /* grand mean scalar */
      n = nrow(x);                          /* number of rows/subjects */
      k = ncol(x);                          /* number of cols. */
      tmeans = x[:,];                       /* between means */
      submeans = x[,:];                     /* subject means */   
      /* SS between */
      SSb = (tmeans-gmean)##2;
      SSb = n*(SSb[+]);
      /* SS within */
      SSw = (tmeans-x)##2;
      SSw = SSw[+,];
      SSw = SSw[+];
      /* SS subjects */
      SSsub = k*((submeans-gmean)##2);
      SSsub = SSsub[+];
      SSerror = SSw-SSsub;                   /* SSerror */
      MSt = SSb/(k-1);                       /* MSt */     
      MSerror = SSerror/((n-1)*(k-1));       /* MSerror */
      F = MSt/MSerror;                       /* F statistic */
      return(F);
   finish;

/* module to independently permute elements of each row of a matrix */
/* see SAS Communities: https://communities.sas.com/message/212330 */
/* see "The DO Loop": http://blogs.sas.com/content/iml/2014/05/29/permute-elements-within-each-row-of-a-matrix.html */

   start PermuteWithinRows(m);
      colIdx = ranperm(1:ncol(m), nrow(m));   /* permute col. indices */ 
      f = (row(m) - 1) * ncol(m);             /* constant for row */ 
      matIdx = f + colIdx;                    /* matrix indices */
      return( shape(m[matIdx], nrow(m)) );
   finish;

   fobt = fmod(xobt);                        /* calculate F for observed data */
   print fobt;
   call symputx('fobt',fobt);           

   call randseed(12345);
   B = 10000;                                /* number or reps */
   fdist = j(B,1);                           /* allocate vector to hold bootstrap dist of F */
   do j = 1 to B;                            /* bootstrap loop */
      x = PermuteWithinRows(xobt);           /* call PermuteWithinRows module */   
      F = fmod(x);                           /* calculate F on permuted data */
      fdist[j,] = F;                         /* store calc’d F in fdist vector */
   end;

/* compute p-value */
   pval = sum(fdist > abs(fobt)) / B;
   print pval[label='P-value'];
   call symputx('p',pval);

/* create dataset */
   create bootf var {fdist};
   append;
   close bootf;
quit;
