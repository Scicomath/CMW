double anisospace()
{
  double diff=0,sum=0;
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	diff+=e[sx][sy]*((sy-Middle)*(sy-Middle)-(sx-Middle)*(sx-Middle));
	sum+=e[sx][sy]*((sx-Middle)*(sx-Middle)+(sy-Middle)*(sy-Middle));
      }

  return diff/sum;
}



double anisomomentum()
{
  double diff=0,sum=0;
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	diff+=(e[sx][sy]+eos(e[sx][sy]))*(u[0][sx][sy]*u[0][sx][sy]-u[1][sx][sy]*u[1][sx][sy]);
	diff+=pixx[sx][sy]-piyy[sx][sy];
	sum+=(e[sx][sy]+eos(e[sx][sy]))*(u[0][sx][sy]*u[0][sx][sy]+u[1][sx][sy]*u[1][sx][sy]);
	sum+=pixx[sx][sy]+piyy[sx][sy];
	sum+=2*eos(e[sx][sy]);
      }
  return diff/sum;
}






double uphi(int sx,int sy)
{
  double temp=(u[1][sx][sy]*sx-u[0][sx][sy]*sy)/sqrt(sx*sx+sy*sy);
  return temp;
}

double ur(int sx,int sy)
{
  double temp=(u[0][sx][sy]*sx+u[1][sx][sy]*sy)/sqrt(sx*sx+sy*sy);
  return temp;
}

int foundit(int sx, int sy)
{
  globali=geti(e[sx][sy]);
  globalx=getx(globali,e[sx][sy]);

  freeze_out << (sx-Middle)/5.06842*AT << "\t";
  freeze_out << (sy-Middle)/5.06842*AT << "\t";
  freeze_out << u[0][sx][sy] << "\t";
  freeze_out << u[1][sx][sy] << "\t";
  freeze_out << pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
  freeze_out << pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
  freeze_out << piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\n";
  return 0;
}

void fancyfoundx(double sx, int sy,streambuf* pbuf)
{

  //fstream str;

  streambuf* important;
  
  important=cout.rdbuf();


  cout.rdbuf(pbuf);

  cout << (sx-Middle)/5.06842*AT << "\t";
  cout << (sy-Middle)/5.06842*AT << "\t";

  double workhorse[NUMT];
  double iarr[NUMT];
  gsl_spline * hspl=gsl_spline_alloc (gsl_interp_cspline, NUMT);
  gsl_interp_accel * hsacc=gsl_interp_accel_alloc ();


  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=u[0][i][sy];
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=u[1][i][sy];
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=pixx[i][sy]/(e[i][sy]+eos(e[i][sy]));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=pixy[i][sy]/(e[i][sy]+eos(e[i][sy]));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=piyy[i][sy]/(e[i][sy]+eos(e[i][sy]));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\n";

  /*
  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=T(i,sy)/AT;
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\n";
  */

  gsl_spline_free (hspl);
  gsl_interp_accel_free (hsacc);

  cout.rdbuf(important);

}

void fancyfoundy(int sx, double sy,streambuf* pbuf)
{

  streambuf* important;
  
  important=cout.rdbuf();


  //fstream str;

  cout.rdbuf(pbuf);

  cout << (sx-Middle)/5.06842*AT << "\t";
  cout << (sy-Middle)/5.06842*AT << "\t";

  double workhorse[NUMT];
  double iarr[NUMT];
  gsl_spline * hspl=gsl_spline_alloc (gsl_interp_cspline, NUMT);
  gsl_interp_accel * hsacc=gsl_interp_accel_alloc ();
  

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=u[0][sx][i];
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=u[1][sx][i];
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=pixx[sx][i]/(e[sx][i]+eos(e[sx][i]));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=pixy[sx][i]/(e[sx][i]+eos(e[sx][i]));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=piyy[sx][i]/(e[sx][i]+eos(e[sx][i]));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\n";

  gsl_spline_free (hspl);
  gsl_interp_accel_free (hsacc);

  cout.rdbuf(important);
}


double rootfunction(double x,void *params)
{
  double temp;
  temp=gsl_spline_eval (workspline, x, wac);
  return (temp-TF);
}

//fancy freeze-out with interpolation
void fancyfreeze()
{
  int dflag=1;
  int debug=0;

  /* too risky
  //check whether there is non-trivial freeze-out
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	globali=geti(e[sx][sy]);
	globalx=getx(globali,e[sx][sy]);
	if (T(sx,sy)>TF)
	  dflag=0;
      }
  */
  if (T(Middle,Middle)>TF)
    dflag=0;

  //if yes, use whole lattice to determine 
  //freeze-out surface as accurate as possible
  //(important for corse lattices)
  if (dflag==0)
    {
      
      if (debug==1)
	printf("Starting freeze-out\n");

      double warr[Middle];
      double iarr[Middle];

      int limiter;
      int status;
      int iter = 0, max_iter = 1000;
      const gsl_root_fsolver_type *TT;
      gsl_root_fsolver *ss;
      double r=0;

      gsl_function F;
      double dummy=0;

      F.function=&rootfunction;
      F.params=&dummy;

      TT = gsl_root_fsolver_brent;
      ss = gsl_root_fsolver_alloc (TT);

      //endpoint
      int ssy=Middle;
      int ssx=Middle;
      for (ssx=Middle;ssx<=NUMT;ssx++)
	{
	  globali=geti(e[ssx][ssy]);
	  globalx=getx(globali,e[ssx][ssy]);
	  warr[ssx-Middle]=T(ssx,ssy);
	  iarr[ssx-Middle]=ssx;
	  if (debug==1)
	    printf("%i %f\n",ssx,warr[ssx-Middle]);
	}

      gsl_spline_init (workspline,iarr,warr,Middle);

      gsl_root_fsolver_set (ss, &F, Middle, NUMT);
      iter=0;
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate (ss);
	  r = gsl_root_fsolver_root (ss);
	  double x_lo = gsl_root_fsolver_x_lower (ss);
          double x_hi = gsl_root_fsolver_x_upper (ss);
	  status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
      fancyfoundx(r,ssy,freeze_out.rdbuf());

      //correct slightly to
      //inhibit bad conversion properties
      limiter=(int) (r-(AT/20));

      printf("freezel=%f\n",(limiter-Middle)/5.06842*AT);

      //if we are very close to the finish
      //decrease step size to capture
      //important final effects
      //note:does not work so great
      //would need more complicated adjustment
      //if (((limiter-Middle)/5.06842*AT)<1.5)
      //	{
      //	  //SNAPUPDATE=(int) (SNAPUPDATE/2.);
      //	  //SNAPUPDATE++;
      //	  UPDATE=(int) (UPDATE/1.2);
      //	  UPDATE++;
      //	}

      if (debug==1)
	printf("E1 done limit=%i, %f\n",limiter,r);



      //by scanning in y:
      //upper right quarter
      for (int sx=limiter;sx>=Middle;sx--)
	{
	  for (int sy=Middle;sy<=NUMT;sy++)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      warr[sy-Middle]=T(sx,sy);
	      iarr[sy-Middle]=sy;
	    }
	  
	  gsl_spline_init (workspline,iarr,warr,Middle);

	  iter = 0;
	  gsl_root_fsolver_set (ss, &F, Middle, NUMT);
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (ss);
	      r = gsl_root_fsolver_root (ss);
	      double x_lo = gsl_root_fsolver_x_lower (ss);
	      double x_hi = gsl_root_fsolver_x_upper (ss);
	      status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  fancyfoundy(sx,r,freeze_out.rdbuf());

	  if (debug==1)
	    printf("fq: %i done\n",sx);
	}



      //upper left quarter
      for (int sx=Middle-1;sx>=2*Middle-limiter;sx--)
	{
	  for (int sy=Middle;sy<=NUMT;sy++)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      warr[sy-Middle]=T(sx,sy);
	      iarr[sy-Middle]=sy;
	    }
	  
	  gsl_spline_init (workspline,iarr,warr,Middle);
	  
	  gsl_root_fsolver_set (ss, &F, Middle, NUMT);
	  iter=0;
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (ss);
	      r = gsl_root_fsolver_root (ss);
	      double x_lo = gsl_root_fsolver_x_lower (ss);
	      double x_hi = gsl_root_fsolver_x_upper (ss);
	      status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  fancyfoundy(sx,r,freeze_out.rdbuf());
	  if (debug==1)
	    printf("sq: %i done\n",sx);
	}

      
      //endpoint
      ssy=Middle;
      ssx=Middle;
      for (ssx=1;ssx<=Middle;ssx++)
	{
	  globali=geti(e[ssx][ssy]);
	  globalx=getx(globali,e[ssx][ssy]);
	  warr[ssx-1]=T(ssx,ssy);
	  iarr[ssx-1]=ssx;
	}

      gsl_spline_init (workspline,iarr,warr,Middle);

      gsl_root_fsolver_set (ss, &F, 1,Middle);
      iter=0;
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate (ss);
	  r = gsl_root_fsolver_root (ss);
	  double x_lo = gsl_root_fsolver_x_lower (ss);
          double x_hi = gsl_root_fsolver_x_upper (ss);
	  status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	}
      while (status == GSL_CONTINUE && iter < max_iter);      
      fancyfoundx(r,ssy,freeze_out.rdbuf());

      if (debug==1)
	printf("E2 done %i, %f\n",2*Middle-limiter,r);


      //lower left quarter
      for (int sx=2*Middle-limiter;sx<=Middle;sx++)
	{
	  for (int sy=1;sy<=Middle;sy++)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      warr[sy-1]=T(sx,sy);
	      iarr[sy-1]=sy;
	    }

	  gsl_spline_init (workspline,iarr,warr,Middle);
	  
	  gsl_root_fsolver_set (ss, &F, 1,Middle);
	  iter=0;
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (ss);
	      r = gsl_root_fsolver_root (ss);
	      double x_lo = gsl_root_fsolver_x_lower (ss);
	      double x_hi = gsl_root_fsolver_x_upper (ss);
	      status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  fancyfoundy(sx,r,freeze_out.rdbuf());
	  if (debug==1)
	    printf("tq: %i done\n",sx);
	}

      //lower right quarter
      for (int sx=Middle+1;sx<=limiter;sx++)
	{
	  for (int sy=1;sy<=Middle;sy++)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      warr[sy-1]=T(sx,sy);
	      iarr[sy-1]=sy;
	    }

	  gsl_spline_init (workspline,iarr,warr,Middle);
	  
	  gsl_root_fsolver_set (ss, &F, 1,Middle);
	  iter=0;
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (ss);
	      r = gsl_root_fsolver_root (ss);
	      double x_lo = gsl_root_fsolver_x_lower (ss);
	      double x_hi = gsl_root_fsolver_x_upper (ss);
	      status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  fancyfoundy(sx,r,freeze_out.rdbuf());
	  if (debug==1)
	    printf("4q: %i done\n",sx);
	}

      
      //cout << " Found " << r << endl;


      
      gsl_root_fsolver_free (ss);
      
      freeze_out << "TIME \t" << t/5.06842*AT << endl;
    }
  else
    reachedTf=1;
}

void stupidfreeze()
{

  if (T(Middle,Middle)<TF)
    {
      reachedTf=1;
      
      for (int sx=1;sx<=NUMT;sx++)
	for (int sy=1;sy<=NUMT;sy++)
	  {
	    freeze_out << (sx-Middle)/5.06842*AT << "\t";
	    freeze_out << (sy-Middle)/5.06842*AT << "\t";
	    freeze_out << u[0][sx][sy] << "\t";
	    freeze_out << u[1][sx][sy] << "\t";
	    freeze_out << pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	    freeze_out << pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	    freeze_out << piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	    freeze_out << T(sx,sy)/AT << "\n";
	  }
      freeze_out << "TIME \t" << t/5.06842*AT << endl;
    }
}

//constant temperature freeze-out
void freezeout()
{
  int dflag=1;

  //endpoint
  int ssy=Middle;
  int ssx=Middle;
  globali=geti(e[ssx][ssy]);
  globalx=getx(globali,e[ssx][ssy]);
  while(T(ssx,ssy)>TF)
    {
      ssx++;
      globali=geti(e[ssx][ssy]);
      globalx=getx(globali,e[ssx][ssy]);
    }
  if (ssx!=Middle)
    dflag=foundit(ssx-1,ssy);


  //by scanning in y:
  //upper right quarter
  for (int sx=NUMT;sx>=Middle;sx--)
    {
      int sy=Middle;
      globali=geti(e[sx][sy]);
      globalx=getx(globali,e[sx][sy]);
      while(T(sx,sy)>TF)
	{
	  sy++;
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	}
      if (sy!=Middle)
	dflag=foundit(sx,sy-1);
    }

  //upper left quarter
  for (int sx=Middle-1;sx>=1;sx--)
    {
      int sy=Middle;
      globali=geti(e[sx][sy]);
      globalx=getx(globali,e[sx][sy]);
      while(T(sx,sy)>TF)
	{
	  sy++;
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	}
      if (sy!=Middle)
	dflag=foundit(sx,sy-1);
    }

    //endpoint
  ssy=Middle;
  ssx=Middle;
  globali=geti(e[ssx][ssy]);
  globalx=getx(globali,e[ssx][ssy]);
  while(T(ssx,ssy)>TF)
    {
      ssx--;
      globali=geti(e[ssx][ssy]);
      globalx=getx(globali,e[ssx][ssy]);
    }
  if (ssx!=Middle)
    dflag=foundit(ssx+1,ssy);

   //lower left quarter
  for (int sx=1;sx<=Middle;sx++)
    {
      int sy=Middle;
      globali=geti(e[sx][sy]);
      globalx=getx(globali,e[sx][sy]);
      while(T(sx,sy)>TF)
	{
	  sy--;
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	}
      if (sy!=Middle)
	dflag=foundit(sx,sy+1);
    }

  //lower right quarter
  for (int sx=Middle+1;sx<=NUMT;sx++)
    {
      int sy=Middle;
      globali=geti(e[sx][sy]);
      globalx=getx(globali,e[sx][sy]);
      while(T(sx,sy)>TF)
	{
	  sy--;
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	}
      if (sy!=Middle)
	dflag=foundit(sx,sy+1);
    }

  
  //if we could find any fluid element with
  //T>TF, then we're done
  if (dflag==1)
    reachedTf=1;
  else
    freeze_out << "TIME \t" << t/5.06842*AT << endl;

}

//freeze-out that doesn't assume a monotonically 
//decreasing temperature from the center out
void blockfreeze()
{
//hydro evolution is finished if everywhere in the time slice T < TF
  if (T(Middle,Middle) < TF) reachedTf = 1;

  //All grid points that are nearest neighbors in x, y, or tau are checked pairwise.
  //A rectangular piece of the freezeout surface is defined to be half-way between
  //any pair whose temperatures straddle TF.
  for (int sx=1;sx<=NUMT;sx++) {
    for (int sy=1;sy<=NUMT;sy++)
    {
      if (sy != NUMT)
      {
	if (((T(sx,sy) > TF) && (T(sx,sy+1) <= TF)) || ((T(sx,sy) <= TF) && (T(sx,sy+1) > TF)))
	{ 
	  reachedTf = 0;
	  int direction = 2;
	  if ((T(sx,sy) <= TF) && (T(sx,sy+1) > TF)) direction = -2;
	  freeze_out << (sx-Middle)/5.06842*AT << "\t";
	  freeze_out << (sy-Middle+0.5)/5.06842*AT << "\t";
	  freeze_out << t/5.06842*AT << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + u[0][sx][sy+1]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + u[1][sx][sy+1]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixx[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1]))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixy[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1])) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
			      + piyy[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1])))<< "\t";
	  freeze_out << 0.5 * (T(sx,sy) + T(sx,sy+1))/AT << "\n";
	}
      }
      if (sx != NUMT)
      {
      	if (((T(sx,sy) > TF) && (T(sx+1,sy) <= TF)) || ((T(sx,sy) <= TF) && (T(sx+1,sy) > TF)))
	{ 
	  reachedTf = 0;
	  int direction = 1;
	  if ((T(sx,sy) <= TF) && (T(sx+1,sy) > TF)) direction = -1;
	  freeze_out << (sx-Middle+0.5)/5.06842*AT << "\t";
	  freeze_out << (sy-Middle)/5.06842*AT << "\t";
	  freeze_out << t/5.06842*AT << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + u[0][sx+1][sy]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + u[1][sx+1][sy]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixx[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy]))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixy[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy])) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
			      + piyy[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy])))<< "\t";
	  freeze_out << 0.5 * (T(sx,sy) + T(sx+1,sy))/AT << "\n";
	}
      }

      if (t != TINIT*5.06842/AT)
      {
	if ((T(sx,sy) <= TF && (Tlast(sx,sy) > TF)) || ((T(sx,sy) > TF) && (Tlast(sx,sy) <= TF)))
	{
	  int direction = 3;
	  if ((T(sx,sy) > TF) && (Tlast(sx,sy) <= TF)) direction = -3;
	  freeze_out << (sx-Middle)/5.06842*AT << "\t";
	  freeze_out << (sy-Middle)/5.06842*AT << "\t";
	  freeze_out << t/5.06842*AT - 0.5 * UPDATE*EPS*AT/5.06842 << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + ulast[0][sx][sy]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + ulast[1][sx][sy]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 	
			      + pixxlast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy]))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) 
			      + pixylast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy])) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy]))
			      + piyylast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy])))<< "\t";
	  freeze_out << 0.5 * (T(sx,sy) + Tlast(sx,sy))/AT << "\n";
	}
      }
    }
  }
}

void outputMeasurements(double t) 
{
  cout.precision(5);
  int dwidth = 13;
  cout.width(dwidth); cout << t/5.06842*AT;
  globali=geti(e[Middle][Middle]);
  globalx=getx(globali,e[Middle][Middle]);
  double ex=anisospace();
  double ep=anisomomentum();
  cout.width(dwidth); cout << T(Middle,Middle)/AT;
  cout << "\t" << ex << "\t" << ep;
  //cout << "\t" << overlapS()/AT/AT;
  cout << endl;
  
  ecces << t/5.06842*AT <<"\t";
  ecces << ex << "\t";
  ecces << ep << "\n";

    //freezeout();
//   if (FREEZE==1)
//     fancyfreeze();
//   else
//     stupidfreeze();
    
  switch (FREEZE){
    case 0:
      stupidfreeze();
      break;
    case 1:
      fancyfreeze();
      break;
    case 2:
      blockfreeze();
      break;
    default:
      blockfreeze();
  }
}

void snapTprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Tprofile_%.2f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/5.06842*AT << "\t";
    out << T(Middle,s)/AT << endl;
  }
  out.close();
}

void snapEDprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/EDprofile_%.2f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    
    out << (s-Middle)/5.06842*AT << "\t";
    out << e[Middle][s]/AT/AT/AT/AT << endl;
  }
  out.close();
}

void snapTcontour(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Tcontour_%.3f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	  out << T(sx,sy)/AT << "\t";
	}
      out << endl;
    }
  out.close();
}

void snapFOdata(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/FOdata_%.3f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  out << u[0][sx][sy] << "\t";
	  out << u[1][sx][sy] << "\t";
	  out << pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	  out << pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	  out << piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy])) << "\t";
	}
      out << endl;
    }
  out.close();
  //putting names into 
  //file to facilitate later 
  //freeze-out calc
  meta << fname << "\n";
}

void snapVprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Vprofile_%.2f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/5.06842*AT << "\t";
    out << u[1][Middle][s]/ut(Middle,s) << endl;
  }
  out.close();
}

void snapVxprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Vxprofile_%.2f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[s][Middle]);
    globalx=getx(globali,e[s][Middle]);
    out << (s-Middle)/5.06842*AT << "\t";
    out << u[0][s][Middle]/ut(s,Middle) << endl;
  }
  out.close();
}

void snapV2profile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/V2profile_%.2f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int s=0;s<NUMT/2./sqrt(2.);s++)
  {
    globali=geti(e[Middle+s][Middle+s]);
    globalx=getx(globali,e[Middle+s][Middle+s]);
    out << s*sqrt(2)/5.06842*AT << "\t";
    out << (u[1][Middle+s][Middle+s]+u[0][Middle+s][Middle+s])/sqrt(2.)/ut(Middle+s,Middle+s) << endl;
  }
  out.close();
}

void snappieeprofile(double time)
{
  fstream out;
  char fname[255];
  double vr=0;
  sprintf(fname,"data/snapshot/Piprofile_%.2f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/5.06842*AT << "\t";
    out << -pi(3,3,Middle,s)*t*t/(e[Middle][s]*4./3.)<< endl;
  }
  out.close();
}

void snappirrprofile(double time)
{
  fstream out;
  char fname[255];
  double vr=0;
  sprintf(fname,"data/snapshot/PiRprofile_%.2f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/5.06842*AT << "\t";
    out << -pi(1,1,Middle,s)/(4./3.*e[Middle][s])<< endl;
  }
  out.close();
}


void snapuphiprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/uphiprofile_%.2f.dat",time/5.06842*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/5.06842*AT << "\t";
    out << u[0][Middle][s]/ut(Middle,s) << endl;
  }
  out.close();
}


void cheat(int sx,int sy)
{
  globali=geti(e[sx][sy]);
  globalx=getx(globali,e[sx][sy]);

  globut=ut(sx,sy);

  thf[0]=theta0(sx,sy);
  thf[1]=theta1(sx,sy);
  thf[3]=theta3(sx,sy);
  
  mypixt=pi(0,2,sx,sy);
  mypiyt=pi(1,2,sx,sy);
  mypitt=pi(2,2,sx,sy);
  mypiee=pi(3,3,sx,sy);
  
  dtpiialpha(dtpixx,0,0,sx,sy);
  dtpiialpha(dtpixy,0,1,sx,sy);
  dtpiialpha(dtpiyy,1,1,sx,sy);
}


void fordebug(int sx, int sy)
{
  //systematically:

  cout << "-----------------------------------------------\n";
  printf("thing \t 25,22 \t 65,68\n");
  //cout << "thing\t" << "25,22" <<"\t" << "65,68" << endl;
  

  printf("e\t %.12g \t %.12g\n",e[25][22],e[65][68]);
  //printf("this is a test %f9\n",1.12345678910);
  printf("u[0]\t %.12g \t %.12g\n",u[0][25][22],u[0][65][68]);
  printf("u[1]\t %.12g \t %.12g\n",u[1][25][22],u[1][65][68]);
  //cout << "u[0]\t" << u[0][25][22] << "\t" <<  u[0][65][68] << endl;
  //cout << "u[1]\t" << u[1][25][22] << "\t" <<  u[1][65][68] << endl;

  cheat(25,22);
  printf("pixx\t %.12g \t",pixx[25][22]);
  cheat(65,68);
  printf("%.12g\n",pixx[65][68]);

  cheat(25,22);
  printf("pixy\t %.12g \t",pixy[25][22]);
  cheat(65,68);
  printf("%.12g\n",pixy[65][68]);

  cheat(25,22);
  printf("piyy\t %.12g \t",piyy[25][22]);
  cheat(65,68);
  printf("%.12g\n",piyy[65][68]);

  cheat(25,22);
  printf("pixt\t %.12g \t",pishell(0,2,25,22));
  cheat(65,68);
  printf("%.12g\n",pishell(0,2,65,68));

  cheat(25,22);
  printf("piyt\t %.12g \t",pishell(1,2,25,22));
  cheat(65,68);
  printf("%.12g\n",pishell(1,2,65,68));

  cheat(25,22);
  printf("pitt\t %.12g \t",pishell(2,2,25,22));
  cheat(65,68);
  printf("%.12g\n",pishell(2,2,65,68));
}

void snapshot(double tt)
{
  double hel2[4];
  double dt[4];
  double tb[4],td[4];
 
  //fordebug(10,10);
  
  snapTcontour(tt);
  snapFOdata(tt);

  snapTprofile(tt);   
  snapEDprofile(tt);
  snapVprofile(tt);
  //snapVxprofile(tt);
  //snapV2profile(tt);
  //snapV2profile(tt);
  //snapuphiprofile(tt);
  snappieeprofile(tt);
  snappirrprofile(tt);
}
