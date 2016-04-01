	//-----------------------------------------------------------
	// sonic/supersonic outflow boundary condition implementation
	// 
	// In 0/p file, sideRight boundary condition must be of 
	// fixed value type 
	//-----------------------------------------------------------		
	label patchID = mesh.boundaryMesh().findPatchID("sideRight");
	
	//const scalarField& p_test = p.boundaryField()[patchID];
	    
	// make a reference to cell adjacent to boundary patch
	scalar T1_interior = thermo1.T()[celln];
	scalar T2_interior = thermo2.T()[celln];
	scalar rho1_interior = rho1[celln];
	scalar rho2_interior = rho2[celln];
	scalar alpha1_interior = alpha1[celln];
	scalar alpha2_interior = alpha2[celln]; 
  	scalar p_interior = p[(celln)];
	scalar U_interior = mag(U1[(celln)])*alpha1[celln]+mag(U2[celln])*alpha2[celln];
 	scalar rho_interior = alpha1[(celln)]*rho1[(celln)]+alpha2[(celln)]*rho2[(celln)];
	
	Info << "pressure1 = " << p_interior << nl << endl;	
	Info << "velocity1 = " << U_interior << nl << endl;
	Info << "density1 = " << rho_interior << nl << endl;	
 	Info << "temperatue1 = " << thermo1.T()[celln] << nl << endl;
	//Info << "enthalpy1 = " << thermo1.he()[celln] << nl << endl; 

  	// declare the varaibles of ghost cell for supersonic and subsonic outflow
	// ghost cell 1
  	scalar p_ghost_update, U_ghost_update, rho_ghost_update, T_ghost_update, rho1_ghost_update, rho2_ghost_update, alpha1_ghost_update;
  	// ghost cell 2
	scalar p_ghost_update2, U_ghost_update2, rho_ghost_update2;
	// ghost cell 1
	scalar U_ghost, p_ghost, rho_ghost, rho1_ghost, rho2_ghost, T1_ghost, T2_ghost, alpha1_ghost, alpha2_ghost;
	// ghost cell 2
	scalar U_ghost2, p_ghost2, rho_ghost2;
	// speed of sound 
	scalar sos, sos1, sos2;

  	// obtaining the ghost cell information

	if (runTime.time().value() == 0.0)
	{	
	    U_ghost = scalar(0);
	    p_ghost = scalar(599000);
	    rho_ghost = scalar(33.295);
      	    //rho1_ghost = scalar(33.295);
      	    //rho2_ghost = scalar(1088.9);
      	    //T1_ghost = scalar(220.00);
      	    //T2_ghost = scalar(220.00);
            //alpha1_ghost = scalar(0.99)*(rho_ghost/rho1_ghost);
            //alpha2_ghost = scalar(1.0) - alpha1_ghost;
	}
  	else
  	{
  	    U_ghost = U_ghost_update;  	    //mag(U.boundaryField()[patchID][0]);
      	    p_ghost = scalar(599000);       //p.boundaryField()[patchID][0];
	    rho_ghost = rho_ghost_update;   //alpha1.boundaryField()[patchID][0]*rho1.boundaryField()[patchID][0]+alpha2.boundaryField()[patchID][0]*rho2.boundaryField()[patchID][0];
	    //rho1_ghost = rho1.boundaryField()[patchID][0];
      	    //rho2_ghost = rho2.boundaryField()[patchID][0];
      	    //T1_ghost = thermo1.T().boundaryField()[patchID][0]; 
      	    //T2_ghost = thermo2.T().boundaryField()[patchID][0];
      	    //alpha1_ghost = alpha1.boundaryField()[patchID][0]; 
      	    //alpha2_ghost = alpha2.boundaryField()[patchID][0];
  	}
	
    	  Info << "pressure2 = " << p.boundaryField()[patchID][0] << nl << endl;
	  Info << "velocity2 = " << U1.boundaryField()[patchID][0] << nl << endl;
	  Info << "velocity2 = " << U2.boundaryField()[patchID][0] << nl << endl;
	  Info << "densityVap2 = " << rho1.boundaryField()[patchID][0] << nl << endl;
	  Info << "densityLiq2 = " << rho2.boundaryField()[patchID][0] << nl << endl;	
          Info << "temperature2 = " << thermo2.T().boundaryField()[patchID][0] << nl << endl;
    //
    // obtaining the speed of sound      
    //
	 
    long ierr = 0; 
    char herr[255];
    double x_cal[] = {1.0}, xliq_cal[] = {0.0}, xvap_cal[] = {0.0}, q_cal = NULL,T_cal = NULL, rho_cal = NULL, rholiq_cal=NULL, rhovap_cal=NULL, p_cal=NULL, h_cal = NULL, s_cal = NULL, e_cal = NULL, cp_cal = NULL, cv_cal = NULL, w_cal = NULL;
    double hjt = NULL;	
    double MW = 44.01; // kg/kmol

    T_cal = T1_interior, rho_cal = rho1_interior/MW;
    THERMdll(T_cal,rho_cal,x_cal,p_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,hjt);
	   	
    if (w_cal = w_cal) sos1 = w_cal;
    else if (w_cal != w_cal) sos1 = sos1;	

    T_cal = T2_interior, rho_cal = rho2_interior/MW;
    THERMdll(T_cal,rho_cal,x_cal,p_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,hjt);

    if (w_cal = w_cal) sos2 = w_cal;
    else if (w_cal != w_cal) sos2 = sos2;	

    // bulk speed of sound (weighted against volume fractions)
    if (alpha1_interior > 0.0 && alpha1_interior < 1.0)
    {
    	sos = Foam::sqrt(1.0/(rho_interior*(alpha1_interior/(rho1_interior*sos1*sos1)+alpha2_interior/(rho2_interior*sos2*sos2)))); 
    }
    else if (alpha1_ghost == 1.0) sos = sos1;
    else if (alpha1_ghost == 0.0) sos = sos2;
  	
    Info << "sound speed =" << sos1 << nl << endl;	
    Info << "sound speed =" << sos2 << nl << endl;	
    Info << "sound speed =" << sos << nl << endl;  


    // following the Mach line and acoustic wave  	
    rho_ghost_update = rho_interior + (p_ghost-p_interior)/(sos*sos);	
		
    U_ghost_update = U_interior - (p_ghost-p_interior)/(rho_interior*sos); 


    // Thermodynamic calculations
		
  	// singlephase or twophase PD-flash calculation
 	p_cal = p_ghost/1000.0; // kpa
  	rho_cal = rho_ghost_update/MW; // mol/L
  	PDFLSHdll(p_cal,rho_cal,x_cal,T_cal,rholiq_cal,rhovap_cal,xliq_cal,xvap_cal,q_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,ierr,herr,255);
   
	T_ghost_update = T_cal;
  
  	if (q_cal >= 0.0 && q_cal <= 1.0) 
  	{
    		alpha1_ghost_update = q_cal*(rho_cal/rhovap_cal);
  	}
  	else if (q_cal > 1.0)
  	{
    		alpha1_ghost_update = 1.0;
  	}
  	else 
  	{
    		alpha1_ghost_update = 0.0;
  	}

  	Info << "p_ghost = " << p_ghost << nl << endl;	
	Info << "U_ghost_update = " << U_ghost_update << nl << endl;
	Info << "rho_ghost_update = " << rho_ghost_update << nl << endl;
	Info << "T_ghost_update = " << T_ghost_update << nl << endl;
  	Info << "alpha1_ghost_update = " << alpha1_ghost_update << nl << endl;
  
  //-----------------------------------------------------------

  //-----------------------------------------------------------
  // subsonic outflow boundary condition
  //-----------------------------------------------------------
  
  // obtaining current timestep deltaT
  scalar dt = runTime.deltaT().value();

  // declare essential variables	
  scalar area1 = 233.0*233.0, area2 = 50.0*50.0, dischargeCoeff = 0.9, areaRatio = dischargeCoeff*(area2/area1);
  scalar rhoU_ghost, rhoU_ghost_update, ddrhoU_ghost;	

  //ghost cell information for subsonic outflow
  if (runTime.time().value() == 0.0)
  {	
	    U_ghost2 = scalar(0);
	    p_ghost2 = scalar(599000);
	    rho_ghost2 = scalar(33.295);
      	    rho1_ghost = scalar(33.295);
      	    rho2_ghost = scalar(1088.9);
      	    T1_ghost = scalar(220.00);
      	    T2_ghost = scalar(220.00);
            alpha1_ghost = scalar(1.0)*(rho_ghost/rho1_ghost);
            alpha2_ghost = scalar(1.0) - alpha1_ghost;
   }
   else
   {
  	    U_ghost2 = U_ghost_update2;
      	    p_ghost2 = p.boundaryField()[patchID][0];
	    rho_ghost2 = alpha1.boundaryField()[patchID][0]*rho1.boundaryField()[patchID][0]+alpha2.boundaryField()[patchID][0]*rho2.boundaryField()[patchID][0];
	    rho1_ghost = rho1.boundaryField()[patchID][0];
      	    rho2_ghost = rho2.boundaryField()[patchID][0];
      	    T1_ghost = thermo1.T().boundaryField()[patchID][0]; 
      	    T2_ghost = thermo2.T().boundaryField()[patchID][0];
      	    alpha1_ghost = alpha1.boundaryField()[patchID][0]; 
      	    alpha2_ghost = alpha2.boundaryField()[patchID][0];
   }	
 
  //calculat rho*U for the new ghost cell (subsonic ghost cell)  
  rhoU_ghost = rho_ghost*U_ghost*dischargeCoeff*areaRatio;	
  rhoU_ghost_update = rho_ghost_update*U_ghost_update*dischargeCoeff*areaRatio;	 	
  ddrhoU_ghost = (rhoU_ghost_update-rhoU_ghost)/dt; 	
	
  //wave amplitude calculation  
  scalar wave22 = (U_ghost2/xDimDim[celln])*(p_ghost2 - p_interior-sos*sos*(rho_ghost2-rho_interior));
  scalar wave32 = ((U_ghost2+sos)/xDimDim[celln])*(p_ghost2-p_interior+rho_ghost2*sos*(U_ghost2-U_interior));
  scalar wave12 = (1.0/(U_ghost2-sos))*(-2.0*U_ghost2*wave22-2.0*ddrhoU_ghost*sos*sos-(U_ghost2+sos)*wave32); 	

  //solving characteristic form governing equations	
  scalar DpDt2 = -0.5*(wave32+wave12);
  scalar DuDt2 = -0.5/(rho_ghost2*sos)*(wave32-wave12);
  scalar DrhoDt2 = 1.0/(sos*sos)*(DpDt2+wave22);

  p_ghost_update2 = p_ghost2 + DpDt2*dt;
  U_ghost_update2 = U_ghost2 + DuDt2*dt;
  rho_ghost_update2 = rho_ghost2 + DrhoDt2*dt;			

  Info << "drhoUdt = " << ddrhoU_ghost << nl << endl;	
  Info << "p' = " << p_ghost_update2 << nl << endl;
  Info << "U' = " << U_ghost_update2 << nl << endl;	  
  Info << "rho' = " << rho_ghost_update2 << nl << endl;	
