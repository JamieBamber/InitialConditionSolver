/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYPOTENTIALFUNCTION_HPP
#define MYPOTENTIALFUNCTION_HPP

void my_potential_function(Real &V_of_phi, Real &dVdphi, const Real &phi_here, Real Lambda, Real mu, string inflation_model)
{
	double mpl=1.0;

	if (inflation_model == "phi43"){

        // phi^4/3
        V_of_phi = Lambda * pow(-phi_here, 4.0/3.0);

	}
	else if (inflation_model == "phi23"){

	// phi^2/3
        double phi_flat = 2e-4;
        double phi_reheat = -1e-4;
        double phi_reheat_min = phi_flat;
        double Lambda_reheat = 23938.2* Lambda;

        Real V_of_phi_inf = Lambda * pow(-phi_here, 2.0/3.0);
        Real dVdphi_inf = -2.0/3.0 * Lambda * pow(-phi_here, -1./3.0);
        
        Real V_of_phi_reheat = Lambda_reheat*pow((phi_here - phi_reheat_min) / mpl, 2.0);
        Real dVdphi_reheat = 2.0 * Lambda_reheat*(phi_here - phi_reheat_min) * pow(mpl, -2.0);

        Real V_of_phi_flat = 0.0;
        Real dVdphi_flat = 0.0;

		if (phi_here < phi_reheat){
                        V_of_phi = V_of_phi_inf;
                        dVdphi = dVdphi_inf;}
		if (phi_here > phi_reheat && phi_here < phi_reheat_min){
                        V_of_phi = V_of_phi_reheat;
                        dVdphi = dVdphi_reheat;}
		if (phi_here > phi_reheat_min){
                        V_of_phi = V_of_phi_flat;
                        dVdphi = dVdphi_flat;}
	}	
	else if (inflation_model == "starobinsky"){

	// Starobinsky
		V_of_phi = Lambda*pow(1. - exp((phi_here/mu*pow(mpl,-1))),2);
                dVdphi = -2.*Lambda/mu*(1 - exp((phi_here/mu*pow(mpl,-1))))*exp((phi_here/mu*pow(mpl,-1)))*pow(mpl,-1);

	}
	else if (inflation_model == "quadratic"){

	// Hilltop Quadratic
	double phi_flat = 0.0;
        double phi_reheat = 0.9*mu;
        double phi_reheat_min = 1.11111*mu;
        double Lambda_reheat = 4.26316*Lambda*pow(mu,-2.0);

        Real V_of_phi_flat = Lambda;
        Real dVdphi_flat = 0.0;

        Real V_of_phi_inf = Lambda*(1 - pow(mu,-2.)*pow(phi_here,2.));
        Real dVdphi_inf = -2.0*Lambda*pow(mu,-2.)*pow(phi_here,1.);

        Real V_of_phi_reheat = Lambda_reheat*pow((phi_here - phi_reheat_min) / mpl, 2.0);
        Real dVdphi_reheat = 2.0 * Lambda_reheat*(phi_here - phi_reheat_min) * pow(mpl, -2.0);

		if (phi_here < phi_flat){
                        V_of_phi = V_of_phi_flat;
                        dVdphi = dVdphi_flat;}
                if (phi_here > phi_flat && phi_here < phi_reheat){
                        V_of_phi = V_of_phi_inf;
                        dVdphi = dVdphi_inf;}
                if (phi_here > phi_reheat){
                        V_of_phi = V_of_phi_reheat;
                        dVdphi = dVdphi_reheat;}

	}
	else if (inflation_model == "cubic"){

	// Hilltop Cubic
        double phi_flat = 0.0;
        double phi_reheat = 0.9*mu;
        double phi_reheat_min = 1.12305*mu;
        double Lambda_reheat = 5.44732*Lambda*pow(mu,-2.0);

        double V_of_phi_flat = Lambda;
        double dVdphi_flat = 0.0;

        Real V_of_phi_inf = Lambda*(1. - pow(mu,-3.)*pow(phi_here,3.));
        Real dVdphi_inf = -3.*Lambda*pow(mu,-3.)*pow(phi_here,2.);

        Real V_of_phi_reheat = Lambda_reheat*pow((phi_here - phi_reheat_min) / mpl, 2.0);
        Real dVdphi_reheat = 2.0 * Lambda_reheat*(phi_here - phi_reheat_min) * pow(mpl, -2.0);

		if (phi_here < phi_flat){
                        V_of_phi = V_of_phi_flat;
                        dVdphi = dVdphi_flat;}
                if (phi_here > phi_flat && phi_here < phi_reheat){
                        V_of_phi = V_of_phi_inf;
                        dVdphi = dVdphi_inf;}
                if (phi_here > phi_reheat){
                        V_of_phi = V_of_phi_reheat;
                        dVdphi = dVdphi_reheat;}

	}
	else if (inflation_model == "quartic"){

	// Hilltop Quartic
        double phi_flat = 0.0;
        double phi_reheat = 0.9*mu;
        double phi_reheat_min = 1.13587*mu;
        double Lambda_reheat = 6.18134 *Lambda*pow(mu,-2.0);

        Real V_of_phi_flat = Lambda;
        Real dVdphi_flat = 0.0;

        Real V_of_phi_inf = Lambda*(1 - pow(mu,-4)*pow(phi_here,4));
        Real dVdphi_inf = -4*Lambda*pow(mu,-4)*pow(phi_here,3);

        Real V_of_phi_reheat = Lambda_reheat*pow((phi_here - phi_reheat_min) / mpl, 2.0);
        Real dVdphi_reheat = 2.0 * Lambda_reheat*(phi_here - phi_reheat_min) * pow(mpl, -2.0);

		if (phi_here < phi_flat){
                        V_of_phi = V_of_phi_flat;
                        dVdphi = dVdphi_flat;}
                if (phi_here > phi_flat && phi_here < phi_reheat){
                        V_of_phi = V_of_phi_inf;
                        dVdphi = dVdphi_inf;}
                if (phi_here > phi_reheat){
                        V_of_phi = V_of_phi_reheat;
                        dVdphi = dVdphi_reheat;}
	}
	else if (inflation_model == "dbrane"){

	// DBrane
        double phi_reheat_min = -0.55584*mu;
        double phi_reheat = -1.2*mu;
        double Lambda_reheat = 1.24776*Lambda*pow(mu,-2.0);

        Real V_of_phi_inf = Lambda*(1. - pow(mu,4.)*pow(phi_here,-4.));
        Real dVdphi_inf = 4.*Lambda*pow(mu,4.)*pow(phi_here,-5.);

        Real V_of_phi_reheat = Lambda_reheat*pow((phi_here - phi_reheat_min) / mpl, 2.0);
        Real dVdphi_reheat = 2.0 * Lambda_reheat*(phi_here - phi_reheat_min) * pow(mpl, -2.0);

        Real V_of_phi_flat = 0.;
        Real dVdphi_flat = 0.;

        if (phi_here < phi_reheat){
                V_of_phi = V_of_phi_inf;
                dVdphi = dVdphi_inf;}
        if (phi_here > phi_reheat && phi_here < phi_reheat_min){
                V_of_phi = V_of_phi_reheat;
                dVdphi = dVdphi_reheat;}
        if (phi_here > phi_reheat_min){
                V_of_phi = V_of_phi_flat;
                dVdphi = dVdphi_flat;}
	}

}



#endif /* MYPOTENTIALFUNCTION_HPP */
