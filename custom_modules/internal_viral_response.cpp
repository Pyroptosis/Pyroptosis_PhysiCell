#include "./internal_viral_response.h" 

using namespace PhysiCell; 

std::string internal_virus_response_version = "0.2.0"; 

Submodel_Information internal_virus_response_model_info; 

void internal_virus_response_model_setup( void )
{
	// set up the model 
		// set version info 
	internal_virus_response_model_info.name = "internal viral response"; 
	internal_virus_response_model_info.version = internal_virus_response_version; 
		// set functions 
	internal_virus_response_model_info.main_function = NULL; 
	internal_virus_response_model_info.phenotype_function = internal_virus_response_model; 
	internal_virus_response_model_info.mechanics_function = NULL; 
	
		// what microenvironment variables do you expect? 
		// what custom data do I need? 
	internal_virus_response_model_info.cell_variables.push_back( "max_infected_apoptosis_rate" ); 
	internal_virus_response_model_info.cell_variables.push_back( "max_apoptosis_half_max" ); 
	internal_virus_response_model_info.cell_variables.push_back( "apoptosis_hill_power" ); 	
	
		// register the submodel  
	internal_virus_response_model_info.register_model();	
		// set functions for the corresponding cell definition 
		
//	pCD = find_cell_definition( "lung epithelium" ); 
//	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
//	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	
	return; 
}

void internal_virus_response_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	static Cell_Definition* pCD = find_cell_definition( "lung epithelium" ); 
	
if( pCell->custom_data["cell_virus_induced_apoptosis_flag"]>=1 )
{
	pCell->custom_data["cell_apo_time"]=pCell->custom_data["cell_apo_time"]+dt;
}
	
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 
	
	

	// now, set apoptosis rate 
	
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "apoptosis" );
	

	
	
	
	// if we're infected, secrete a chemokine for the immune model
	static int nAV = pCell->custom_data.find_variable_index( "assembled_virion" ); 	
	double AV = pCell->custom_data[nAV]; 
	
	static int nP = pCell->custom_data.find_variable_index( "viral_protein"); 
	double P = pCell->custom_data[nP];

	static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
	

	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
		
	static int nR = pCell->custom_data.find_variable_index("viral_RNA");
	double R = pCell->custom_data[nR];
	
	if( R >= 1.00 - 1e-16 ) 
	{
		pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 1.0; 
	}

	if( pCell->custom_data["infected_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
	{
		double rate = AV; 
		rate /= pCell->custom_data["max_apoptosis_half_max"];
		if( rate > 1.0 )
		{ rate = 1.0; }
		rate *= pCell->custom_data[ "infected_cell_chemokine_secretion_rate" ];

		phenotype.secretion.secretion_rates[chemokine_index] = rate; 
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0;

		// (Adrianne) adding pro-inflammatory cytokine secretion by infected cells
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = pCell->custom_data["activated_cytokine_secretion_rate"];
	}
	

	if( pCell->custom_data["cell_pyroptosis_flag"]>=1 )
	{
		pyroptosis_cascade( pCell, phenotype, dt ); 
		return;
	};
	
	// Sara&Fiona: Internalise pro-pyroptotic cytokine from the microenvironment
	static int pro_pyroptotic_cytokine = microenvironment.find_density_index("pro-pyroptosis cytokine"); 
	double pyroptotic_cytokine_concentration = phenotype.molecular.internalized_total_substrates[pro_pyroptotic_cytokine]; 


    // (Sara&Fiona) pyroptosis cascade in the cell is initiated if cell's viral_RNA is >1 (i.e. >=3). This is arbitraty to check things work.
	if( R>2 && pCell->custom_data["cell_pyroptosis_flag"]==0 && pCell->custom_data["cell_virus_induced_apoptosis_flag"]==0 && pCell->custom_data["cell_bystander_pyroptosis_flag"]==0)
	{
		// set the probability (in 0,100) that a cell with a death-sentence pyroptoses (not apoptoses)
		int cell_death_pyroptosis_probability = 100; 
		// randomise a number in 0,100 that determines the cell death mode (pyroptosis or apoptosis)
		int cell_death_dice =  rand() % 101;
		if(cell_death_dice < cell_death_pyroptosis_probability) 
		{
			pCell->custom_data["cell_pyroptosis_flag"]=1; 
			pCell->custom_data["cell_pyroptosis_time"]=1;
			//cell pyroptoses
		}
		else 
		{
			pCell->custom_data["cell_virus_induced_apoptosis_flag"]=1; //cell apoptoses
			phenotype.death.rates[apoptosis_model_index] = 1; 
			
			
			
		}


		
		return;
	}
    // (Sara&Fiona) % base case 200
    else if(R<2 && pyroptotic_cytokine_concentration>200.0 && pCell->custom_data["cell_pyroptosis_flag"]==0 && pCell->custom_data["cell_virus_induced_apoptosis_flag"]==0 && pCell->custom_data["cell_bystander_pyroptosis_flag"]==0) 
    {
		pCell->custom_data["cell_pyroptosis_flag"]=2; // Pyroptosis cascade is initiated
		pCell->custom_data["cell_bystander_pyroptosis_flag"]=1; // Pyroptosis cascade is initiated
		pCell->custom_data["cell_pyroptosis_time"]=1;
		
		//printf("Pyro bystander effect!\n");
		return;
    }
	
	return; 
	

}



// Sara&Fiona: code for pyroptosis
void pyroptosis_cascade( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	pCell->custom_data["cell_pyroptosis_time"]=pCell->custom_data["cell_pyroptosis_time"]+dt;
	// Sara&Fiona: Pyroptosis code starts here.
	
	// Intracellular components
	double nfkb_n = pCell->custom_data.find_variable_index( "nuclear_NFkB" ); 
	double nlrp3_i = pCell->custom_data.find_variable_index( "inactive_NLRP3" ); 
	double nlrp3_a = pCell->custom_data.find_variable_index( "active_NLRP3" ); 
	double nlrp3_b = pCell->custom_data.find_variable_index( "bound_NLRP3" ); 
	double asc_b = pCell->custom_data.find_variable_index( "bound_ASC" );
	double caspase1_b = pCell->custom_data.find_variable_index( "bound_caspase1" );
	double gsdmd_c = pCell->custom_data.find_variable_index( "cleaved_gasderminD" );
	double il_1b_p = pCell->custom_data.find_variable_index( "pro_IL_1b" );
	double il_1b_c = pCell->custom_data.find_variable_index( "cytoplasmic_IL_1b" );
	double il_1b_e = pCell->custom_data.find_variable_index( "external_IL_1b" );
	double il_18_c = pCell->custom_data.find_variable_index( "cytoplasmic_IL_18" );
	double il_18_e = pCell->custom_data.find_variable_index( "external_IL_18" );
	double volume_c = pCell->custom_data.find_variable_index( "cytoplasmic_volume" );



    double k_nlrp3_ita =  0.7; //k1 in matlab
    double k_nlrp3_atb = 1; //k2 in matlab
    double k_asc_ftb = 0.04; //k3 in matlab
    double k_c1_ftb = 0.03; //k4 in matlab
    double k_il1b_cte = 1; //k5 in matlab
    double k_il18_cte = 1; //k6 in matlab
    double k_vol_c = 0.2; //k7 in matlab
    // Decay constants
    double d_nlrp3 = 0.002; // delta1 matlab
    double d_il =0.004; //delta2 in matlab
    // Hill function rates
    double a_nlrp3 = 0.07; //alpha1 in matlab
    double a_gsdmd =0.1; // alpha2 in matlab
	double a_il1b_p = 0.06; //alpha3 in matlab
    double a_il1b_c = 1; //alpha4 in matlab
    double a_il18 = 1;//alpha5 in matlab
	double hm_nfkb = 0.3; // NF_50 in matlab
    double hm_c1 = 0.3; // C1_50 in matlab
    double hex_nfkb = 2.0; //gammaNF in matlab
    double hex_c1 = 2.0 ; //gammaC1 in matlab

	double no_a=-1; // a in matlab
	double no_b=2;// b in matlab
	double no_c=1000; // c in matlab

	double nf_hh=0.55; //hh in matlab
	double nf_orig=0.25; //nfkb0 in matlab (can this be called in from external)
	double nf_ss=0.8; // ss in matlab
	double nf_tt=10; // tt in matlab
	double t_i=pCell->custom_data["cell_pyroptosis_time"]-dt;	// t in matlab

	
	
	//Update nuclear NFkB (updated backward) // check how to call in initial condition
	pCell->custom_data[nfkb_n]= (nf_hh)*exp(-(log((t_i)/nf_tt)*log((t_i)/nf_tt))/nf_ss)+nf_orig;
	

	//Set Hill function 1
	double hill_nfkb = (pow((pCell->custom_data[nfkb_n]-nf_orig),hex_nfkb))/(pow(hm_nfkb,hex_nfkb)+pow((pCell->custom_data[nfkb_n]-nf_orig),hex_nfkb));
	
	//Update NLRP3 (inactive, active and bound) (updated backward)
 	pCell->custom_data[nlrp3_i] = ((pCell->custom_data[nlrp3_i])+dt*a_nlrp3*hill_nfkb)/(1+dt*k_nlrp3_ita+dt*d_nlrp3);		
	
	
	pCell->custom_data[nlrp3_a]=(-(1+dt*d_nlrp3)+sqrt(((1+dt*d_nlrp3)*(1+dt*d_nlrp3))+4*(dt*k_nlrp3_atb)*(pCell->custom_data[nlrp3_a]+dt*k_nlrp3_ita*pCell->custom_data[nlrp3_i])))/(2*dt*k_nlrp3_atb);
	
	
	pCell->custom_data[nlrp3_b] = pCell->custom_data[nlrp3_b] + dt * k_nlrp3_atb * pCell->custom_data[nlrp3_a]* pCell->custom_data[nlrp3_a];

	
	double F_no = 1/(1+pow(((pCell->custom_data[nlrp3_b]-no_a)/no_b),-no_c)); // define new function for ASC binding
	//Update bound ASC (updated backward)
	pCell->custom_data[asc_b]=(pCell->custom_data[asc_b]+dt*k_asc_ftb*F_no*pCell->custom_data[nlrp3_b])/(1+dt*k_asc_ftb*F_no*pCell->custom_data[nlrp3_b]);
	
		
	//Update bound caspase1 (updated backward)
	
	pCell->custom_data[caspase1_b]=(pCell->custom_data[caspase1_b]+dt*k_c1_ftb*pCell->custom_data[asc_b])/(1+dt*k_c1_ftb*pCell->custom_data[asc_b]);
	
	//Set Hill function 2
	double hill_c1 = (pow(pCell->custom_data[caspase1_b],hex_c1))/(pow(hm_c1,hex_c1)+pow(pCell->custom_data[caspase1_b],hex_c1));
		
	
	//Update cleaved GSDMD (updated backward)
	pCell->custom_data[gsdmd_c] = (pCell->custom_data[gsdmd_c]+dt*a_gsdmd*hill_c1)/(1+dt*a_gsdmd*hill_c1);

	//Set G function (same now that total GSDMD concentration is 1 au of concentration)
	double g_gsdmd = pCell->custom_data[gsdmd_c]/1;
	
	
	//Update IL1b (pro, cytoplasmic, external)	We want to relate this to secreted cytokine IL1b (updated backward)
	pCell->custom_data[il_1b_p] = (pCell->custom_data[il_1b_p]+dt*a_il1b_p*hill_nfkb)/(1+dt*a_il1b_c*hill_c1+dt*d_il);
	pCell->custom_data[il_1b_c] = (pCell->custom_data[il_1b_c]+dt*a_il1b_c*hill_c1*(pCell->custom_data[il_1b_p]))/(1+dt*d_il+dt*k_il1b_cte*g_gsdmd);		
	pCell->custom_data[il_1b_e] = pCell->custom_data[il_1b_e] + dt * (k_il1b_cte*g_gsdmd*pCell->custom_data[il_1b_c]);

	// checked to here	
	
	//Update IL18 (cytoplasmic, external)(updated backward)
	pCell->custom_data[il_18_c] = (pCell->custom_data[il_18_c]+dt*a_il18*hill_c1*(1-pCell->custom_data[il_18_e]))/((1+dt*a_il18*hill_c1)*(1+dt*k_il18_cte*g_gsdmd));
	pCell->custom_data[il_18_e] = pCell->custom_data[il_18_e] +  dt * k_il18_cte*g_gsdmd*pCell->custom_data[il_18_c];

	//Update cytoplasmic volume (updated backward)
	pCell->custom_data[volume_c] = pCell->custom_data[volume_c]/(1-dt * k_vol_c * g_gsdmd);

	// (Yafei) need to update the real radius 
	phenotype.volume.total = pCell->custom_data[volume_c]; 


	//Temporary: "super fast" apoptosis occurs when cell should burst. 
	//To do: We actually want the cell to rupture once a cytoplasmic critical volume is reached (e.g. 1.5 of initial cytoplasmic volume from in vitro data). 
	
	static double initial_total_volume = 2494;
	if( pCell->custom_data[volume_c]>1.5*2494)
	{
		static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "apoptosis" );
		pCell->custom_data["cell_pyroptosis_flag"]=3;  

		//std::cout<<pCell->custom_data["cell_pyroptosis_time"]<<std::endl;
		//The cell's 'apoptosis death rate' is set to be "super high" 
		phenotype.death.rates[apoptosis_model_index] = 1; 
	}

	// (Adrianne) update cell pro-inflammatory secretion rate based on IL18 secretion rate - need to double check unit conversion
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = pCell->custom_data["activated_cytokine_secretion_rate"]+k_il18_cte*g_gsdmd*pCell->custom_data[il_18_c];
    // (Sara and Fiona)
	static int propyroptotic_cytokine_index = microenvironment.find_density_index("pro-pyroptosis cytokine");
	pCell->phenotype.secretion.secretion_rates[propyroptotic_cytokine_index] = k_il1b_cte*g_gsdmd*pCell->custom_data[il_1b_c];

	//To do: We also want IL-1beta (pCell->custom_data[il_1b_e]) to be secreted. 
	//IL-1beta can induce pyroptosis in other cells (pyroptosis bystander effect). 
	
	return;
}
