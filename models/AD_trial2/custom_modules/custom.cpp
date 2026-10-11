/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cctype>
#include <mutex>
#include <unordered_set>

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	// Revert to default internal phenotype update so that death / volume progression occurs.
	// (Previously we overrode with an empty phenotype_function, which froze volume changes
	//  and could prevent volume-based removal behavior.)
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		int n_cells = 0;
		if( parameters.ints.find_index("number_of_cells") != -1 )
		{
			n_cells = parameters.ints("number_of_cells");
		}
		else
		{
			std::cout << "Note: parameter 'number_of_cells' not found; using 0 (will rely on CSV load instead)." << std::endl;
			n_cells = 0;
		}

		for( int n = 0 ; n < n_cells ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();

	// Safety: ensure every existing cell (including those loaded from CSV) has the custom rule
	// In some setups, function pointers may not propagate back to already-created cells.
	int reassigned = 0;
	for( Cell* pC : *all_cells )
	{
		if( pC->functions.custom_cell_rule != custom_function )
		{
			pC->functions.custom_cell_rule = custom_function;
			reassigned++;
		}
	}
	if( reassigned > 0 )
	{
		std::cout << "setup_tissue: reassigned custom_cell_rule to " << reassigned << " loaded cells" << std::endl;
	}
	
	return; 
}

static inline void trim_inplace(std::string &s)
{
	size_t a = s.find_first_not_of(" \t\r\n");
	if (a == std::string::npos) { s.clear(); return; }
	size_t b = s.find_last_not_of(" \t\r\n");
	s = s.substr(a, b - a + 1);
}

// Loadable coloring function: looks up a color map file first, falls back to paint_by_number
std::vector<std::string> my_coloring_function( Cell* pCell )
{
	static bool loaded = false;
	static std::unordered_map<std::string,std::string> cmap;

	if( ! loaded )
	{
		std::string fname = "cell_color_map.txt"; // simple key=value per line: celltype=color
		std::ifstream ifs( fname );
		static bool loaded_msg = false;
		if( ! ifs )
		{
			// try the output folder (PhysiCell_settings.folder)
			std::string alt = PhysiCell_settings.folder + "/cell_color_map.txt";
			ifs.open( alt );
			if( ifs && ! loaded_msg )

				{
					std::cout << "my_coloring_function: loaded color map from " << alt << "\n";
					loaded_msg = true;
			}
		}
		else if( ! loaded_msg )
		{
			std::cout << "my_coloring_function: loaded color map from ./cell_color_map.txt\n";
			loaded_msg = true;
		}

	if( ifs )
		{
			std::string line;
			while( std::getline( ifs, line ) )
			{
				// strip comments
				size_t com = line.find('#'); if( com != std::string::npos ) line = line.substr(0, com);
				trim_inplace(line);
				if( line.empty() ) continue;
				size_t eq = line.find('='); if( eq == std::string::npos ) continue;
						std::string name = line.substr(0, eq);
						std::string color = line.substr(eq+1);
						trim_inplace(name); trim_inplace(color);
						if( ! name.empty() && ! color.empty() ) {
							std::string name_l = name;
							std::transform(name_l.begin(), name_l.end(), name_l.begin(), [](unsigned char c){ return std::tolower(c); });
							cmap[name_l] = color;
						}
			}
			ifs.close();
		}
		else
		{
			std::cout << "my_coloring_function: no cell_color_map.txt found in current dir or " << PhysiCell_settings.folder << "\n";
		}

		// sensible defaults if the user didn't provide mapping
			if( cmap.find("s.aurenus") == cmap.end() ) cmap["s.aurenus"] = "rgb(212,175,55)";
		    if( cmap.find("s.epidermis") == cmap.end() ) cmap["s.epidermis"] = "rgb(35, 74, 231)";
			if( cmap.find("sebocyte") == cmap.end() ) cmap["sebocyte"] = "rgb(225, 113, 21)"; 
			if( cmap.find("s.hominis") == cmap.end() ) cmap["s.hominis"] = "rgb(207, 16, 16)";
			if( cmap.find("other commensals") == cmap.end() ) cmap["other commensals"] = "rgb(38, 243, 27)";


		loaded = true;
	}

	// safe default: if no mapping, fallback to paint_by_number
	if( pCell == nullptr ) return paint_by_number_cell_coloring(pCell);

	// ALWAYS override with death colors first so dead cells appear brown and not the living type color
	// (PhysiCell Studio default behavior). We intentionally avoid user map so death is obvious.
	if( pCell->phenotype.death.dead )
	{
		// Use a standard necrotic brown. Could differentiate apoptosis vs necrosis if desired.
		return { "saddlebrown", "black", "saddlebrown", "black" };
	}

	std::string type_name = "";
	if( pCell->type >= 0 && pCell->type < (int)cell_definitions_by_index.size() )
	{
		type_name = cell_definitions_by_index[ pCell->type ]->name;
	}
	// perform case-insensitive lookup
	std::string type_name_l = type_name;
	std::transform(type_name_l.begin(), type_name_l.end(), type_name_l.begin(), [](unsigned char c){ return std::tolower(c); });

	auto it = cmap.find( type_name_l );
	if( it != cmap.end() )
	{
		std::string interior_color = it->second;
		std::vector<std::string> output = { interior_color, "black", interior_color, "black" };

		// print using-color message only once per cell type to avoid flooding stdout
		static std::mutex report_mutex;
		static std::unordered_set<std::string> reported_types;
		{
			std::lock_guard<std::mutex> g(report_mutex);
			if( reported_types.find(type_name_l) == reported_types.end() )
			{
				reported_types.insert(type_name_l);
				std::cout << "my_coloring_function: using color '" << interior_color << "' for '" << type_name << "'\n";
			}
		}

		return output;
	}

	// otherwise use existing paint_by_number behavior
	return paint_by_number_cell_coloring(pCell);
}

// No custom phenotype update needed (left intentionally empty / unused).
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	// Simplified logic: if immediate_remove_small_microbes == true -> any dead cell removed immediately.
	// Otherwise: remove only when volume.total < removal_volume_threshold.
	struct RemovalConfig { double threshold = 1.0; bool immediate_remove = true; };
	static RemovalConfig cfg; static std::once_flag once;
	std::call_once( once, [](){
		int idx = parameters.doubles.find_index("removal_volume_threshold"); if( idx != -1 ) cfg.threshold = parameters.doubles("removal_volume_threshold");
		int bidx = parameters.bools.find_index("immediate_remove_small_microbes"); if( bidx != -1 ) cfg.immediate_remove = parameters.bools("immediate_remove_small_microbes");
		std::cout << "[custom_function init] threshold=" << cfg.threshold
		          << " immediate_remove_small_microbes=" << (cfg.immediate_remove?"true":"false") << "\n"; } );

	if( !phenotype.death.dead ) return; // alive
	if( cfg.immediate_remove ) { pCell->flag_for_removal(); return; }
	if( phenotype.volume.total < cfg.threshold ) { pCell->flag_for_removal(); return; }
	return; 
} 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

// baseline and interaction reactions among substrates
void update_substrate_reactions( double dt )
{
	// parameter names expected in XML (config/settings)
	static bool params_loaded = false;
	static double AMP_baseline_rate = 0.0; // 1/min per volume
	static double Cytokine_baseline_rate = 0.0;
	static double toxin_to_cytokine_rate = 0.0; // 1/min * [toxin]
	static double cytokine_to_AMP_degradation_rate = 0.0; // 1/min * [cytokine]

	if( ! params_loaded )
	{
		if( parameters.doubles.find_index("AMP_baseline_rate") != -1 )
		{ AMP_baseline_rate = parameters.doubles("AMP_baseline_rate"); }
		if( parameters.doubles.find_index("Cytokine_baseline_rate") != -1 )
		{ Cytokine_baseline_rate = parameters.doubles("Cytokine_baseline_rate"); }
		if( parameters.doubles.find_index("toxin_to_cytokine_rate") != -1 )
		{ toxin_to_cytokine_rate = parameters.doubles("toxin_to_cytokine_rate"); }
		if( parameters.doubles.find_index("cytokine_to_AMP_degradation_rate") != -1 )
		{ cytokine_to_AMP_degradation_rate = parameters.doubles("cytokine_to_AMP_degradation_rate"); }

		params_loaded = true;
	}

	// substrate indices: use names from microenvironment (from settings)
	int idx_AMP = microenvironment.find_density_index("AMP");
	int idx_cyt = microenvironment.find_density_index("Th2_cytokine");
	int idx_toxin = microenvironment.find_density_index("toxin");

	if( idx_AMP < 0 || idx_cyt < 0 || idx_toxin < 0 )
	{ return; }

	int nvoxels = microenvironment.number_of_voxels();

	for( int n = 0 ; n < nvoxels ; n++ )
	{
		std::vector<double> & densities = microenvironment.density_vector( n );

		double toxin = densities[idx_toxin];
		double cytokine = densities[idx_cyt];
		double amp = densities[idx_AMP];

		// baseline production: add amount = rate * dt
		densities[idx_AMP] += AMP_baseline_rate * dt;
		densities[idx_cyt] += Cytokine_baseline_rate * dt;

		// toxin increases cytokine: simple linear production proportional to toxin
		if( toxin > 0 )
		{
			double add_cyt = toxin_to_cytokine_rate * toxin * dt;
			densities[idx_cyt] += add_cyt;
		}

		// cytokine decreases AMP: simple linear degradation proportional to cytokine
		if( cytokine > 0 && amp > 0 )
		{
			double degrade_amp = cytokine_to_AMP_degradation_rate * cytokine * dt;
			// ensure concentration doesn't go negative
			densities[idx_AMP] -= degrade_amp;
			if( densities[idx_AMP] < 0 ) densities[idx_AMP] = 0.0;
		}
	}

	// after reactions, adjust diffusivity of AMP and bacteriocin based on EPS
	adjust_diffusivity_with_EPS( dt );
}

// Approximate a local, per-voxel reduction in effective diffusivity by
// adding a local EPS-dependent binding/immobilization (a local sink) for
// AMP and bacteriocin. This reduces the mobile concentration locally and
// therefore mimics reduced spread without changing global BioFVM D values.
void adjust_diffusivity_with_EPS( double dt )
{
	int idx_AMP = microenvironment.find_density_index("AMP");
	int idx_bac = microenvironment.find_density_index("bacteriocin");
	int idx_EPS = microenvironment.find_density_index("EPS");

	if( idx_AMP < 0 || idx_bac < 0 || idx_EPS < 0 )
	{
		return;
	}

	// Parameters controlling how strongly EPS immobilizes substrates.
	// k_bind has units 1/(concentration * time) so that uptake = k_bind * EPS * conc * dt
	double k_bind_amp = 0.5;     // default if not set in XML
	double k_bind_bac = 0.5;     // default if not set in XML
	double max_fraction_removed = 0.9; // cap per time-step fraction removed to keep stability

	if( parameters.doubles.find_index("k_bind_amp") != -1 )
	{ k_bind_amp = parameters.doubles("k_bind_amp"); }
	if( parameters.doubles.find_index("k_bind_bac") != -1 )
	{ k_bind_bac = parameters.doubles("k_bind_bac"); }
	if( parameters.doubles.find_index("max_fraction_removed") != -1 )
	{ max_fraction_removed = parameters.doubles("max_fraction_removed"); }

	int nvoxels = microenvironment.number_of_voxels();

	int idx_AMP_bound = microenvironment.find_density_index("AMP_bound");
	int idx_bac_bound = microenvironment.find_density_index("bacteriocin_bound");

    

	for( int n = 0 ; n < nvoxels ; n++ )
	{
		std::vector<double> & densities = microenvironment.density_vector( n );

		double eps = densities[idx_EPS];

		if( eps <= 0 ) continue; // no EPS -> no effect

		// AMP
		double amp = densities[idx_AMP];
		if( amp > 0 )
		{
			double uptake = k_bind_amp * eps * amp * dt;
			// cap uptake so we don't remove more than a fraction of the voxel's content
			double max_uptake = max_fraction_removed * amp;
			if( uptake > max_uptake ) uptake = max_uptake;
			densities[idx_AMP] -= uptake;
			if( densities[idx_AMP] < 0 ) densities[idx_AMP] = 0.0;
			if( idx_AMP_bound >= 0 )
			{
				densities[idx_AMP_bound] += uptake; // accumulate bound fraction
			}
		}

		// bacteriocin
		double bac = densities[idx_bac];
		if( bac > 0 )
		{
			double uptake = k_bind_bac * eps * bac * dt;
			double max_uptake = max_fraction_removed * bac;
			if( uptake > max_uptake ) uptake = max_uptake;
			densities[idx_bac] -= uptake;
			if( densities[idx_bac] < 0 ) densities[idx_bac] = 0.0;
			if( idx_bac_bound >= 0 )
			{
				densities[idx_bac_bound] += uptake; // accumulate bound fraction
			}
		}
	}

    
}