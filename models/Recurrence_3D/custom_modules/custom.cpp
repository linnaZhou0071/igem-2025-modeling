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

	// Enumerate loaded user double parameters for debugging
	std::cout << "[debug] Checking key trapping params: ";
	int didx = parameters.doubles.find_index("k_trap_dermis");
	int eidx = parameters.doubles.find_index("k_trap_epidermis");
	std::cout << "k_trap_dermis="; if( didx!=-1 ) std::cout << parameters.doubles("k_trap_dermis"); else std::cout << "(missing)";
	std::cout << " k_trap_epidermis="; if( eidx!=-1 ) std::cout << parameters.doubles("k_trap_epidermis"); else std::cout << "(missing)";
	std::cout << std::endl;
	// No fallback insertion here because API lacks add_variable; rely on XML correctness.
	
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

	// After initialization, populate ECM.dermis (z in [0,200)) and ECM.epidermis (z in [200,260]) mask fields.
	// These are non-diffusive / non-decaying fields we keep static at 1 inside their regions, 0 elsewhere.
	int idx_dermis = microenvironment.find_density_index("ECM.dermis");
	int idx_epidermis = microenvironment.find_density_index("ECM.epidermis");
	if( idx_dermis >= 0 || idx_epidermis >= 0 )
	{
		// Precompute region bounds
		double z_dermis_min = 0.0; double z_dermis_max = 200.0; // [0,200)
		double z_epidermis_min = 200.0; double z_epidermis_max = 260.0; // [200,260]
		int nvox = microenvironment.number_of_voxels();
		for( int n=0; n<nvox; n++ )
		{
			Voxel& vox = microenvironment.mesh.voxels[n];
			double zc = vox.center[2];
			std::vector<double>& dens = microenvironment.density_vector(n);
			if( idx_dermis >= 0 )
			{
				if( zc >= z_dermis_min && zc < z_dermis_max ) dens[idx_dermis] = 1.0; else dens[idx_dermis] = 0.0;
			}
			if( idx_epidermis >= 0 )
			{
				if( zc >= z_epidermis_min && zc <= z_epidermis_max ) dens[idx_epidermis] = 1.0; else dens[idx_epidermis] = 0.0;
			}
		}
		// Carve sebaceous gland void: now a 2-phase geometry
		// Phase 1 (z_base -> join_z): frustum (linear radius from r_base to r_join)
		// Phase 2 (z_base -> join_z): spherical cap closure below frustum to mimic bulbous closed base
		// Above join_z up to z_top: continue linear taper to r_top (opening)
		// Parameters from user_parameters
		struct GlandGeom { 
			double cx=0, cy=0; 
			double z_base=40, z_top=260; 
			double r_base=80, r_top=40; 
			double z_join=120; 
			double cap_R=90; // sphere radius for cap closure
			bool taper=true; 
		};
		GlandGeom gg;
		int id;
		if( (id=parameters.doubles.find_index("gland_center_x"))!=-1 ) gg.cx = parameters.doubles("gland_center_x");
		if( (id=parameters.doubles.find_index("gland_center_y"))!=-1 ) gg.cy = parameters.doubles("gland_center_y");
		if( (id=parameters.doubles.find_index("gland_base_z"))!=-1 ) gg.z_base = parameters.doubles("gland_base_z");
		if( (id=parameters.doubles.find_index("gland_top_z"))!=-1 ) gg.z_top = parameters.doubles("gland_top_z");
		if( (id=parameters.doubles.find_index("gland_radius_base"))!=-1 ) gg.r_base = parameters.doubles("gland_radius_base");
		if( (id=parameters.doubles.find_index("gland_radius_top"))!=-1 ) gg.r_top = parameters.doubles("gland_radius_top");
		int ib = parameters.bools.find_index("gland_taper"); if( ib != -1 ) gg.taper = parameters.bools("gland_taper");
		if( (id=parameters.doubles.find_index("gland_join_z"))!=-1 ) gg.z_join = parameters.doubles("gland_join_z");
		if( (id=parameters.doubles.find_index("gland_cap_radius"))!=-1 ) gg.cap_R = parameters.doubles("gland_cap_radius");
		if( gg.z_top <= gg.z_base ){ gg.z_top = gg.z_base + 1.0; }
		if( gg.z_join <= gg.z_base ) gg.z_join = gg.z_base + 0.5*(gg.z_top-gg.z_base); // fallback
		if( gg.z_join > gg.z_top ) gg.z_join = gg.z_top - 0.5*(gg.z_top-gg.z_base);
		if( gg.cap_R < 1.0 ) gg.cap_R = gg.r_base; // sanity
		// Precompute inverse height for taper
		double inv_h_total = 1.0 / (gg.z_top - gg.z_base);
		double inv_h_join = 1.0 / (gg.z_join - gg.z_base);
		int cleared = 0;
		for( int n=0; n<nvox; n++ )
		{
			Voxel& vox = microenvironment.mesh.voxels[n];
			double x=vox.center[0], y=vox.center[1], z=vox.center[2];
			if( z < gg.z_base || z > gg.z_top ) continue; // outside vertical span
			// radius logic:
			double r_linear = gg.r_base; // linear interpolation whole height (for above join)
			if( gg.taper )
			{
				double f_lin = (z - gg.z_base) * inv_h_total; if( f_lin<0 ) f_lin=0; if( f_lin>1 ) f_lin=1;
				r_linear = gg.r_base + (gg.r_top - gg.r_base) * f_lin;
			}
			double r_here = r_linear;
			if( z <= gg.z_join )
			{
				// Below join: combine frustum shrink (within base->join) with spherical cap giving extra outward bulge
				double f_join = (z - gg.z_base) * inv_h_join; if( f_join<0 ) f_join=0; if( f_join>1 ) f_join=1;
				// Frustum radius at z within base->join (assume target radius at join is r_linear_at_join)
				double r_join_linear = gg.r_base + ( (r_linear) - gg.r_base ) * f_join; // reuse r_linear projected
				// Spherical cap center (approx) placed below base so that at z=z_base radius roughly r_base
				double z_cap_center = gg.z_base - (gg.cap_R - 0.0); // sphere center below base
				double dz_cap = z - z_cap_center; if( dz_cap < 0 ) dz_cap = 0; // inside sphere above center
				double r_cap = sqrt( std::max(0.0, gg.cap_R*gg.cap_R - dz_cap*dz_cap ) );
				// Blend: take max of frustum and spherical cap to get bulbous shape
				r_here = std::max(r_join_linear , r_cap);
			}
			double dx = x - gg.cx; double dy = y - gg.cy;
			double r2 = dx*dx + dy*dy;
			if( r2 <= r_here*r_here )
			{
				// Inside gland volume: zero ECM masks
				std::vector<double>& dens = microenvironment.density_vector(n);
				if( idx_dermis >=0 ) dens[idx_dermis] = 0.0;
				if( idx_epidermis >=0 ) dens[idx_epidermis] = 0.0;
				cleared++;
			}
		}
		std::cout << "setup_microenvironment: ECM.dermis / ECM.epidermis masks initialized with sebaceous gland void (cleared voxels=" << cleared << ")." << std::endl;
	}

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

// Coloring disabled: use default paint_by_number only.
std::vector<std::string> my_coloring_function( Cell* pCell )
{
	return paint_by_number_cell_coloring( pCell );
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
// Forward declaration for ECM trapping helper (defined later in file)
void adjust_trapping_with_ECM_regions( double dt );

void update_substrate_reactions( double dt )
{
	// One-time debug of trapping parameters presence
	static bool printed_params = false;
	if( !printed_params )
	{
		int idx_dermis_param = parameters.doubles.find_index("k_trap_dermis");
		int idx_epi_param = parameters.doubles.find_index("k_trap_epidermis");
		int idx_epi_typo = parameters.doubles.find_index("k_trap_epidermi"); // possible user typo
		std::cout << "[debug] k_trap_dermis index=" << idx_dermis_param
		          << " k_trap_epidermis index=" << idx_epi_param
		          << " (typo k_trap_epidermi index=" << idx_epi_typo << ")" << std::endl;
		if( idx_dermis_param != -1 )
			std::cout << "[debug] k_trap_dermis value=" << parameters.doubles("k_trap_dermis") << std::endl;
		if( idx_epi_param != -1 )
			std::cout << "[debug] k_trap_epidermis value=" << parameters.doubles("k_trap_epidermis") << std::endl;
		else if( idx_epi_typo != -1 )
			std::cout << "[debug] Using typo parameter k_trap_epidermi value=" << parameters.doubles("k_trap_epidermi") << std::endl;
		printed_params = true;
	}
	// EPS-based immobilization (legacy simplified logic)
	adjust_diffusivity_with_EPS( dt );
	// ECM-based trapping (new): uniformly trap all mobile substrates in dermis / epidermis masks.
	adjust_trapping_with_ECM_regions( dt );
}

// Approximate a local, per-voxel reduction in effective diffusivity by
// adding a local EPS-dependent binding/immobilization (a local sink) for
// AMP and bacteriocin. This reduces the mobile concentration locally and
// therefore mimics reduced spread without changing global BioFVM D values.
void adjust_diffusivity_with_EPS( double dt )
{
	int idx_AMP = microenvironment.find_density_index("AMP");
	int idx_EPS = microenvironment.find_density_index("EPS");
	if( idx_AMP < 0 || idx_EPS < 0 ) { return; }

	// Parameters for AMP immobilization by EPS
	double k_bind_amp = 0.5; // default
	double max_fraction_removed = 0.9; // stability cap
	if( parameters.doubles.find_index("k_bind_amp") != -1 )
		k_bind_amp = parameters.doubles("k_bind_amp");
	if( parameters.doubles.find_index("max_fraction_removed") != -1 )
		max_fraction_removed = parameters.doubles("max_fraction_removed");

	int idx_AMP_bound = microenvironment.find_density_index("AMP_bound");
	int nvoxels = microenvironment.number_of_voxels();
	for( int n = 0; n < nvoxels; n++ )
	{
		std::vector<double>& densities = microenvironment.density_vector(n);
		double eps = densities[idx_EPS];
		if( eps <= 0 ) continue;
		double amp = densities[idx_AMP];
		if( amp <= 0 ) continue;
		double uptake = k_bind_amp * eps * amp * dt;
		double max_uptake = max_fraction_removed * amp;
		if( uptake > max_uptake ) uptake = max_uptake;
		densities[idx_AMP] -= uptake;
		if( densities[idx_AMP] < 0 ) densities[idx_AMP] = 0.0;
		if( idx_AMP_bound >= 0 ) densities[idx_AMP_bound] += uptake;
	}
}

// Generic trapping in ECM.dermis / ECM.epidermis regions for all substrates (except already-bound / non-diffusive fields).
// We model trapping as: dS = -k_trap_region * region_mask * S * dt, optionally accumulating in a corresponding bound field if it exists.
// Currently we use same k across all substrates inside a region (separate k for dermis vs epidermis from user parameters).
void adjust_trapping_with_ECM_regions( double dt )
{
	int idx_dermis = microenvironment.find_density_index("ECM.dermis");
	int idx_epidermis = microenvironment.find_density_index("ECM.epidermis");
	if( idx_dermis < 0 && idx_epidermis < 0 ) return; // nothing to do

	// Fetch unified region-level rates
	int idx_kd = parameters.doubles.find_index("k_trap_dermis");
	int idx_ke = parameters.doubles.find_index("k_trap_epidermis");
	double k_dermis = (idx_kd!=-1)? parameters.doubles("k_trap_dermis") : 0.5;
	double k_epidermis = (idx_ke!=-1)? parameters.doubles("k_trap_epidermis") : 0.5;
	double max_fraction_trapped = 0.9; int idx_cap = parameters.doubles.find_index("max_fraction_trapped_ECM"); if( idx_cap!=-1 ) max_fraction_trapped = parameters.doubles("max_fraction_trapped_ECM"); if( max_fraction_trapped<=0 ) max_fraction_trapped=1.0; if( max_fraction_trapped>1 ) max_fraction_trapped=1.0;

	// List substrates to apply trapping (skip masks and bound immobile pools except allow EPS)
	std::vector<std::string> substrate_names = { "sebum","AMP","Th2_cytokine","AIP_sa","EPS","toxin" };
	struct SubIdx { int i; };
	std::vector<int> indices; indices.reserve(substrate_names.size());
	for( auto &nm : substrate_names ){
		int si = microenvironment.find_density_index(nm); if( si<0 ) continue; if( microenvironment.diffusion_coefficients[si]<=0 && nm!="EPS") continue; indices.push_back(si);
	}
	int nvox = microenvironment.number_of_voxels();
	for( int n=0; n<nvox; ++n ){
		std::vector<double>& dens = microenvironment.density_vector(n);
		double dermis_mask    = (idx_dermis>=0)? dens[idx_dermis] : 0.0;
		double epidermis_mask = (idx_epidermis>=0)? dens[idx_epidermis] : 0.0;
		if( dermis_mask<=0 && epidermis_mask<=0 ) continue;
		double k_local = 0.0; if( dermis_mask>0 ) k_local += k_dermis; if( epidermis_mask>0 ) k_local += k_epidermis; if( k_local<=0 ) continue;
		for( int si : indices ){
			double &S = dens[si]; if( S<=0 ) continue; double removal = k_local * S * dt; double max_removal = max_fraction_trapped * S; if( removal>max_removal ) removal=max_removal; if( removal>S ) removal=S; S -= removal; }
	}
}