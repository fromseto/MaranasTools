/*
A KBase module: MaranasTools
This sample module contains one small method - filter_contigs.
*/

module MaranasTools {
    /* A boolean - 0=false, 1=true
        @range (0, 1)
    */
    typedef int boolean;

    /* An X/Y/Z style reference to an FBA model.
    */
    typedef string model_upa;

    /* The id of a compound that exists either in the model or in the biochemistry.
    */
    typedef string compound_id;
	
	/* A string representing a reaction id.
    */
    typedef string reaction_id;
	
	/* A string representing a Media id.
    */
    typedef string media_id;
    
    
    
    /*
     model - the FBA model to use as a basis for modification
     start_compound - the initial compound to be used as a source for the pathway
     target_compound - the target compound to maximize yield for in the pathway
     max_steps - the maximum number of steps to allow in the optimized pathway - any pathway
                 created that has more than this number of steps is disqualified
     use_heterologous_steps - allows adding
     dG_threshold - a threshold free energy value to further constrain the path optimization
     */
    typedef structure {
        model_upa model;
        compound_id start_compound;
        compound_id target_compound;
        int max_steps;
        boolean use_heterologous_steps;
        float dG_threshold; /* advanced */
    } OptStoicParams;

    /*
     report_name - name of the report object that gets generated.
     report_ref - UPA of the report object that gets generated.
     */
    typedef structure {
        string report_name;
        string report_ref;
    } OptStoicOutput;

	/*
     model - the community FBA model
     
     */
    typedef structure {
        model_upa model;
        media_id medium;
        
        list<reaction_id> reaction_fva;
        float gr_tol; /* tolerance for convergence of max growth rate */
        list<float> opt_gr_percent; /* % of max growth rate at which FVA is performed */
        
        boolean use_random_uptake_bounds;
        float total_carb_uptake; /* total specific uptake rate of carbon in C-mmol/gDW/hr, used only if use_random_uptake_bounds is true */
        list<string> custom_bound_list;
        list<compound_id> media_supplement_list;
        boolean return_min_sum_flux;
        /* Need: additional constraints on fluxes and organism abundances,
           e.g., for implementing maximum total carbon uptake rate, molecular crowding constraints. */
    } SteadyComParams;
    
    typedef structure {
        string report_name;
        string report_ref;
        /* Need: maximum growth rate, flux distribution, organism abundance, 
        FVA ranges (at different % of maximum growth rate) */
    } SteadyComOutput;

    
    funcdef run_optstoic(OptStoicParams params) returns (OptStoicOutput output) authentication required;

	funcdef run_steadycom(SteadyComParams params) returns (SteadyComOutput output) authentication required;
};
