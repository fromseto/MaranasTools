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
        string workspace_name;
    } OptStoicParams;

    /*
     report_name - name of the report object that gets generated.
     report_ref - UPA of the report object that gets generated.
     */
    typedef structure {
        string report_name;
        string report_ref;
    } OptStoicOutput;

    typedef structure {
        string model_upa;
        float fixed_gr;
    } ModelInput;

    typedef structure {
        list<ModelInput> model_inputs;
        string medium_upa;
        string flux_output;
        string workspace_name;
    } SteadyComParams;

    typedef structure {
        string report_name;
        string report_ref;
        string flux_output;
    } SteadyComOutput;

    funcdef run_optstoic(OptStoicParams params) returns (OptStoicOutput output) authentication required;
    funcdef run_steadycom(SteadyComParams params) returns (SteadyComOutput) authentication required;
};
