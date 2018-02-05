
package us.kbase.maranastools;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: SteadyComParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "model_inputs",
    "medium_upa",
    "flux_output",
    "workspace_name"
})
public class SteadyComParams {

    @JsonProperty("model_inputs")
    private List<ModelInput> modelInputs;
    @JsonProperty("medium_upa")
    private String mediumUpa;
    @JsonProperty("flux_output")
    private String fluxOutput;
    @JsonProperty("workspace_name")
    private String workspaceName;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("model_inputs")
    public List<ModelInput> getModelInputs() {
        return modelInputs;
    }

    @JsonProperty("model_inputs")
    public void setModelInputs(List<ModelInput> modelInputs) {
        this.modelInputs = modelInputs;
    }

    public SteadyComParams withModelInputs(List<ModelInput> modelInputs) {
        this.modelInputs = modelInputs;
        return this;
    }

    @JsonProperty("medium_upa")
    public String getMediumUpa() {
        return mediumUpa;
    }

    @JsonProperty("medium_upa")
    public void setMediumUpa(String mediumUpa) {
        this.mediumUpa = mediumUpa;
    }

    public SteadyComParams withMediumUpa(String mediumUpa) {
        this.mediumUpa = mediumUpa;
        return this;
    }

    @JsonProperty("flux_output")
    public String getFluxOutput() {
        return fluxOutput;
    }

    @JsonProperty("flux_output")
    public void setFluxOutput(String fluxOutput) {
        this.fluxOutput = fluxOutput;
    }

    public SteadyComParams withFluxOutput(String fluxOutput) {
        this.fluxOutput = fluxOutput;
        return this;
    }

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public SteadyComParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((("SteadyComParams"+" [modelInputs=")+ modelInputs)+", mediumUpa=")+ mediumUpa)+", fluxOutput=")+ fluxOutput)+", workspaceName=")+ workspaceName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
