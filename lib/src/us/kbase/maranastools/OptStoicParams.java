
package us.kbase.maranastools;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: OptStoicParams</p>
 * <pre>
 * model - the FBA model to use as a basis for modification
 * start_compound - the initial compound to be used as a source for the pathway
 * target_compound - the target compound to maximize yield for in the pathway
 * max_steps - the maximum number of steps to allow in the optimized pathway - any pathway
 *             created that has more than this number of steps is disqualified
 * use_heterologous_steps - allows adding
 * dG_threshold - a threshold free energy value to further constrain the path optimization
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "model",
    "start_compound",
    "target_compound",
    "max_steps",
    "use_heterologous_steps",
    "dG_threshold",
    "workspace_name"
})
public class OptStoicParams {

    @JsonProperty("model")
    private String model;
    @JsonProperty("start_compound")
    private String startCompound;
    @JsonProperty("target_compound")
    private String targetCompound;
    @JsonProperty("max_steps")
    private Long maxSteps;
    @JsonProperty("use_heterologous_steps")
    private Long useHeterologousSteps;
    @JsonProperty("dG_threshold")
    private Double dGThreshold;
    @JsonProperty("workspace_name")
    private String workspaceName;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("model")
    public String getModel() {
        return model;
    }

    @JsonProperty("model")
    public void setModel(String model) {
        this.model = model;
    }

    public OptStoicParams withModel(String model) {
        this.model = model;
        return this;
    }

    @JsonProperty("start_compound")
    public String getStartCompound() {
        return startCompound;
    }

    @JsonProperty("start_compound")
    public void setStartCompound(String startCompound) {
        this.startCompound = startCompound;
    }

    public OptStoicParams withStartCompound(String startCompound) {
        this.startCompound = startCompound;
        return this;
    }

    @JsonProperty("target_compound")
    public String getTargetCompound() {
        return targetCompound;
    }

    @JsonProperty("target_compound")
    public void setTargetCompound(String targetCompound) {
        this.targetCompound = targetCompound;
    }

    public OptStoicParams withTargetCompound(String targetCompound) {
        this.targetCompound = targetCompound;
        return this;
    }

    @JsonProperty("max_steps")
    public Long getMaxSteps() {
        return maxSteps;
    }

    @JsonProperty("max_steps")
    public void setMaxSteps(Long maxSteps) {
        this.maxSteps = maxSteps;
    }

    public OptStoicParams withMaxSteps(Long maxSteps) {
        this.maxSteps = maxSteps;
        return this;
    }

    @JsonProperty("use_heterologous_steps")
    public Long getUseHeterologousSteps() {
        return useHeterologousSteps;
    }

    @JsonProperty("use_heterologous_steps")
    public void setUseHeterologousSteps(Long useHeterologousSteps) {
        this.useHeterologousSteps = useHeterologousSteps;
    }

    public OptStoicParams withUseHeterologousSteps(Long useHeterologousSteps) {
        this.useHeterologousSteps = useHeterologousSteps;
        return this;
    }

    @JsonProperty("dG_threshold")
    public Double getDGThreshold() {
        return dGThreshold;
    }

    @JsonProperty("dG_threshold")
    public void setDGThreshold(Double dGThreshold) {
        this.dGThreshold = dGThreshold;
    }

    public OptStoicParams withDGThreshold(Double dGThreshold) {
        this.dGThreshold = dGThreshold;
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

    public OptStoicParams withWorkspaceName(String workspaceName) {
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
        return ((((((((((((((((("OptStoicParams"+" [model=")+ model)+", startCompound=")+ startCompound)+", targetCompound=")+ targetCompound)+", maxSteps=")+ maxSteps)+", useHeterologousSteps=")+ useHeterologousSteps)+", dGThreshold=")+ dGThreshold)+", workspaceName=")+ workspaceName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
