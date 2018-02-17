
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
 * <p>Original spec-file type: OptStoicParams</p>
 * <pre>
 * compound_id start_compound;
 * compound_id target_compound;
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "model",
    "reactant_stoichs",
    "product_stoichs",
    "integer_stoich",
    "objective",
    "exclude_compound_ids",
    "use_heterologous_steps",
    "num_pathways",
    "dG_threshold",
    "workspace_name"
})
public class OptStoicParams {

    @JsonProperty("model")
    private String model;
    @JsonProperty("reactant_stoichs")
    private List<ReactantStoich> reactantStoichs;
    @JsonProperty("product_stoichs")
    private List<ProductStoich> productStoichs;
    @JsonProperty("integer_stoich")
    private Long integerStoich;
    @JsonProperty("objective")
    private String objective;
    @JsonProperty("exclude_compound_ids")
    private String excludeCompoundIds;
    @JsonProperty("use_heterologous_steps")
    private Long useHeterologousSteps;
    @JsonProperty("num_pathways")
    private Long numPathways;
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

    @JsonProperty("reactant_stoichs")
    public List<ReactantStoich> getReactantStoichs() {
        return reactantStoichs;
    }

    @JsonProperty("reactant_stoichs")
    public void setReactantStoichs(List<ReactantStoich> reactantStoichs) {
        this.reactantStoichs = reactantStoichs;
    }

    public OptStoicParams withReactantStoichs(List<ReactantStoich> reactantStoichs) {
        this.reactantStoichs = reactantStoichs;
        return this;
    }

    @JsonProperty("product_stoichs")
    public List<ProductStoich> getProductStoichs() {
        return productStoichs;
    }

    @JsonProperty("product_stoichs")
    public void setProductStoichs(List<ProductStoich> productStoichs) {
        this.productStoichs = productStoichs;
    }

    public OptStoicParams withProductStoichs(List<ProductStoich> productStoichs) {
        this.productStoichs = productStoichs;
        return this;
    }

    @JsonProperty("integer_stoich")
    public Long getIntegerStoich() {
        return integerStoich;
    }

    @JsonProperty("integer_stoich")
    public void setIntegerStoich(Long integerStoich) {
        this.integerStoich = integerStoich;
    }

    public OptStoicParams withIntegerStoich(Long integerStoich) {
        this.integerStoich = integerStoich;
        return this;
    }

    @JsonProperty("objective")
    public String getObjective() {
        return objective;
    }

    @JsonProperty("objective")
    public void setObjective(String objective) {
        this.objective = objective;
    }

    public OptStoicParams withObjective(String objective) {
        this.objective = objective;
        return this;
    }

    @JsonProperty("exclude_compound_ids")
    public String getExcludeCompoundIds() {
        return excludeCompoundIds;
    }

    @JsonProperty("exclude_compound_ids")
    public void setExcludeCompoundIds(String excludeCompoundIds) {
        this.excludeCompoundIds = excludeCompoundIds;
    }

    public OptStoicParams withExcludeCompoundIds(String excludeCompoundIds) {
        this.excludeCompoundIds = excludeCompoundIds;
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

    @JsonProperty("num_pathways")
    public Long getNumPathways() {
        return numPathways;
    }

    @JsonProperty("num_pathways")
    public void setNumPathways(Long numPathways) {
        this.numPathways = numPathways;
    }

    public OptStoicParams withNumPathways(Long numPathways) {
        this.numPathways = numPathways;
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
        return ((((((((((((((((((((((("OptStoicParams"+" [model=")+ model)+", reactantStoichs=")+ reactantStoichs)+", productStoichs=")+ productStoichs)+", integerStoich=")+ integerStoich)+", objective=")+ objective)+", excludeCompoundIds=")+ excludeCompoundIds)+", useHeterologousSteps=")+ useHeterologousSteps)+", numPathways=")+ numPathways)+", dGThreshold=")+ dGThreshold)+", workspaceName=")+ workspaceName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
