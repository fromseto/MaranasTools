
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
 * <p>Original spec-file type: ModelInput</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "model_upa",
    "fixed_gr"
})
public class ModelInput {

    @JsonProperty("model_upa")
    private String modelUpa;
    @JsonProperty("fixed_gr")
    private Double fixedGr;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("model_upa")
    public String getModelUpa() {
        return modelUpa;
    }

    @JsonProperty("model_upa")
    public void setModelUpa(String modelUpa) {
        this.modelUpa = modelUpa;
    }

    public ModelInput withModelUpa(String modelUpa) {
        this.modelUpa = modelUpa;
        return this;
    }

    @JsonProperty("fixed_gr")
    public Double getFixedGr() {
        return fixedGr;
    }

    @JsonProperty("fixed_gr")
    public void setFixedGr(Double fixedGr) {
        this.fixedGr = fixedGr;
    }

    public ModelInput withFixedGr(Double fixedGr) {
        this.fixedGr = fixedGr;
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
        return ((((((("ModelInput"+" [modelUpa=")+ modelUpa)+", fixedGr=")+ fixedGr)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
