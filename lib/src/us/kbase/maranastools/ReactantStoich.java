
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
 * <p>Original spec-file type: Reactant_stoich</p>
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
    "start_compound",
    "fixed_stoich"
})
public class ReactantStoich {

    @JsonProperty("start_compound")
    private String startCompound;
    @JsonProperty("fixed_stoich")
    private Double fixedStoich;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("start_compound")
    public String getStartCompound() {
        return startCompound;
    }

    @JsonProperty("start_compound")
    public void setStartCompound(String startCompound) {
        this.startCompound = startCompound;
    }

    public ReactantStoich withStartCompound(String startCompound) {
        this.startCompound = startCompound;
        return this;
    }

    @JsonProperty("fixed_stoich")
    public Double getFixedStoich() {
        return fixedStoich;
    }

    @JsonProperty("fixed_stoich")
    public void setFixedStoich(Double fixedStoich) {
        this.fixedStoich = fixedStoich;
    }

    public ReactantStoich withFixedStoich(Double fixedStoich) {
        this.fixedStoich = fixedStoich;
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
        return ((((((("ReactantStoich"+" [startCompound=")+ startCompound)+", fixedStoich=")+ fixedStoich)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
