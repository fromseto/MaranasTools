
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
 * <p>Original spec-file type: Product_stoich</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "target_compound",
    "fixed_stoich"
})
public class ProductStoich {

    @JsonProperty("target_compound")
    private String targetCompound;
    @JsonProperty("fixed_stoich")
    private Double fixedStoich;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("target_compound")
    public String getTargetCompound() {
        return targetCompound;
    }

    @JsonProperty("target_compound")
    public void setTargetCompound(String targetCompound) {
        this.targetCompound = targetCompound;
    }

    public ProductStoich withTargetCompound(String targetCompound) {
        this.targetCompound = targetCompound;
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

    public ProductStoich withFixedStoich(Double fixedStoich) {
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
        return ((((((("ProductStoich"+" [targetCompound=")+ targetCompound)+", fixedStoich=")+ fixedStoich)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
