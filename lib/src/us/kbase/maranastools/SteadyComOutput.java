
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
 * <p>Original spec-file type: SteadyComOutput</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "report_name",
    "report_ref",
    "flux_output"
})
public class SteadyComOutput {

    @JsonProperty("report_name")
    private String reportName;
    @JsonProperty("report_ref")
    private String reportRef;
    @JsonProperty("flux_output")
    private String fluxOutput;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("report_name")
    public String getReportName() {
        return reportName;
    }

    @JsonProperty("report_name")
    public void setReportName(String reportName) {
        this.reportName = reportName;
    }

    public SteadyComOutput withReportName(String reportName) {
        this.reportName = reportName;
        return this;
    }

    @JsonProperty("report_ref")
    public String getReportRef() {
        return reportRef;
    }

    @JsonProperty("report_ref")
    public void setReportRef(String reportRef) {
        this.reportRef = reportRef;
    }

    public SteadyComOutput withReportRef(String reportRef) {
        this.reportRef = reportRef;
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

    public SteadyComOutput withFluxOutput(String fluxOutput) {
        this.fluxOutput = fluxOutput;
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
        return ((((((((("SteadyComOutput"+" [reportName=")+ reportName)+", reportRef=")+ reportRef)+", fluxOutput=")+ fluxOutput)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
