# Model metadata


<mark style="background-color: #FFE331">**Below is a template of the README.md file for the model-metadata folder of your hub. Italics in brackets are placeholders for information about your hub. **</mark>


This folder contains metadata files for the models submitting to the  *[hub name]*. The specification for these files has been adapted to be consistent with [model metadata guidelines in the hubverse documentation](https://docs.hubverse.io/en/latest/user-guide/model-metadata.html).

Each model is required to have metadata in 
[yaml format](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html).

These instructions provide detail about the [data
format](#Data-format) as well as [validation](#Data-validation) that
you can do prior to a pull request with a metadata file.

# Data format

## Required variables

This section describes each of the variables (keys) in the yaml document.
Please order the variables in this order.

### team_name
The name of your team that is less than 50 characters.

### team_abbr
The name of your team that is less than 16 characters.

### model_name
The name of your model that is less than 50 characters.

### model_abbr
An abbreviated name for your model that is less than 16 alphanumeric characters. 

### model_contributors

A list of all individuals involved in the forecasting effort.
A names, affiliations, and email address is required for each contributor. Individuals may also include an optional orcid identifiers.
All email addresses provided will be added to an email distribution list for model contributors.

The syntax of this field should be 
```
model_contributors: [
  {
    "name": "Modeler Name 1",
    "affiliation": "Institution Name 1",
    "email": "modeler1@example.com",
    "orcid": "1234-1234-1234-1234"
  },
  {
    "name": "Modeler Name 2",
    "affiliation": "Institution Name 2",
    "email": "modeler2@example.com",
    "orcid": "1234-1234-1234-1234"
  }
]
```

### license

One of the [accepted licenses](https://github.com/cdcepi/FluSight-forecast-hub/blob/673e983fee54f3a21448071ac46a9f78d27dd164/hub-config/model-metadata-schema.json#L69-L75).

We encourage teams to submit as a "cc-by-4.0" to allow the broadest possible uses
including private vaccine production (which would be excluded by the "cc-by-nc-4.0" license). 

### designated_model 

A team-specified boolean indicator (`true` or `false`) for whether the model should be considered eligible for inclusion in a Hub ensemble and public visualization. A team may specify up to two models as a designated_model for inclusion. Models which have a designated_model value of 'False' will still be included in internal forecasting hub evaluations.

### data_inputs

List or description of the data sources used to inform the model. Particularly those used beyond the target data of confirmed influenza hospital admissions.

### methods

A brief description of your forecasting methodology that is less than 200 
characters.

### methods_long

A full description of the methods used by this model. Among other details, this should include whether spatial correlation is considered and how the model accounts for uncertainty. If the model is modified, this field can also be used to provide the date of the modification and a description of the change.

### ensemble_of_models

A boolean value (`true` or `false`) that indicates whether a model is an ensemble of any separate component models.

### ensemble_of_hub_models

A boolean value (`true` or `false`) that indicates whether a model is an ensemble specifically of other models submitted to the forecasting hub.

## Optional

### model_version
An identifier of the version of the model

### website_url

A url to a website that has additional data about your model. 
We encourage teams to submit the most user-friendly version of your 
model, e.g. a dashboard, or similar, that displays your model forecasts. 

### repo_url

A github (or similar) repository url containing code for the model. 

### citation

One or more citations to manuscripts or preprints with additional model details. For example, "Gibson GC , Reich NG , Sheldon D. Real-time mechanistic bayesian forecasts of Covid-19 mortality. medRxiv. 2020. https://doi.org/10.1101/2020.12.22.20248736".

### team_funding 

Any information about funding source(s) for the team or members of the team that would be natural to include on any resulting publications. For example, "National Institutes of General Medical Sciences (R01GM123456). The content is solely the responsibility of the authors and does not necessarily represent the official views of NIGMS."

# Data validation

Optionally, you may validate a model metadata file locally before submitting it to the hub in a pull request. Note that this is not required, since the validations will also run on the pull request. To run the validations locally, follow these steps:

1. Create a fork of the `*[insert hub name]*` repository and then clone the fork to your computer.
2. Create a draft of the model metadata file for your model and place it in the `model-metadata` folder of this clone.
3. Install the hubValidations package for R by running the following command from within an R session:
``` r
remotes::install_github("hubverse-org/hubValidations")
```
4. Validate your draft metadata file by running the following command in an R session:
``` r
hubValidations::validate_model_metadata(
    hub_path="<path to your clone of the hub repository>",
    file_path="<name of your metadata file>")
```

For example, if your working directory is the root of the hub repository, you can use a command similar to the following:
``` r
hubValidations::validate_model_metadata(hub_path=".", file_path="UMass-trends_ensemble.yml")
```

If all is well, you should see output similar to the following:
```
✔ model-metadata-schema.json: File exists at path hub-config/model-metadata-schema.json.
✔ UMass-trends_ensemble.yml: File exists at path model-metadata/UMass-trends_ensemble.yml.
✔ UMass-trends_ensemble.yml: Metadata file extension is "yml" or "yaml".
✔ UMass-trends_ensemble.yml: Metadata file directory name matches "model-metadata".
✔ UMass-trends_ensemble.yml: Metadata file contents are consistent with schema specifications.
✔ UMass-trends_ensemble.yml: Metadata file name matches the `model_id` specified within the metadata file.
```

If there are any errors, you will see a message describing the problem.

<mark style="background-color: #FFE331">**As an example, here is a link to the [Flusight-Forecast_Hub model-metadata README](https://github.com/cdcepi/FluSight-forecast-hub/blob/master/model-metadata/README.md).**</mark>
