import pandas as pd
from snakemake.utils import validate


#### Config file and sample spreadsheet ####
configfile: "config/config.yaml"


# for testing: config = { "samples" : "config/samples.csv" }

samples = pd.read_csv(config["samples"]).set_index("Sample_ID", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")
