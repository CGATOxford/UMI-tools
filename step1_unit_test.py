import umi_tools.network as network
from itertools import product
input_data = {"ACGT": 456,
              "AAAT": 90,
              "ACAT": 72,
              "TCGT": 2,
              "CCGT": 2,
              "ACAG": 1}
output_data = {"unique": [['ACAG'], ['ACGT'], ['ACAT'], ['CCGT'], ['TCGT'], ['AAAT']],
               "percentile": [['ACAG'], ['ACGT'], ['ACAT'], ['CCGT'], ['TCGT'], ['AAAT']],
               "cluster": [['ACGT', 'AAAT', 'ACAT', 'CCGT', 'TCGT', 'ACAG']],
               "adjacency": [['ACGT', 'CCGT', 'TCGT'], ['AAAT'], ['ACAT', 'ACAG']],
               "directional": [['ACGT', 'ACAT', 'TCGT', 'CCGT', 'ACAG'], ['AAAT']]}


methods = ["unique", "percentile", "cluster", "adjacency", "directional"]

for method in methods:

    clusterer = network.UMIClusterer(method)
    clusters = clusterer(input_data.keys(), input_data, threshold=1)

    assert clusters == output_data[method], \
        "failed on method %s\n %s is not %s" % (method, clusters, output_data[method])
