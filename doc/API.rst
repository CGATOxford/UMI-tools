Making use of our Alogrithmns: The API
========================================

The UMIClusterer class
----------------------

There are occasions when you would like to make use of the deduplication/clustering algorithmns that we use here, but can't quite get the commandline tools to do exactly what you want. To help you use our methods in your tool we export the key class from the UMI-Tools package: `umi_tools.UMIClusterer`. It allows access to all the algorithmns detailed in the UMI-Tools paper.

To use the class, first create an instance, specifying the deduplication method. The default is the `directional` method, that is also the default for the command likne tools::

  from umi_tools import UMIClusterer
  clusterer = UMIClusterer(cluster_method="directional")

The returned instance is a *functor*. That is a class instance that can be called like a function. The function takes as its main argument a dictionary where the keys are UMI sequences and the value is the count of that UMI sequence. For example you might have reads with three different UMIs at one position: "ATAT", "GTAT", "CCAT". You may seen the first one 10 times, the second 5 times and the third 3 times::

  umis = {b"ATAT": 10,
          b"GTAT": 5,
	  b"CCAT": 3}

We pass this to our cluster, and it will return to us a list of lists, with each sub-list being a cluster of UMIs which we predict arose from the same original molecule::

  clustered_umis = clusterer(umis, threshold=1)
  print clustered_umis


  [["ATAT", "GTAT"], ["CCAT"]]

`ATAT` and `GTAT` are only a single base apart, and the count of
`ATAT` is >= 2n-1, where n is the count of `GTAT`. `CCAT` is two base
edits away from the others and so is in a seperate cluster. The order
of the UMIs in the groups is meaningful: we predict that `ATAT` was
the "true" UMI and so it is listed first.

Thus if you were in the bussiness of deduplicating reads, you'd keep one read associated with the `ATAT` and `CCAT` UMIs, and discard the reads associated with `GTAT`. Or if you were, for example, building a tools to geneate concensus read sequences, then you'd build one consensus from the reads with `ATAT` and `GTAT` and a different consensus read from the reads with the `CCAT` UMI.

Note that the `threshold` argument to the `UMIClusterer` allows the use of different edit distance threshold (although it is optional and the default is 1).

The `UMIClusterer` class is the only part of the API we publically export, and therefore guarentee won't change between major versions. You are welcome to peruse the rest of the API, especially if you are interested in contributing, but we can't guarentee that it will stay stable.
