The network-based deduplication methods
=======================================

``dedup`` enables deduplication of UMIs by a variety of schemes. These are explained in the `UMI-tools publication <http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract>`_ (open access). A modified Figure 1e from the publication is reproduced below.

.. figure:: https://user-images.githubusercontent.com/6096414/93078070-bcfaf980-f681-11ea-87a0-db5ccb979875.png
   :width: 600
   :alt: UMI-tools figure 1e

   Schematic representation of UMI deduplication by the 6 methods.

All methods generate read groups which are inferred to represent duplicate reads from a single unique molecule prior to PCR. Where the read group contains more than one UMI, the UMI with the highest frequency is selected as the representative UMI for the group. Ties are broken randomly.

The simplest methods, *unique* and *percentile*, group reads with the exact same UMI. 

The network-based methods, *cluster*, *adjacency* and *directional*, build networks where nodes are UMIs and edges connect UMIs with an edit distance <= threshold (usually 1). The groups of reads are then defined from the network in a method-specific manner.

**cluster**: Form networks of connected UMIs (based on hamming distance threshold). Each connected component is a read group. In the above example, all the UMIs are contained in a single connected component and thus there is one read group containing all reads, with ACGT as the 'selected' UMI.

**adjacency**: Form networks as above. For each connected component, select the node (UMI) with the highest counts. Visit all nodes one edge away. If all nodes have been visited, stop. Otherwise, select the top 2 nodes with highest counts and visit all nodes one edge away. Repeat process until all nodes have been visited. Nodes which are not one of the n selected nodes are placed in a read group with the selected node with the highest counts for which there is an edge.

In the example above, ACGT (456) is selected first. Visiting all nodes one edge away leaves AAAT and ACAG unvisted. AAAT (90) is then additionally selected. Visiting all nodes one edge away from the selected nodes now leaves ACAG unvisited. Finally, ACAT (72) is additionally selected. Visiting all nodes one edge away now leaves no nodes unvisited. TCGT and CCGT are assigned to the ACGT read group, AAAT is in read group by itself and ACAG is assigned to the ACAT read group.

**directional** (default): Form networks with edges defined based on hamming distance threshold and node A counts >= (2 * node B counts) - 1. Each connected component is a read group, with the node with the highest counts selected as the top node for the component. In the example above, the directional edges yield two connected components. One with AAAT by itself and the other with the remaining UMIs with ACGT as the selected node.
