'''
network.py - Network methods for dealing with UMIs
=========================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

'''
import collections
import itertools
import sys
import numpy as np

import pyximport
pyximport.install(build_in_temp=False)

try:
    from umi_tools._dedup_umi import edit_distance
    import umi_tools.Utilities as U

except:
    from _dedup_umi import edit_distance
    import Utilities as U

sys.setrecursionlimit(10000)


def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched


def recursive_search(node, adj_list):
    children = adj_list[node]
    children = [x for x in children if x not in recursive_search.component]
    for child in children:
        recursive_search.component.update((child,))
        recursive_search.component.update(
            recursive_search(child, adj_list))
    return recursive_search.component


def breadth_first_search_recursive(node, adj_list):
    try:
        recursive_search.component = set((node,))
        return recursive_search(node, adj_list)

    except RecursionError as error:
        U.info('Recursion Error: %s' % error)
        return breadth_first_search(node, adj_list)


def remove_umis(adj_list, cluster, nodes):
    '''removes the specified nodes from the cluster and returns
    the remaining nodes '''

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove


class UMIClusterer:
    '''A functor that clusters a dictionary of UMIs and their counts.
    The primary return value is either a list of representative UMIs
    or a list of lists where each inner list represents the contents of
    one cluster.

    Optionally:

      - identify the parent UMIs and return:
         - selected reads
         - umis
         - counts

    The initiation of the functor defines the methods:

      ** get_adj_list ** - returns the edges connecting the UMIs

      ** connected_components ** - returns clusters of connected components
                                   using the edges in the adjacency list

      ** get_best ** - returns the parent UMI(s) in the connected_components

      ** reduce_clusters ** - loops through the connected components in a
                              cluster and returns the unique reads. Optionally
                              returns lists of umis and counts per umi also

    Note: The get_adj_list and connected_components methods are not required by
    all custering methods. Where there are not required, the methods return
    None or the input parameters.

    '''

    # "get_best" methods #

    def _get_best_min_account(self, cluster, adj_list, counts):
        ''' return the min UMI(s) need to account for cluster'''
        if len(cluster) == 1:
            return list(cluster)

        sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                              reverse=True)

        for i in range(len(sorted_nodes) - 1):
            if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
                return sorted_nodes[:i+1]

    def _get_best_percentile(self, cluster, counts):
        ''' return all UMIs with counts >1% of the
        median counts in the cluster '''

        if len(cluster) == 1:
            return list(cluster)
        else:
            threshold = np.median(list(counts.values()))/100
            return [read for read in cluster if counts[read] > threshold]

    # "get_adj_list" methods #

    def _get_adj_list_adjacency(self, umis, counts, threshold):
        ''' identify all umis within hamming distance threshold'''

        adj_list = {umi: [] for umi in umis}
        for umi1, umi2 in itertools.combinations(umis, 2):
            if edit_distance(umi1, umi2) <= threshold:
                adj_list[umi1].append(umi2)
                adj_list[umi2].append(umi1)

        return adj_list

    def _get_adj_list_directional(self, umis, counts, threshold=1):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (2 * second umi counts)-1'''

        adj_list = {umi: [] for umi in umis}
        for umi1, umi2 in itertools.combinations(umis, 2):
            if edit_distance(umi1, umi2) <= threshold:
                if counts[umi1] >= (counts[umi2]*2)-1:
                    adj_list[umi1].append(umi2)
                if counts[umi2] >= (counts[umi1]*2)-1:
                    adj_list[umi2].append(umi1)

        return adj_list

    def _get_adj_list_null(self, umis, counts, threshold):
        ''' for methods which don't use a adjacency dictionary'''
        return None

    # "get_connected_components" methods #

    def _get_connected_components_adjacency(self, umis, graph, counts):
        ''' find the connected UMIs within an adjacency dictionary'''

        # TS: TO DO: Work out why recursive function doesn't lead to same
        # final output. Then uncomment below

        #if len(graph) < 10000:
        #    self.search = breadth_first_search_recursive
        #else:
        #    self.search = breadth_first_search

        found = set()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                #component = self.search(node, graph)
                component = breadth_first_search(node, graph)
                found.update(component)
                components.append(component)

        return components

    def _get_connected_components_null(self, umis, adj_list, counts):
        ''' for methods which don't use a adjacency dictionary'''
        return umis

    # "group" methods #

    def _group_unique(self, clusters, adj_list, counts):
        ''' return groups for unique method'''
        if len(clusters) == 1:
            groups = [clusters]
        else:
            groups = [[x] for x in clusters]

        return groups

    def _group_directional(self, clusters, adj_list, counts):
        ''' return groups for directional method'''

        observed = set()
        groups = []

        for cluster in clusters:
            if len(cluster) == 1:
                groups.append(list(cluster))
                observed.update(cluster)
            else:
                cluster = sorted(cluster, key=lambda x: counts[x],
                                 reverse=True)
                # need to remove any node which has already been observed
                temp_cluster = []
                for node in cluster:
                    if node not in observed:
                        temp_cluster.append(node)
                        observed.add(node)
                groups.append(temp_cluster)

        return groups

    def _group_adjacency(self, clusters, adj_list, counts):
        ''' return groups for adjacency method'''

        groups = []

        for cluster in clusters:
            if len(cluster) == 1:
                groups.append(list(cluster))

            else:
                observed = set()

                lead_umis = self._get_best_min_account(cluster,
                                                       adj_list, counts)
                observed.update(lead_umis)

                for lead_umi in lead_umis:
                    connected_nodes = set(adj_list[lead_umi])
                    groups.append([lead_umi] +
                                  list(connected_nodes - observed))
                    observed.update(connected_nodes)

        return groups

    def _group_cluster(self, clusters, adj_list, counts):
        ''' return groups for cluster or directional methods'''

        groups = []
        for cluster in clusters:
            groups.append(sorted(cluster, key=lambda x: counts[x],
                                 reverse=True))

        return groups

    def _group_percentile(self, clusters, adj_list, counts):
        ''' Return "groups" for the the percentile method. Note
        that grouping isn't really compatible with the percentile
        method. This just returns the retained UMIs in a structure similar
        to other methods '''

        retained_umis = self._get_best_percentile(clusters, counts)
        groups = [[x] for x in retained_umis]

        return groups

    def __init__(self, cluster_method="directional"):
        ''' select the required class methods for the cluster_method'''

        if cluster_method == "adjacency":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_adjacency

        elif cluster_method == "directional":
            self.get_adj_list = self._get_adj_list_directional
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_directional

        elif cluster_method == "cluster":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_cluster

        elif cluster_method == "percentile":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            # percentile method incompatible with defining UMI groups
            self.get_groups = self._group_percentile

        if cluster_method == "unique":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            self.get_groups = self._group_unique

    def __call__(self, umis, counts, threshold):
        '''Counts is a directionary that maps UMIs to their counts'''

        len_umis = [len(x) for x in umis]
        assert max(len_umis) == min(len_umis), (
            "not all umis are the same length(!):  %d - %d" % (
                min(len_umis), max(len_umis)))

        adj_list = self.get_adj_list(umis, counts, threshold)

        clusters = self.get_connected_components(umis, adj_list, counts)

        final_umis = [list(x) for x in
                      self.get_groups(clusters, adj_list, counts)]

        return final_umis


class ReadClusterer:
    '''This is a wrapper for applying the UMI methods to bundles of BAM reads.
    It is currently a pretty transparent wrapper on UMIClusterer. Basically
    taking a read bundle, extracting the UMIs and Counts, running UMIClusterer
    and returning the results along with annotated reads'''

    def __init__(self, cluster_method="directional"):

        self.UMIClusterer = UMIClusterer(cluster_method=cluster_method)

    def __call__(self, bundle, threshold, stats=False, further_stats=False,
                 deduplicate=True):
        '''Process the the bundled reads according to the method specified
        in the constructor. Note that in this implementation, stats and
        further_stats have no effect and their corresponding return values
        are always None. Only present to maintain signature with previous
        versions. Return signature is:
 
        reads, final_umis, umi_counts, topologies, nodes

        Meaning of these depends on whether deduplicate is True or not.
        If True:
        reads:        predicted best reads for deduplicated position
        final_umis:   list of predicted parent UMIs
        umi_counts:   Some of read counts for reads represented by the
                      corresponding UMI
        topologies:   Always None
        nodes:        Always None

        If False:
        reads:        identical to bundle as called
        final_umis:   list of lists of UMIs that are predicted to arise
                      from same biological molecule.
        umi_counts:   dictionary mapping UMIs to counts.
        topologies:   Always None
        nodes:        Alyways None

        '''

        umis = bundle.keys()
        counts = {umi: bundle[umi]["count"] for umi in umis}

        clusters = self.UMIClusterer(umis, counts, threshold)

        if deduplicate:
            final_umis = [cluster[0] for cluster in clusters]
            umi_counts = [sum(counts[umi] for umi in cluster)
                          for cluster in clusters]
            reads = [bundle[umi]["read"] for umi in final_umis]
        else:
            reads = bundle
            umi_counts = counts
            final_umis = clusters

        topologies = None
        nodes = None

        return (reads, final_umis, umi_counts, topologies, nodes)
